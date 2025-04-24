#include "QuadElement.h"
#include <cmath>
constexpr double QuadElement::gp[QuadElement::ngp];
constexpr double QuadElement::gw[QuadElement::ngp];
QuadElement::QuadElement(const std::array<Eigen::Vector3d, 4>& nodes,
                        const Material& material, double thickness)
    : nodeCoords(nodes), material(material), thickness(thickness) {
    computeStiffnessMatrix();
    computeMassMatrix();
}
void QuadElement::computeStiffnessMatrix() {
    Ke = Eigen::MatrixXd::Zero(8, 8);
    double E = material.youngsModulus;
    double v = material.poissonRatio;
    Eigen::Matrix3d D;
    D << 1, v, 0,
         v, 1, 0,
         0, 0, (1-v)/2;
    D *= E / (1 - v*v);
    for (int i = 0; i < ngp; ++i) {
        for (int j = 0; j < ngp; ++j) {
            Eigen::Matrix<double, 3, 8> B = computeBMatrix(gp[i], gp[j]);
            Eigen::Matrix3d J = jacobian(gp[i], gp[j]);
            double detJ = J.determinant();
            if (detJ <= 0) {
                throw std::runtime_error("Invalid Jacobian determinant");
            }
            Ke += B.transpose() * D * B * detJ * thickness * gw[i] * gw[j];
        }
    }
}
void QuadElement::computeMassMatrix() {
    Me.setZero(8, 8);
    double rho = material.density;
    for (int i = 0; i < ngp; ++i) {
        for (int j = 0; j < ngp; ++j) {
            double xi = gp[i], eta = gp[j];
            double w = gw[i] * gw[j];
            Eigen::Vector4d N = shapeFunctions(xi, eta);
            Eigen::Matrix3d J = jacobian(xi, eta);
            double detJ = J.determinant();
            Eigen::MatrixXd Nm(2,8);
            Nm.setZero();
            for (int k = 0; k < 4; ++k) {
                Nm(0,2*k) = N(k);
                Nm(1,2*k+1) = N(k);
            }
            Me += Nm.transpose() * Nm * detJ * thickness * rho * w;
        }
    }
}
void QuadElement::computeStress(const Eigen::VectorXd& U) {
    stresses.setZero();
    Eigen::Matrix3d D;
    double E = material.youngsModulus;
    double v = material.poissonRatio;
    D << 1, v, 0,
         v, 1, 0,
         0, 0, (1-v)/2;
    D *= E / (1 - v*v);
    Eigen::Vector3d integratedStress = Eigen::Vector3d::Zero();
    for (int i = 0; i < ngp; ++i) {
        for (int j = 0; j < ngp; ++j) {
            double xi = gp[i], eta = gp[j];
            Eigen::Matrix<double, 3, 8> B = computeBMatrix(xi, eta);
            integratedStress += D * B * U;
        }
    }
    stresses = integratedStress / (ngp * ngp);
}
double QuadElement::getVonMisesStress() const {
    double sxx = stresses(0);
    double syy = stresses(1);
    double sxy = stresses(2);
    return sqrt(sxx*sxx + syy*syy - sxx*syy + 3*sxy*sxy);
}
bool QuadElement::checkFailure(double safetyFactor) const {
    if (material.name.find("Carbon Fiber") != std::string::npos) {
        // Tsai-Wu failure criterion
        double Xt = material.yieldStrength;
        double Xc = 0.9 * Xt;
        double Yt = 40e6;
        double Yc = 120e6;
        double S = 70e6;
        double F1 = 1/Xt - 1/Xc;
        double F2 = 1/Yt - 1/Yc;
        double F11 = 1/(Xt*Xc);
        double F22 = 1/(Yt*Yc);
        double F66 = 1/(S*S);
        double F12 = -0.5 * sqrt(F11*F22);
        double s1 = stresses(0);
        double s2 = stresses(1);
        double t12 = stresses(2);
        double FI = F1*s1 + F2*s2 + F11*s1*s1 + F22*s2*s2 + F66*t12*t12 + 2*F12*s1*s2;
        return (FI * safetyFactor) >= 1.0;
    } else {
        //// Von Mises for metals
        return (getVonMisesStress() * safetyFactor) >= material.yieldStrength;
    }
}
double QuadElement::getCharacteristicLength() const {
    double l1 = (nodeCoords[1] - nodeCoords[0]).norm();
    double l2 = (nodeCoords[2] - nodeCoords[1]).norm();
    double l3 = (nodeCoords[3] - nodeCoords[2]).norm();
    double l4 = (nodeCoords[0] - nodeCoords[3]).norm();
    return 0.25 * (l1 + l2 + l3 + l4);
}
// Shape function and derivative
Eigen::Vector4d QuadElement::shapeFunctions(double xi, double eta) const {
    return Eigen::Vector4d {
        0.25*(1-xi)*(1-eta),
        0.25*(1+xi)*(1-eta),
        0.25*(1+xi)*(1+eta),
        0.25*(1-xi)*(1+eta)
    };
}
Eigen::Matrix<double, 4, 2> QuadElement::shapeFunctionDerivatives(double xi, double eta) const {
    Eigen::Matrix<double, 4, 2> dN;
    dN << -0.25*(1-eta), -0.25*(1-xi),
           0.25*(1-eta), -0.25*(1+xi),
           0.25*(1+eta),  0.25*(1+xi),
          -0.25*(1+eta),  0.25*(1-xi);
    return dN;
}
Eigen::Matrix3d QuadElement::jacobian(double xi, double eta) const {
    Eigen::Matrix<double, 4, 2> dN = shapeFunctionDerivatives(xi, eta);
    Eigen::Matrix3d J = Eigen::Matrix3d::Zero();
    for (int i = 0; i < 4; ++i) {
        J(0,0) += nodeCoords[i].x() * dN(i,0);
        J(0,1) += nodeCoords[i].x() * dN(i,1);
        J(1,0) += nodeCoords[i].y() * dN(i,0);
        J(1,1) += nodeCoords[i].y() * dN(i,1);
    }
    J(2,2) = 1.0;
    return J;
}
Eigen::Matrix<double, 3, 8> QuadElement::computeBMatrix(double xi, double eta) const {
    Eigen::Matrix<double, 3, 8> B;
    B.setZero();
    Eigen::Matrix<double, 4, 2> dN = shapeFunctionDerivatives(xi, eta);
    Eigen::Matrix3d J = jacobian(xi, eta);
    Eigen::Matrix3d invJ = J.inverse();
    for (int i = 0; i < 4; ++i) {
        double dNdx = invJ(0,0)*dN(i,0) + invJ(0,1)*dN(i,1);
        double dNdy = invJ(1,0)*dN(i,0) + invJ(1,1)*dN(i,1);
        B(0, 2*i) = dNdx;
        B(1, 2*i+1) = dNdy;
        B(2, 2*i) = dNdy;
        B(2, 2*i+1) = dNdx;
    }
    return B;
}
