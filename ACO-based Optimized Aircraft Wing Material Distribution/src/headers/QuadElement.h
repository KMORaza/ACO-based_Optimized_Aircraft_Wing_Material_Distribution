#ifndef QUAD_ELEMENT_H
#define QUAD_ELEMENT_H
#include <array>
#include <Eigen/Dense>
#include "MaterialProperties.h"

class QuadElement {
public:
    QuadElement(const std::array<Eigen::Vector3d, 4>& nodes,
               const Material& material, double thickness) {
        if (thickness <= 0) throw std::invalid_argument("Thickness must be positive");
        if (material.youngsModulus <= 0) throw std::invalid_argument("Invalid material");
        const double area = computeArea(nodes);
        if (area <= 1e-12) throw std::invalid_argument("Degenerate element");
        this->nodeCoords = nodes;
        this->material = material;
        this->thickness = thickness;
        computeStiffnessMatrix();
        computeMassMatrix();
    }
    void computeStiffnessMatrix();
    void computeMassMatrix();
    void computeStress(const Eigen::VectorXd& U);
    const Eigen::MatrixXd& getStiffnessMatrix() const { return Ke; }
    const Eigen::MatrixXd& getMassMatrix() const { return Me; }
    double getVonMisesStress() const;
    bool checkFailure(double safetyFactor) const;
    const std::array<Eigen::Vector3d, 4>& getNodes() const { return nodeCoords; }
    double getThickness() const { return thickness; }
    const Material& getMaterial() const { return material; }
    double getCharacteristicLength() const;

private:
    std::array<Eigen::Vector3d, 4> nodeCoords;
    Material material;
    double thickness;
    Eigen::MatrixXd Ke; // 8x8 stiffness matrix
    Eigen::MatrixXd Me; // 8x8 mass matrix
    Eigen::Vector3d stresses; // σx, σy, τxy
    //// Shape functions and derivatives
    Eigen::Vector4d shapeFunctions(double xi, double eta) const;
    Eigen::Matrix<double, 4, 2> shapeFunctionDerivatives(double xi, double eta) const;
    Eigen::Matrix3d jacobian(double xi, double eta) const;
    Eigen::Matrix<double, 3, 8> computeBMatrix(double xi, double eta) const;
    //// Integration
    static constexpr int ngp = 2;
    static constexpr double gp[ngp] = {-0.577350269, 0.577350269};
    static constexpr double gw[ngp] = {1.0, 1.0};
};
private:
    double computeArea(const std::array<Eigen::Vector3d, 4>& nodes) const {
        return 0.5 * ((nodes[1] - nodes[0]).cross(nodes[2] - nodes[0]).norm() +
               0.5 * ((nodes[2] - nodes[0]).cross(nodes[3] - nodes[0]).norm();
    }
#endif // QUAD_ELEMENT_H
