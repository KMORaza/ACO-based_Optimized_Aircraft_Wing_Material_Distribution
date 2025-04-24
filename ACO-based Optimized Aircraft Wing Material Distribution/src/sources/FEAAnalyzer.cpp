#include "FEAAnalyzer.h"
#include <Eigen/Sparse>
#include <Eigen/IterativeLinearSolvers>
#include <iostream>
#include <cmath>

FEAAnalyzer::FEAAnalyzer(WingStructure& wing) : wing(wing) {}
FEAResult FEAAnalyzer::analyze(const std::vector<Eigen::VectorXd>& loadCases) {
    currentResult.hasConverged = true;
    currentResult.violatesConstraints = false;
    currentResult.maxStress = 0.0;
    currentResult.maxDisplacement = 0.0;
    currentResult.weight = wing.computeWeight();
    currentResult.elementStresses.clear();
    currentResult.elementStresses.reserve(wing.getTotalElements());
    try {
        for (const auto& load : loadCases) {
            wing.applyLoad(load);
            solveLinearSystem();
            computeStresses();
            if (!checkConstraints()) {
                currentResult.violatesConstraints = true;
            }
            double currentDisp = wing.getDisplacements().lpNorm<Eigen::Infinity>();
            if (currentDisp > currentResult.maxDisplacement) {
                currentResult.maxDisplacement = currentDisp;
            }
        }
    } catch (const std::exception& e) {
        std::cerr << "FEA Error: " << e.what() << std::endl;
        currentResult.hasConverged = false;
        currentResult.violatesConstraints = true;
    }
    return currentResult;
}
void FEAAnalyzer::solveLinearSystem() {
    Eigen::SparseMatrix<double> K = wing.getGlobalStiffnessMatrix().sparseView();
    Eigen::VectorXd F = wing.getLoads();
    if (K.rows() == 0 || K.cols() == 0) {
        throw std::runtime_error("Empty stiffness matrix");
    }
    if (K.rows() != K.cols()) {
        throw std::runtime_error("Non-square stiffness matrix");
    }
    if (K.rows() != F.size()) {
        throw std::runtime_error("Matrix/vector size mismatch");
    }
    applyBoundaryConditions(K, F);
    Eigen::ConjugateGradient<Eigen::SparseMatrix<double>,
                           Eigen::Lower|Eigen::Upper,
                           Eigen::IncompleteCholesky<double>> solver;
    solver.compute(K);
    if (solver.info() != Eigen::Success) {
        std::cerr << "Matrix decomposition failed with error: "
                 << solver.info() << std::endl;
        throw std::runtime_error("Matrix decomposition failed");
    }
    Eigen::VectorXd U = solver.solve(F);
    wing.setDisplacements(U);
}
void FEAAnalyzer::applyBoundaryConditions(Eigen::SparseMatrix<double>& K,
                                        Eigen::VectorXd& F) const {
    const auto& bc = wing.getBoundaryConditions();
    int totalDOF = K.rows();
    for (size_t i = 0; i < bc.size(); ++i) {
        if (i >= totalDOF) {
            throw std::runtime_error("Boundary condition index out of range");
        }
        if (bc[i]) {
            for (int k = 0; k < K.outerSize(); ++k) {
                for (Eigen::SparseMatrix<double>::InnerIterator it(K, k); it; ++it) {
                    if (it.row() == i || it.col() == i) {
                        if (it.row() == it.col()) {
                            it.valueRef() = 1.0;
                        } else {
                            it.valueRef() = 0.0;
                        }
                    }
                }
            }
            F[i] = 0.0;
        }
    }
}
void FEAAnalyzer::computeStresses() {
    const Eigen::VectorXd& U = wing.getDisplacements();
    currentResult.elementStresses.resize(wing.getTotalElements());
    for (int e = 0; e < wing.getTotalElements(); ++e) {
        const auto& elem = wing.getElement(e);
        Eigen::VectorXd Ue(8);
        const auto& nodes = elem.getNodes();
        for (int i = 0; i < 4; ++i) {
            int nodeIdx = std::distance(&wing.getNode(0), &nodes[i]);
            Ue(2*i) = U(2*nodeIdx);
            Ue(2*i+1) = U(2*nodeIdx+1);
        }
        const_cast<QuadElement&>(elem).computeStress(Ue);
        double vmStress = elem.getVonMisesStress();
        currentResult.elementStresses[e] = vmStress;
        if (vmStress > currentResult.maxStress) {
            currentResult.maxStress = vmStress;
        }
    }
}
bool FEAAnalyzer::checkConstraints() const {
    const double SAFETY_FACTOR = 1.5;
    for (int e = 0; e < wing.getTotalElements(); ++e) {
        const auto& elem = wing.getElement(e);
        if (elem.checkFailure(SAFETY_FACTOR)) {
            return false;
        }
        if (computeBucklingFactor(e) < SAFETY_FACTOR) {
            return false;
        }
    }
    return true;
}
double FEAAnalyzer::computeBucklingFactor(int e) const {
    const auto& elem = wing.getElement(e);
    double t = elem.getThickness();
    double a = elem.getCharacteristicLength();
    double E = elem.getMaterial().youngsModulus;
    double v = elem.getMaterial().poissonRatio;
    double sigma_vm = currentResult.elementStresses[e];
    if (sigma_vm <= 0) return std::numeric_limits<double>::max();
    double k = 4.0;
    double sigma_cr = (k * M_PI*M_PI * E) / (12*(1-v*v)) * pow(t/a, 2);
    return sigma_cr / sigma_vm;
}


