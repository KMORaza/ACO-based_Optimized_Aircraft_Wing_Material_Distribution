#ifndef FEA_ANALYZER_H
#define FEA_ANALYZER_H
#include <vector>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include "WingStructure.h"

struct FEAResult {
    double maxStress;
    double maxDisplacement;
    double weight;
    bool hasConverged;
    bool violatesConstraints;
    std::vector<double> elementStresses;
};
class FEAAnalyzer {
public:
    FEAAnalyzer(WingStructure& wing);
    FEAResult analyze(const std::vector<Eigen::VectorXd>& loadCases);
    double computeBucklingFactor(int elementIndex) const;
private:
    WingStructure& wing;
    mutable FEAResult currentResult;
    void solveLinearSystem();
    void computeStresses();
    bool checkConstraints() const;
    void applyBoundaryConditions(Eigen::SparseMatrix<double>& K, Eigen::VectorXd& F) const;
};
#endif // FEA_ANALYZER_H
