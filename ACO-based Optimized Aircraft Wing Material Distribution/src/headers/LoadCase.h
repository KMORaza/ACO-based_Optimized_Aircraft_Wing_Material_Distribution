#ifndef LOAD_CASE_H
#define LOAD_CASE_H

#include <vector>
#include <Eigen/Dense>
#include "WingStructure.h"

class LoadCase {
public:
    static std::vector<Eigen::VectorXd> createStandardLoadCases(const WingStructure& wing);

private:
    static Eigen::VectorXd createSteadyLoad(const WingStructure& wing);
    static Eigen::VectorXd createGustLoad(const WingStructure& wing);
    static Eigen::VectorXd createManeuverLoad(const WingStructure& wing);
};

#endif // LOAD_CASE_H
