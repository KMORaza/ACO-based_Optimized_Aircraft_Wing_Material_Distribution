#include "LoadCase.h"
#include <cmath>

std::vector<Eigen::VectorXd> LoadCase::createStandardLoadCases(const WingStructure& wing) {
    return {
        createSteadyLoad(wing),
        createGustLoad(wing),
        createManeuverLoad(wing)
    };
}
Eigen::VectorXd LoadCase::createSteadyLoad(const WingStructure& wing) {
    Eigen::VectorXd F(wing.getTotalNodes() * 2);
    F.setZero();
    double totalLift = 15000 * wing.getSpan(); // N per unit span
    for (int i = 0; i < wing.getTotalNodes(); ++i) {
        const auto& node = wing.getNode(i);
        double y = node.y();
        double eta = 2*y/wing.getSpan() - 1; // -1 to 1
        //// Elliptic lift distribution
        double lift = totalLift * 4/(M_PI*wing.getSpan()) * sqrt(1 - eta*eta);
        //// Apply in local z-direction
        F(2*i+1) = lift;
    }
    return F;
}
Eigen::VectorXd LoadCase::createGustLoad(const WingStructure& wing) {
    Eigen::VectorXd F = createSteadyLoad(wing);
    double gustFactor = 1.5;
    for (int i = 0; i < wing.getTotalNodes(); ++i) {
        const auto& node = wing.getNode(i);
        double y = node.y();
        //// 1-cos gust profile
        double gust = 0.5 * (1 - cos(2*M_PI*y/wing.getSpan()));
        F(2*i+1) *= (1 + gustFactor * gust);
    }
    return F;
}
Eigen::VectorXd LoadCase::createManeuverLoad(const WingStructure& wing) {
    Eigen::VectorXd F = createSteadyLoad(wing);
    double maneuverFactor = 2.5; // 2.5g maneuver
    for (int i = 0; i < wing.getTotalNodes(); ++i) {
        const auto& node = wing.getNode(i);
        double y = node.y();
        //// Increased load at root
        double distribution = 1.0 + 0.5*(1 - y/wing.getSpan());
        F(2*i+1) *= maneuverFactor * distribution;
    }
    return F;
}
