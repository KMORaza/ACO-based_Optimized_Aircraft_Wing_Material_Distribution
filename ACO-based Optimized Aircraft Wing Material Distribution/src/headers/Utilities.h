#ifndef UTILITIES_H
#define UTILITIES_H

#include <Eigen/Dense>
#include <vector>

namespace WingUtils {
    Eigen::VectorXd createLoadVector(const std::vector<Eigen::Vector3d>& nodes,
                                    double loadMagnitude, int direction);
    void exportToCSV(const std::vector<double>& data, const std::string& filename);
}

#endif // UTILITIES_H
