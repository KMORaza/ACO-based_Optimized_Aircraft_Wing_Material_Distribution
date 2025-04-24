#ifndef WING_STRUCTURE_H
#define WING_STRUCTURE_H
#include <vector>
#include <memory>
#include <Eigen/Dense>
#include "MaterialProperties.h"
#include "QuadElement.h"

class WingStructure {
public:
    WingStructure(double span, double rootChord, double tipChord,
                 double sweepAngle, double dihedralAngle,
                 int elementsAlongSpan, int elementsAlongChord);
    void generateWingGeometry();
    void setMaterial(const Material& material);
    void setElementThickness(int index, double thickness);
    void applyLoad(const Eigen::VectorXd& loadVector);
    void applyBoundaryConditions(const std::vector<bool>& fixedDOFs);
    void setDisplacements(const Eigen::VectorXd& newDisplacements);
    double computeWeight() const;
    double computeMaxStress() const;
    double computeMaxDisplacement() const;
    double getChordAtLocation(double y) const;
    const Eigen::MatrixXd& getGlobalStiffnessMatrix() const {
        return globalK;
        }
    const Eigen::VectorXd& getDisplacements() const {
        return displacements;
        }
    const Eigen::VectorXd& getLoads() const {
        return loads;
        }
    const std::vector<bool>& getBoundaryConditions() const {
        return boundaryConditions;
        }
    int getTotalElements() const {
        return elements.size();
        }
    int getTotalNodes() const {
        return nodes.size();
        }
    double getSpan() const {
        return span;
        }
    const Eigen::Vector3d& getNode(int index) const {
        return nodes[index];
        }
    const QuadElement& getElement(int index) const {
        return *elements[index];
        }
private:
    void assembleGlobalMatrices();
    double span, rootChord, tipChord, sweepAngle, dihedralAngle;
    int elementsAlongSpan, elementsAlongChord;
    std::vector<Eigen::Vector3d> nodes;
    std::vector<std::unique_ptr<QuadElement>> elements;
    Material material;
    Eigen::MatrixXd globalK;
    Eigen::MatrixXd globalM;
    Eigen::VectorXd displacements;
    Eigen::VectorXd loads;
    std::vector<bool> boundaryConditions;
};
#endif // WING_STRUCTURE_H
