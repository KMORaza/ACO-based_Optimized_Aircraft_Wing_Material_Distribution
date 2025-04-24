#include "WingStructure.h"
#include <cmath>
#include <array>
#include <algorithm>

WingStructure::WingStructure(double span, double rootChord, double tipChord,
                           double sweepAngle, double dihedralAngle,
                           int elementsAlongSpan, int elementsAlongChord)
    : span(span), rootChord(rootChord), tipChord(tipChord),
      sweepAngle(sweepAngle), dihedralAngle(dihedralAngle),
      elementsAlongSpan(elementsAlongSpan), elementsAlongChord(elementsAlongChord) {
    generateWingGeometry();
}

void WingStructure::generateWingGeometry() {
    int totalNodes = (elementsAlongSpan + 1) * (elementsAlongChord + 1);
    nodes.reserve(totalNodes);
    double sweepTan = tan(sweepAngle * M_PI / 180.0);
    double dihedralTan = tan(dihedralAngle * M_PI / 180.0);
    for (int i = 0; i <= elementsAlongSpan; ++i) {
        double eta = double(i) / elementsAlongSpan;
        double y = eta * span;
        double chord = rootChord - (rootChord - tipChord) * eta;
        double xOffset = y * sweepTan;
        double zOffset = y * dihedralTan;
        for (int j = 0; j <= elementsAlongChord; ++j) {
            double xi = double(j) / elementsAlongChord;
            double x = xOffset + xi * chord;
            double z = 0.12 * chord * (0.2969*sqrt(xi) - 0.1260*xi - 0.3516*xi*xi + 0.2843*xi*xi*xi - 0.1015*xi*xi*xi*xi);
            nodes.emplace_back(x, y, zOffset + z);
        }
    }
    elements.reserve(elementsAlongSpan * elementsAlongChord);
    for (int i = 0; i < elementsAlongSpan; ++i) {
        for (int j = 0; j < elementsAlongChord; ++j) {
            int n1 = i * (elementsAlongChord + 1) + j;
            int n2 = n1 + 1;
            int n3 = (i + 1) * (elementsAlongChord + 1) + j + 1;
            int n4 = (i + 1) * (elementsAlongChord + 1) + j;

            std::array<Eigen::Vector3d, 4> elemNodes = {
                nodes[n1], nodes[n2], nodes[n3], nodes[n4]
            };
            elements.push_back(std::make_unique<QuadElement>(elemNodes, material, 0.01));
        }
    }
    assembleGlobalMatrices();
}

void WingStructure::assembleGlobalMatrices() {
    const size_t totalDOF = nodes.size() * 2;
    globalK = Eigen::MatrixXd::Zero(totalDOF, totalDOF);
    globalM = Eigen::MatrixXd::Zero(totalDOF, totalDOF);
    for (size_t e = 0; e < elements.size(); ++e) {
        const auto& Ke = elements[e]->getStiffnessMatrix();
        const auto& Me = elements[e]->getMassMatrix();
        const auto& elemNodes = elements[e]->getNodes();
        for (const auto& node : elemNodes) {
            if (std::find(nodes.begin(), nodes.end(), node) == nodes.end()) {
                throw std::runtime_error("Element node not found in global nodes");
            }
        }
        std::array<size_t, 8> dofs;
        for (size_t i = 0; i < 4; ++i) {
            const auto it = std::find(nodes.begin(), nodes.end(), elemNodes[i]);
            const size_t nodeIdx = std::distance(nodes.begin(), it);
            dofs[2*i] = 2*nodeIdx;
            dofs[2*i+1] = 2*nodeIdx + 1;
            if (dofs[2*i] >= totalDOF || dofs[2*i+1] >= totalDOF) {
                throw std::runtime_error("Invalid DOF index in element " + std::to_string(e));
            }
        }
        for (size_t i = 0; i < 8; ++i) {
            for (size_t j = 0; j < 8; ++j) {
                globalK(dofs[i], dofs[j]) += Ke(i, j);
                globalM(dofs[i], dofs[j]) += Me(i, j);
            }
        }
    }
}
void WingStructure::applyLoad(const Eigen::VectorXd& loadVector) {
    loads = loadVector;
}
void WingStructure::setDisplacements(const Eigen::VectorXd& newDisplacements) {
    displacements = newDisplacements;
}
void WingStructure::setMaterial(const Material& material) {
    this->material = material;
    assembleGlobalMatrices();
}
void WingStructure::applyBoundaryConditions(const std::vector<bool>& fixedDOFs) {
    boundaryConditions = fixedDOFs;
}
void WingStructure::setElementThickness(int index, double thickness) {
    if (index >= 0 && static_cast<size_t>(index) < elements.size()) {
        const auto& oldElem = elements[index];
        elements[index] = std::make_unique<QuadElement>(oldElem->getNodes(),
                                                      oldElem->getMaterial(),
                                                      thickness);
        assembleGlobalMatrices();
    }
}
double WingStructure::getChordAtLocation(double y) const {
    return rootChord - (rootChord - tipChord) * (y / span);
}
double WingStructure::computeWeight() const {
    double weight = 0.0;
    for (const auto& elem : elements) {
        const auto& nodes = elem->getNodes();
        double area = 0.5 * ((nodes[1] - nodes[0]).cross(nodes[2] - nodes[0])).norm() +
                     0.5 * ((nodes[2] - nodes[0]).cross(nodes[3] - nodes[0])).norm();
        weight += area * elem->getThickness() * elem->getMaterial().density;
    }
    return weight;
}
