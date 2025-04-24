#ifndef MATERIAL_PROPERTIES_H
#define MATERIAL_PROPERTIES_H
#include <vector>
#include <string>
#include <map>
#include <Eigen/Dense>

struct Material {
    std::string name;
    double density;         // kg/m^3
    double youngsModulus;   // Pa
    double poissonRatio;
    double yieldStrength;   // Pa
    double thermalExpansion; // 1/K
    Eigen::Matrix3d stiffnessMatrix;
    Eigen::Matrix3d complianceMatrix;
    void calculateDerivedProperties();
};
class MaterialDatabase {
public:
    MaterialDatabase();
    const Material& getMaterial(const std::string& name) const;
    std::vector<std::string> getAvailableMaterials() const;
private:
    std::map<std::string, Material> materials;
    void initializeMetals();
    void initializeComposites();
};
#endif // MATERIAL_PROPERTIES_H
