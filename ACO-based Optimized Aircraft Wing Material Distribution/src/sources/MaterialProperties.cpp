#include "MaterialProperties.h"
#include <stdexcept>

void Material::calculateDerivedProperties() {
    double E = youngsModulus;
    double v = poissonRatio;
    //// Isotropic stiffness matrix
    stiffnessMatrix.setZero();
    stiffnessMatrix(0,0) = stiffnessMatrix(1,1) = stiffnessMatrix(2,2) = 1 - v;
    stiffnessMatrix(0,1) = stiffnessMatrix(0,2) = stiffnessMatrix(1,0) =
    stiffnessMatrix(1,2) = stiffnessMatrix(2,0) = stiffnessMatrix(2,1) = v;
    stiffnessMatrix(3,3) = stiffnessMatrix(4,4) = stiffnessMatrix(5,5) = (1 - 2*v)/2;
    stiffnessMatrix *= E / ((1 + v) * (1 - 2*v));
    complianceMatrix = stiffnessMatrix.inverse();
}
MaterialDatabase::MaterialDatabase() {
    initializeMetals();
    initializeComposites();
}
void MaterialDatabase::initializeMetals() {
    //// Aluminum 7075-T6
    Material aluminum;
    aluminum.name = "Aluminum 7075-T6";
    aluminum.density = 2810.0;
    aluminum.youngsModulus = 71.7e9;
    aluminum.poissonRatio = 0.33;
    aluminum.yieldStrength = 503e6;
    aluminum.thermalExpansion = 23.6e-6;
    aluminum.calculateDerivedProperties();
    materials[aluminum.name] = aluminum;
    //// Titanium 6Al-4V
    Material titanium;
    titanium.name = "Titanium 6Al-4V";
    titanium.density = 4430.0;
    titanium.youngsModulus = 113.8e9;
    titanium.poissonRatio = 0.342;
    titanium.yieldStrength = 880e6;
    titanium.thermalExpansion = 8.6e-6;
    titanium.calculateDerivedProperties();
    materials[titanium.name] = titanium;
}
void MaterialDatabase::initializeComposites() {
    //// Carbon Fiber Epoxy (AS4/3501-6)
    Material carbonFiber;
    carbonFiber.name = "Carbon Fiber AS4/3501-6";
    carbonFiber.density = 1600.0;
    carbonFiber.poissonRatio = 0.28;
    //// Orthotropic properties
    carbonFiber.youngsModulus = 0;
    double E1 = 126e9;
    double E2 = 11e9;
    double G12 = 6.6e9;
    double v12 = 0.28;
    double v21 = v12 * E2 / E1;
    carbonFiber.stiffnessMatrix.setZero();
    carbonFiber.stiffnessMatrix(0,0) = E1 / (1 - v12*v21);
    carbonFiber.stiffnessMatrix(1,1) = E2 / (1 - v12*v21);
    carbonFiber.stiffnessMatrix(0,1) = v12 * E2 / (1 - v12*v21);
    carbonFiber.stiffnessMatrix(1,0) = carbonFiber.stiffnessMatrix(0,1);
    carbonFiber.stiffnessMatrix(2,2) = G12;
    carbonFiber.complianceMatrix = carbonFiber.stiffnessMatrix.inverse();
    carbonFiber.yieldStrength = 600e6;
    materials[carbonFiber.name] = carbonFiber;
}
const Material& MaterialDatabase::getMaterial(const std::string& name) const {
    auto it = materials.find(name);
    if (it == materials.end()) {
        throw std::runtime_error("Material not found: " + name);
    }
    return it->second;
}
std::vector<std::string> MaterialDatabase::getAvailableMaterials() const {
    std::vector<std::string> names;
    for (const auto& pair : materials) {
        names.push_back(pair.first);
    }
    return names;
}
