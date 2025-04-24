#include <iostream>
#include <chrono>
#include "WingStructure.h"
#include "FEAAnalyzer.h"
#include "AntColonyOptimizer.h"
#include "LoadCase.h"
#include "MaterialProperties.h"

int main() {
    try {
        MaterialDatabase materials;

        // wing structure (10m span, 2m root chord, 1m tip chord, 25 degree sweep, 5 degree dihedral)
        WingStructure wing(10.0, 2.0, 1.0, 25.0, 5.0, 20, 10);
        wing.setMaterial(materials.getMaterial("Aluminum 7075-T6"));

        std::vector<bool> fixedDOFs(wing.getTotalNodes() * 2, false);
        for (int i = 0; i < 10; ++i) {
            fixedDOFs[2*i] = true;
            fixedDOFs[2*i+1] = true;
        }
        wing.applyBoundaryConditions(fixedDOFs);

        FEAAnalyzer analyzer(wing);
        AntColonyOptimizer::Parameters acoParams;
        acoParams.colonySize = 30;
        acoParams.maxIterations = 50;
        acoParams.evaporationRate = 0.3;
        acoParams.minThickness = 0.003;
        acoParams.maxThickness = 0.015;
        acoParams.thicknessDiscretization = 13;

        auto start = std::chrono::high_resolution_clock::now();
        AntColonyOptimizer optimizer(wing, analyzer, acoParams);
        optimizer.optimize();
        auto end = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> elapsed = end - start;
        std::cout << "Optimization complete in " << elapsed.count() << " seconds\n";
        std::cout << "Final wing weight: " << optimizer.getBestFitness() << " kg\n";
        optimizer.saveResults("optimization_history.csv");

    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << std::endl;
        return 1;
    }
    return 0;
}
