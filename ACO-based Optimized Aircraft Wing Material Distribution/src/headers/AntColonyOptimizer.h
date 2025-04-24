#ifndef ANT_COLONY_OPTIMIZER_H
#define ANT_COLONY_OPTIMIZER_H
#include <vector>
#include <random>
#include <algorithm>
#include <iostream>
#include "WingStructure.h"
#include "FEAAnalyzer.h"

class AntColonyOptimizer {
public:
    struct Parameters {
        int colonySize = 30;
        int maxIterations = 100;
        double evaporationRate = 0.5;
        double alpha = 1.0; // Pheromone influence
        double beta = 2.0;  // Heuristic influence
        double q0 = 0.7;    // Exploitation probability
        double initialPheromone = 1.0;
        double minThickness = 0.002; // 2mm
        double maxThickness = 0.020; // 20mm
        int thicknessDiscretization = 10;
    };
    AntColonyOptimizer(WingStructure& wing, FEAAnalyzer& analyzer,
                      const Parameters& params);
    void optimize();
    void saveResults(const std::string& filename) const;
    double getBestFitness() const {
        return bestFitness;
        }
    const std::vector<double>& getBestSolution() const {
        return bestSolution;
        }
private:
    struct Ant {
        std::vector<double> thicknesses;
        double fitness;
        bool feasible;
    };
    void initialize();
    void constructSolutions();
    void evaluateAnt(Ant& ant);
    void updatePheromones();
    double heuristicValue(double thickness) const;
    WingStructure& wing;
    FEAAnalyzer& analyzer;
    Parameters params;
    std::vector<std::vector<double>> pheromones;
    std::vector<Ant> ants;
    std::vector<double> bestSolution;
    double bestFitness;
    std::mt19937 rng;
    std::uniform_real_distribution<double> realDist;
    std::vector<double> iterationFitness;
};

#endif // ANT_COLONY_OPTIMIZER_H
