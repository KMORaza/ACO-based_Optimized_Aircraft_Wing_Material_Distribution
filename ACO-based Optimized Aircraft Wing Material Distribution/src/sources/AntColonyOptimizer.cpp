#include "AntColonyOptimizer.h"
#include "LoadCase.h"
#include <fstream>
#include <algorithm>

AntColonyOptimizer::AntColonyOptimizer(WingStructure& wing, FEAAnalyzer& analyzer,
                                     const Parameters& params)
    : wing(wing), analyzer(analyzer), params(params),
      rng(std::random_device{}()), realDist(0.0, 1.0) {
    initialize();
}
template<typename T>
const T& clamp(const T& value, const T& low, const T& high) {
    return (value < low) ? low : (high < value) ? high : value;
}
void AntColonyOptimizer::initialize() {
    int numElements = wing.getTotalElements();
    pheromones.resize(numElements, std::vector<double>(params.thicknessDiscretization, params.initialPheromone));
    ants.resize(params.colonySize);
    bestFitness = std::numeric_limits<double>::max();
}
void AntColonyOptimizer::optimize() {
    for (int iter = 0; iter < params.maxIterations; ++iter) {
        constructSolutions();
        updatePheromones();
        iterationFitness.push_back(bestFitness);
        std::cout << "Iteration " << iter << ": Best weight = " << bestFitness << " kg\n";
    }
}
void AntColonyOptimizer::saveResults(const std::string& filename) const {
    std::ofstream out(filename);
    out << "Iteration,BestFitness\n";
    for (size_t i = 0; i < iterationFitness.size(); ++i) {
        out << i << "," << iterationFitness[i] << "\n";
    }
}
void AntColonyOptimizer::constructSolutions() {
    double thicknessStep = (params.maxThickness - params.minThickness) / (params.thicknessDiscretization - 1);
    for (auto& ant : ants) {
        ant.thicknesses.resize(wing.getTotalElements());
        ant.feasible = true;
        for (int e = 0; e < wing.getTotalElements(); ++e) {
            if (realDist(rng) < params.q0) {
                int bestIndex = std::distance(pheromones[e].begin(),
                    std::max_element(pheromones[e].begin(), pheromones[e].end()));
                ant.thicknesses[e] = params.minThickness + bestIndex * thicknessStep;
            } else {
                std::vector<double> probabilities(params.thicknessDiscretization);
                double sum = 0.0;

                for (int i = 0; i < params.thicknessDiscretization; ++i) {
                    double thickness = params.minThickness + i * thicknessStep;
                    probabilities[i] = pow(pheromones[e][i], params.alpha) *
                                      pow(heuristicValue(thickness), params.beta);
                    sum += probabilities[i];
                }
                double r = realDist(rng) * sum;
                double cumsum = 0.0;
                int selection = 0;

                for (; selection < params.thicknessDiscretization; ++selection) {
                    cumsum += probabilities[selection];
                    if (cumsum >= r) break;
                }
                ant.thicknesses[e] = params.minThickness + selection * thicknessStep;
            }
        }
        evaluateAnt(ant);
    }
}
void AntColonyOptimizer::evaluateAnt(Ant& ant) {
    for (int i = 0; i < static_cast<int>(wing.getTotalElements()); ++i) {
        wing.setElementThickness(i, ant.thicknesses[i]);
    }
    std::vector<Eigen::VectorXd> loadCases = LoadCase::createStandardLoadCases(wing);
    FEAResult result = analyzer.analyze(loadCases);
    ant.feasible = result.hasConverged && !result.violatesConstraints;
    ant.fitness = result.weight;
    if (ant.feasible && ant.fitness < bestFitness) {
        bestFitness = ant.fitness;
        bestSolution = ant.thicknesses;
    }
}
void AntColonyOptimizer::updatePheromones() {
    double thicknessStep = (params.maxThickness - params.minThickness) / (params.thicknessDiscretization - 1);
    for (auto& row : pheromones) {
        for (auto& val : row) {
            val *= (1.0 - params.evaporationRate);
        }
    }
    for (const auto& ant : ants) {
        if (ant.feasible) {
            double deposit = 1.0 / ant.fitness;
            for (int e = 0; e < wing.getTotalElements(); ++e) {
                int index = static_cast<int>((ant.thicknesses[e] - params.minThickness) / thicknessStep + 0.5);
                index = clamp(index, 0, params.thicknessDiscretization - 1);
                pheromones[e][index] += deposit;
            }
        }
    }
    if (bestFitness < std::numeric_limits<double>::max()) {
        double eliteDeposit = 2.0 / bestFitness;
        for (int e = 0; e < wing.getTotalElements(); ++e) {
            int index = static_cast<int>((bestSolution[e] - params.minThickness) / thicknessStep + 0.5);
            index = clamp(index, 0, params.thicknessDiscretization - 1);
            pheromones[e][index] += eliteDeposit;
        }
    }
}
double AntColonyOptimizer::heuristicValue(double thickness) const {
    return 1.0 / thickness;
}
