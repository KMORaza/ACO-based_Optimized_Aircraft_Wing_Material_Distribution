Aircraft Wing Material Distribution Optimization
This C++ project optimizes the material distribution (thickness) of an aircraft wing to minimize weight while ensuring structural integrity under various load conditions (steady flight, gust, and maneuver). It uses Ant Colony Optimization (ACO) for the optimization process and Finite Element Analysis (FEA) to evaluate structural responses. The wingLicenseThis project is licensed under the MIT License. See the LICENSE file for details.
Overview
The program models an aircraft wing as a 2D thin structure, discretized into quadrilateral finite elements. Each element’s thickness is optimized to minimize weight, subject to stress, displacement, and buckling constraints. ACO searches for the optimal thickness distribution, guided by pheromone trails and a heuristic favoring thinner elements. FEA computes stresses, displacements, and buckling factors for realistic load cases.
Key Components

WingStructure: Defines wing geometry (tapered, swept, dihedral) and mesh.
QuadElement: Computes stiffness, mass, and stresses for 4-node quadrilateral elements.
MaterialProperties: Manages material data (e.g., Aluminum 7075-T6, Carbon Fiber).
FEAAnalyzer: Performs FEA to evaluate structural responses.
LoadCase: Generates load distributions (steady, gust, maneuver).
AntColonyOptimizer: Optimizes thicknesses using ACO.
main: Sets up the problem and runs the optimization.

Mathematical Foundations
Below are the key equations used in the code, with LaTeX rendered as images via codecogs and plain-text fallbacks for clarity.
1. Wing Geometry and Weight

Chord Length: = \text{rootChord} - (\text{rootChord} - \text{tipChord}) \cdot \frac{y}{\text{span}})(Plain text: chord(y) = rootChord - (rootChord - tipChord) * (y / span))

Weight:![Weight](https://latex.codecogs.com/png.latex?W = \sum_{\text{elements}} (\text{Area} \cdot t \cdot \rho))(Plain text: W = sum(Area * t * ρ) over all elements)where ( t ) is thickness, ( \rho ) is density.


2. Finite Element Analysis

Stiffness Matrix:![Stiffness](https://latex.codecogs.com/png.latex?\mathbf{K}_e = \int_{-1}^{1} \int_{-1}^{1} \mathbf{B}^T \mathbf{D} \mathbf{B} \cdot t \cdot \det(\mathbf{J}) , d\xi , d\eta)(Plain text: Ke = integral(B^T * D * B * t * det(J) dξ dη) from ξ=-1 to 1, η=-1 to 1)where ( \mathbf{B} ) is the strain-displacement matrix, ( \mathbf{D} ) is the constitutive matrix.

Von Mises Stress (metals):![Von Mises](https://latex.codecogs.com/png.latex?\sigma_{\text{VM}} = \sqrt{\sigma_x^2 + \sigma_y^2 - \sigma_x \sigma_y + 3 \tau_{xy}^2})(Plain text: σ_VM = sqrt(σ_x^2 + σ_y^2 - σ_x * σ_y + 3 * τ_xy^2))

Tsai-Wu Criterion (composites):![Tsai-Wu](https://latex.codecogs.com/png.latex?FI = F_1 \sigma_1 + F_2 \sigma_2 + F_{11} \sigma_1^2 + F_{22} \sigma_2^2 + F_{66} \tau_{12}^2 + 2 F_{12} \sigma_1 \sigma_2)(Plain text: FI = F1 * σ1 + F2 * σ2 + F11 * σ1^2 + F22 * σ2^2 + F66 * τ12^2 + 2 * F12 * σ1 * σ2)Failure if ( FI \cdot 1.5 \geq 1 ).

Buckling Factor:![Buckling](https://latex.codecogs.com/png.latex?\sigma_{\text{cr}} = \frac{k \pi^2 E}{12 (1-\nu^2)} \left( \frac{t}{a} \right)^2, \quad \text{Factor} = \frac{\sigma_{\text{cr}}}{\sigma_{\text{VM}}})(Plain text: σ_cr = (k * π^2 * E) / (12 * (1-ν^2)) * (t / a)^2, Factor = σ_cr / σ_VM)where ( k = 4.0 ), ( a ) is element length.


3. Load Cases

Steady Load (elliptic distribution): = \frac{4 L_{\text{total}}}{\pi \cdot \text{span}} \sqrt{1 - \left( \frac{2y}{\text{span}} - 1 \right)^2})(Plain text: L(y) = (4 * L_total) / (π * span) * sqrt(1 - ((2y / span) - 1)^2)))

4. Ant Colony Optimization

Thickness Probability:![Probability](https://latex.codecogs.com/png.latex?P(e, i) = \frac{\tau(e, i)^{\alpha} \cdot \eta(t_i)^{\beta}}{\sum_j \tau(e, j)^{\alpha} \cdot \eta(t_j)^{\beta}})(Plain text: P(e, i) = (τ(e, i)^α * η(t_i)^β) / sum(τ(e, j)^α * η(t_j)^β))where ( \tau ) is pheromone, ( \eta(t_i) = 1/t_i ), ( \alpha = 1.0 ), ( \beta = 2.0 ).

Pheromone Update:![Pheromone](https://latex.codecogs.com/png.latex?\tau(e, i) \leftarrow (1 - \rho) \cdot \tau(e, i) + \frac{1}{\text{fitness}})(Plain text: τ(e, i) = (1 - ρ) * τ(e, i) + 1/fitness)where ( \rho = 0.3 ), fitness is weight.


How It Works

Setup: A wing (10m span, 2m root chord, 1m tip chord, 25° sweep, 5° dihedral) is meshed into 20×10 elements using Aluminum 7075-T6. Root nodes are fixed.
Load Cases: Three loads (steady, gust, maneuver) are applied, with elliptic distributions.
ACO Loop: For 50 iterations, 30 ants construct thickness distributions (3mm–15mm, 13 steps). FEA evaluates each solution’s weight and constraints.
FEA: Solves ( \mathbf{K} \mathbf{U} = \mathbf{F} ), computes stresses and buckling factors, and checks constraints (safety factor = 1.5).
Output: The best thickness distribution and weight are saved to optimization_history.csv.

Usage

Dependencies: Requires Eigen library for linear algebra.
Build: Compile with a C++17-compliant compiler (e.g., g++).
Run: Execute main.cpp to optimize the wing and generate results.
Output: Check optimization_history.csv for iteration-wise best weights.

Limitations

Composite Bug: QuadElement uses isotropic stiffness for composites (fix by using material.stiffnessMatrix).
Performance: Repeated FEA calls are slow; parallelization could help.
Constraints: Infeasible solutions lack penalties, slowing convergence.

Contributing
Contributions are welcome! Please submit pull requests or open issues for bugs, improvements, or new features.
License
MIT License. See LICENSE for details.
