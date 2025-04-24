# Explanation of Aircraft Wing Material Distribution Optimization Code

This document explains the C++ code for optimizing the material distribution (thickness) of an aircraft wing using Ant Colony Optimization (ACO) and Finite Element Analysis (FEA). The system minimizes the wing’s weight while satisfying stress, displacement, and buckling constraints under multiple load cases (steady flight, gust, and maneuver). Below, we describe the code’s structure, its key components, and the mathematical equations and calculations used to develop the solution.

---

## System Overview

The code models an aircraft wing as a 2D thin structure, discretizing it into quadrilateral finite elements. Each element’s thickness is optimized to minimize the total weight, subject to structural constraints evaluated via FEA. The ACO algorithm searches for the optimal thickness distribution by simulating ant behavior, using pheromone trails and heuristics. The system includes:

- **WingStructure**: Defines the wing’s geometry (tapered, swept, dihedral) and mesh, managing nodes, elements, and global matrices.
- **QuadElement**: Represents a 4-node quadrilateral element, computing stiffness, mass, and stresses.
- **MaterialProperties**: Manages material data (e.g., Aluminum 7075-T6, Carbon Fiber) with isotropic/orthotropic properties.
- **FEAAnalyzer**: Performs FEA to compute stresses, displacements, and buckling factors.
- **LoadCase**: Generates load distributions for steady, gust, and maneuver conditions.
- **AntColonyOptimizer**: Optimizes thicknesses using ACO.
- **main**: Sets up the problem and runs the optimization.

The wing is modeled as a cantilevered structure fixed at the root, with loads applied in the z-direction (lift). The optimization balances weight minimization with structural integrity.

---

## Key Components and Their Equations

### 1. Wing Geometry and Mesh (WingStructure.cpp)

The wing is defined by its span, root chord, tip chord, sweep angle, and dihedral angle. The mesh is generated as a grid of quadrilateral elements (`elementsAlongSpan × elementsAlongChord`).

- **Node Coordinates**:
  - Nodes are placed based on spanwise position \( \eta = y / \text{span} \) and chordwise position \( \xi = x / \text{chord} \).
  - Chord at position \( y \):  
    \[
    \text{chord}(y) = \text{rootChord} - (\text{rootChord} - \text{tipChord}) \cdot \frac{y}{\text{span}}
    \]
  - Sweep offset: \( x_{\text{offset}} = y \cdot \tan(\text{sweepAngle}) \).
  - Dihedral offset: \( z_{\text{offset}} = y \cdot \tan(\text{dihedralAngle}) \).
  - Airfoil shape (NACA profile):  
    \[
    z(\xi) = 0.12 \cdot \text{chord} \cdot (0.2969\sqrt{\xi} - 0.1260\xi - 0.3516\xi^2 + 0.2843\xi^3 - 0.1015\xi^4)
    \]

- **Element Formation**:
  - Each quadrilateral element connects four nodes, forming a 4-node quad with 8 degrees of freedom (DOFs) (2 per node: x, y displacements).

- **Weight Calculation**:
  - Element area is computed using the cross-product of diagonal vectors:  
    \[
    \text{Area} = 0.5 \cdot \left\| (\mathbf{n}_1 - \mathbf{n}_0) \times (\mathbf{n}_2 - \mathbf{n}_0) \right\| + 0.5 \cdot \left\| (\mathbf{n}_2 - \mathbf{n}_0) \times (\mathbf{n}_3 - \mathbf{n}_0) \right\|
    \]
  - Total weight:  
    \[
    W = \sum_{\text{elements}} (\text{Area} \cdot t \cdot \rho)
    \]
    where \( t \) is thickness, \( \rho \) is material density.

### 2. Finite Element Analysis (QuadElement.cpp, FEAAnalyzer.cpp)

FEA computes the wing’s response (displacements, stresses) under applied loads. Each `QuadElement` represents a 4-node quadrilateral element, and `FEAAnalyzer` solves the global system.

- **Element Stiffness Matrix**:
  - Shape functions for a quad element at natural coordinates \( (\xi, \eta) \):  
    \[
    N_1 = 0.25(1-\xi)(1-\eta), \quad N_2 = 0.25(1+\xi)(1-\eta), \quad N_3 = 0.25(1+\xi)(1+\eta), \quad N_4 = 0.25(1-\xi)(1+\eta)
    \]
  - Jacobian matrix \( \mathbf{J} \):  
    \[
    \mathbf{J} = \begin{bmatrix}
    \sum_i x_i \frac{\partial N_i}{\partial \xi} & \sum_i x_i \frac{\partial N_i}{\partial \eta} \\
    \sum_i y_i \frac{\partial N_i}{\partial \xi} & \sum_i y_i \frac{\partial N_i}{\partial \eta} \\
    0 & 0
    \end{bmatrix}, \quad \mathbf{J}(2,2) = 1
    \]
  - Strain-displacement matrix \( \mathbf{B} \):  
    \[
    \mathbf{B} = \begin{bmatrix}
    \frac{\partial N_i}{\partial x} & 0 & \cdots \\
    0 & \frac{\partial N_i}{\partial y} & \cdots \\
    \frac{\partial N_i}{\partial y} & \frac{\partial N_i}{\partial x} & \cdots
    \end{bmatrix}
    \]
    where \( \frac{\partial N_i}{\partial x} = \mathbf{J}^{-1} \cdot \frac{\partial N_i}{\partial \xi} \).
  - Constitutive matrix \( \mathbf{D} \) (isotropic, plane stress):  
    \[
    \mathbf{D} = \frac{E}{1-\nu^2} \begin{bmatrix}
    1 & \nu & 0 \\
    \nu & 1 & 0 \\
    0 & 0 & \frac{1-\nu}{2}
    \end{bmatrix}
    \]
    **Note**: The code incorrectly uses this for composites; it should use `material.stiffnessMatrix`.
  - Stiffness matrix \( \mathbf{K}_e \):  
    \[
    \mathbf{K}_e = \int_{-1}^{1} \int_{-1}^{1} \mathbf{B}^T \mathbf{D} \mathbf{B} \cdot t \cdot \det(\mathbf{J}) \, d\xi \, d\eta
    \]
    Computed using 2×2 Gaussian quadrature with points \( \xi, \eta = \pm 0.57735 \), weights \( w = 1.0 \).

- **Element Mass Matrix**:
  - Mass matrix \( \mathbf{M}_e \):  
    \[
    \mathbf{M}_e = \int_{-1}^{1} \int_{-1}^{1} \mathbf{N}^T \mathbf{N} \cdot t \cdot \rho \cdot \det(\mathbf{J}) \, d\xi \, d\eta
    \]
    where \( \mathbf{N} \) is the shape function matrix, integrated via Gaussian quadrature.

- **Global System**:
  - Global stiffness matrix \( \mathbf{K} \): Assembled by summing element contributions at corresponding DOFs.
  - Equilibrium equation:  
    \[
    \mathbf{K} \mathbf{U} = \mathbf{F}
    \]
    where \( \mathbf{U} \) is the displacement vector, \( \mathbf{F} \) is the load vector.
  - Boundary conditions: Fixed DOFs are enforced by setting \( K_{ii} = 1 \), \( K_{ij} = 0 \), \( F_i = 0 \).
  - Solver: `Eigen::ConjugateGradient` with incomplete Cholesky preconditioner solves the sparse system.

- **Stress Calculation**:
  - Element stresses:  
    \[
    \boldsymbol{\sigma} = \mathbf{D} \mathbf{B} \mathbf{U}_e
    \]
    where \( \mathbf{U}_e \) is the element displacement vector (8 DOFs).
  - Von Mises stress (for metals):  
    \[
    \sigma_{\text{VM}} = \sqrt{\sigma_x^2 + \sigma_y^2 - \sigma_x \sigma_y + 3 \tau_{xy}^2}
    \]
  - Tsai-Wu failure criterion (for composites):  
    \[
    FI = F_1 \sigma_1 + F_2 \sigma_2 + F_{11} \sigma_1^2 + F_{22} \sigma_2^2 + F_{66} \tau_{12}^2 + 2 F_{12} \sigma_1 \sigma_2
    \]
    where \( F_1, F_2, F_{11}, F_{22}, F_{66}, F_{12} \) are strength-based coefficients, and failure occurs if \( FI \cdot \text{safetyFactor} \geq 1 \).

- **Buckling Check**:
  - Critical buckling stress for a plate:  
    \[
    \sigma_{\text{cr}} = \frac{k \pi^2 E}{12 (1-\nu^2)} \left( \frac{t}{a} \right)^2
    \]
    where \( k = 4.0 \), \( a \) is the characteristic length (average edge length), \( E \) is Young’s modulus, \( \nu \) is Poisson’s ratio.
  - Buckling factor:  
    \[
    \text{Factor} = \frac{\sigma_{\text{cr}}}{\sigma_{\text{VM}}}
    \]
    Failure if \( \text{Factor} < \text{safetyFactor} = 1.5 \).

### 3. Material Properties (MaterialProperties.cpp)

Materials (e.g., Aluminum 7075-T6, Carbon Fiber AS4/3501-6) are defined with density, Young’s modulus, Poisson’s ratio, yield strength, and stiffness matrices.

- **Isotropic Stiffness Matrix** (metals):  
  \[
  \mathbf{C} = \frac{E}{(1+\nu)(1-2\nu)} \begin{bmatrix}
  1-\nu & \nu & \nu & 0 & 0 & 0 \\
  \nu & 1-\nu & \nu & 0 & 0 & 0 \\
  \nu & \nu & 1-\nu & 0 & 0 & 0 \\
  0 & 0 & 0 & \frac{1-2\nu}{2} & 0 & 0 \\
  0 & 0 & 0 & 0 & \frac{1-2\nu}{2} & 0 \\
  0 & 0 & 0 & 0 & 0 & \frac{1-2\nu}{2}
  \end{bmatrix}
  \]
  Compliance matrix: \( \mathbf{S} = \mathbf{C}^{-1} \).

- **Orthotropic Stiffness Matrix** (composites):  
  \[
  C_{11} = \frac{E_1}{1-\nu_{12}\nu_{21}}, \quad C_{22} = \frac{E_2}{1-\nu_{12}\nu_{21}}, \quad C_{12} = \frac{\nu_{12} E_2}{1-\nu_{12}\nu_{21}}, \quad C_{66} = G_{12}
  \]
  where \( E_1, E_2, G_{12}, \nu_{12} \) are fiber-direction properties, and \( \nu_{21} = \nu_{12} \frac{E_2}{E_1} \).

### 4. Load Cases (LoadCase.cpp)

Three load cases simulate realistic flight conditions:

- **Steady Load**:
  - Elliptic lift distribution:  
    \[
    L(y) = \frac{4 L_{\text{total}}}{\pi \cdot \text{span}} \sqrt{1 - \left( \frac{2y}{\text{span}} - 1 \right)^2}
    \]
    where \( L_{\text{total}} = 15000 \cdot \text{span} \) N.

- **Gust Load**:
  - Steady load scaled by a gust profile:  
    \[
    L_{\text{gust}}(y) = L(y) \cdot \left( 1 + 1.5 \cdot \frac{1 - \cos\left( \frac{2\pi y}{\text{span}} \right)}{2} \right)
    \]

- **Maneuver Load**:
  - Steady load scaled by 2.5g with root-weighted distribution:  
    \[
    L_{\text{maneuver}}(y) = L(y) \cdot 2.5 \cdot \left( 1 + 0.5 \cdot \left( 1 - \frac{y}{\text{span}} \right) \right)
    \]

Loads are applied as nodal forces in the z-direction (y-index in the load vector).

### 5. Ant Colony Optimization (AntColonyOptimizer.cpp)

ACO optimizes element thicknesses to minimize weight. Each ant represents a thickness distribution, and the algorithm iterates to find the best solution.

- **Solution Construction**:
  - Thicknesses range from 3mm to 15mm, discretized into 13 steps:  
    \[
    t_i = t_{\text{min}} + i \cdot \frac{t_{\text{max}} - t_{\text{min}}}{\text{discretization} - 1}
    \]
  - Probability of selecting thickness \( t_i \) for element \( e \):  
    \[
    P(e, i) = \frac{\tau(e, i)^{\alpha} \cdot \eta(t_i)^{\beta}}{\sum_j \tau(e, j)^{\alpha} \cdot \eta(t_j)^{\beta}}
    \]
    where \( \tau(e, i) \) is the pheromone, \( \eta(t_i) = \frac{1}{t_i} \) is the heuristic (favoring thinner elements), \( \alpha = 1.0 \), \( \beta = 2.0 \).
  - Exploitation (probability \( q_0 = 0.7 \)): Choose the thickness with the highest \( \tau \cdot \eta \).
  - Exploration: Use roulette wheel selection based on \( P(e, i) \).

- **Fitness Evaluation**:
  - Fitness = wing weight if feasible (FEA converges, no constraint violations).
  - Constraints:  
    - Von Mises stress \( \sigma_{\text{VM}} \cdot 1.5 < \sigma_{\text{yield}} \) (metals).
    - Tsai-Wu index \( FI \cdot 1.5 < 1 \) (composites).
    - Buckling factor \( \frac{\sigma_{\text{cr}}}{\sigma_{\text{VM}}} \geq 1.5 \).

- **Pheromone Update**:
  - Evaporation:  
    \[
    \tau(e, i) \leftarrow (1 - \rho) \cdot \tau(e, i), \quad \rho = 0.3
    \]
  - Deposit (for feasible ants):  
    \[
    \tau(e, i) \leftarrow \tau(e, i) + \frac{1}{\text{fitness}}
    \]
  - Elite deposit (best solution):  
    \[
    \tau(e, i) \leftarrow \tau(e, i) + \frac{2}{\text{bestFitness}}
    \]

---

## How the Solution is Developed

1. **Initialization** (`main.cpp`):
   - A wing is created (10m span, 2m root chord, 1m tip chord, 25° sweep, 5° dihedral) with a 20×10 element mesh (200 elements, 231 nodes, 462 DOFs).
   - Material: Aluminum 7075-T6 (\( E = 71.7 \, \text{GPa}, \sigma_{\text{yield}} = 503 \, \text{MPa}, \rho = 2810 \, \text{kg/m}^3 \)).
   - Boundary conditions: 10 root nodes fixed in x and y directions.
   - ACO parameters: 30 ants, 50 iterations, thickness range 3mm–15mm (13 steps).

2. **Mesh and Load Setup**:
   - `WingStructure::generateWingGeometry` creates nodes and elements, applying sweep, dihedral, and airfoil shape.
   - `LoadCase::createStandardLoadCases` generates three load vectors, distributed elliptically and modified for gusts/manuvers.

3. **ACO Optimization**:
   - **Iteration Loop** (50 iterations):
     - Each ant constructs a thickness distribution (`constructSolutions`).
     - `evaluateAnt` sets thicknesses, runs FEA (`FEAAnalyzer::analyze`), and computes fitness (weight if feasible).
     - Pheromones are updated to favor low-weight, feasible solutions.
   - The best solution (lowest weight satisfying constraints) is tracked.

4. **FEA Evaluation**:
   - For each ant and load case:
     - Assemble global stiffness matrix (summing element \( \mathbf{K}_e \)).
     - Solve \( \mathbf{K} \mathbf{U} = \mathbf{F} \) for displacements.
     - Compute stresses (\( \boldsymbol{\sigma} = \mathbf{D} \mathbf{B} \mathbf{U}_e \)) and check failure (von Mises or Tsai-Wu).
     - Compute buckling factors for each element.
   - Constraints are checked: stress and buckling must satisfy a safety factor of 1.5.

5. **Output**:
   - The best thickness distribution and weight are reported.
   - Optimization history (iteration vs. best weight) is saved to `optimization_history.csv`.

---

## Key Calculations in Context

- **Weight Minimization**:
  - The objective is to minimize:  
    \[
    W = \sum_{\text{elements}} (\text{Area} \cdot t \cdot \rho)
    \]
  - ACO searches for the thickness \( t \) per element that minimizes \( W \) while satisfying constraints.

- **Structural Integrity**:
  - Stress constraints ensure the material does not fail:  
    \[
    \sigma_{\text{VM}} \cdot 1.5 \leq \sigma_{\text{yield}} \quad \text{(metals)}, \quad FI \cdot 1.5 < 1 \quad \text{(composites)}
    \]
  - Buckling constraints prevent instability:  
    \[
    \frac{\sigma_{\text{cr}}}{\sigma_{\text{VM}}} \geq 1.5
    \]

- **Load Response**:
  - Displacements \( \mathbf{U} \) are computed to ensure deflections are within acceptable limits (implicitly checked via FEA convergence).
  - Loads are distributed to simulate real flight conditions, ensuring the wing can withstand steady, gust, and maneuver forces.

- **Optimization Search**:
  - ACO balances exploration (trying new thicknesses) and exploitation (favoring known good solutions) using pheromones and heuristics.
  - The heuristic \( \eta(t) = \frac{1}{t} \) drives the solution toward thinner elements, reducing weight, while pheromones guide the search toward feasible regions.

---

## Limitations and Considerations

- **Composite Materials**: The code incorrectly uses isotropic \( \mathbf{D} \) for composites in `QuadElement`. It should use `material.stiffnessMatrix` for orthotropic properties.
- **Buckling Model**: The buckling factor uses a simplified \( k = 4.0 \), which may not account for complex boundary conditions or geometries.
- **Constraint Handling**: Infeasible solutions are discarded without penalty, potentially slowing convergence. A penalty term could improve performance.
- **Performance**: Repeated FEA evaluations (1500 calls for 30 ants × 50 iterations) are computationally expensive. Parallelization or caching could reduce runtime.

---

## Conclusion

The code effectively combines FEA and ACO to optimize aircraft wing material distribution. It uses well-established equations for stiffness, stress, buckling, and load distribution, ensuring engineering realism. The ACO algorithm efficiently searches the thickness space, guided by pheromones and a weight-minimizing heuristic. By addressing the identified limitations (e.g., composite stiffness, constraint handling), the system can be made more robust and efficient, suitable for real-world aerospace applications.
