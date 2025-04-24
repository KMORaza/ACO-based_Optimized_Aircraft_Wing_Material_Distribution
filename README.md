## Optimizing Aircraft Wing Material Distribution using Ant Colony Optimization (ACO)

### Overview
* Optimize the thickness distribution of an aircraft wing to minimize its weight while ensuring structural integrity under specified load conditions. The wing is made of Aluminum 7075-T6 and modeled as a mesh of quadrilateral elements.
* Use ant colony optimization (ACO) to determine the optimal thickness for each element, subject to constraints on stress and buckling, evaluated via finite element analysis (FEA). The wing experiences steady, gust, and maneuver loads.

### Code Components
|  Implementation |   |
|---|---|
| __`WingStructure`__ | Represents the wing’s geometry and mesh, managing nodes, elements, and global matrices|
| __`QuadElement`__ | Implements 4-node quadrilateral elements for FEA, computing stiffness, mass, and stresses.|
| __`MaterialProperties`__ | Manages material data (e.g., Aluminum 7075-T6, Carbon Fiber) with isotropic and orthotropic properties.|
| __`FEAAnalyzer`__ | Performs FEA to compute stresses, displacements, and buckling factors under various loads.|
| __`LoadCase`__ | Generates realistic load distributions for steady flight, gusts, and maneuvers.|
| __`AntColonyOptimizer`__ | Optimizes element thicknesses using ACO to minimize weight while satisfying constraints.|
| __`main`__ | Sets up the wing, materials, and optimization parameters, running the ACO and reporting results.|

### Structure 
*  __`WingStructure`__ 
   * Geometry Generation: The wing is modeled as a tapered, swept, dihedral structure with airfoil profile (z-coordinate computed using a polynomial). The node placement accounts for sweep and dihedral, and elements are correctly formed as quadrilaterals.
   * Matrix Assembly: The global stiffness (globalK) and mass (globalM) matrices are assembled by mapping element DOFs to global DOFs. Error checks for invalid node indices and DOF mappings ensure robustness.
   * Weight Computation: The weight is calculated as the sum of element volumes (area × thickness × density), which is correct for thin structures.
* __`QuadElement`__
   * Stiffness and Mass Matrices: Uses 2×2 Gaussian quadrature for numerical integration, which is appropriate for linear quadrilateral elements. The stiffness matrix uses a plane stress constitutive matrix (`D`), and the mass matrix accounts for material density.
   * Stress Calculation: Computes stresses using the strain-displacement matrix (`B`) and averages over Gauss points. The von Mises stress is correctly implemented for isotropic materials.
   * Failure Criteria: Supports von Mises for metals and Tsai-Wu for composites, which is appropriate for the material types. The Tsai-Wu implementation includes realistic strength parameters (e.g., Xt, Xc, Yt, Yc, S).
* __`MaterialProperties`__
   * Material Database: Initializes realistic materials (Aluminum 7075-T6, Titanium 6Al-4V, Carbon Fiber AS4/3501-6) with accurate properties (density, Young’s modulus, etc.).
   * Stiffness Matrix: Correctly computes isotropic stiffness/compliance matrices for metals and orthotropic matrices for composites. The Carbon Fiber matrix uses realistic E1, E2, G12, and v12 values.
* __`FEAAnalyzer`__
   * Linear Solver: Uses `Eigen::ConjugateGradient` with an incomplete Cholesky preconditioner, suitable for sparse, symmetric positive-definite systems. Boundary conditions are applied by zeroing rows/columns and setting diagonal terms to 1, which is standard.
   * Stress and Constraints: Computes element stresses using local displacements and checks failure (via `QuadElement::checkFailure`) and buckling (via `computeBucklingFactor`). The buckling factor uses a simplified plate buckling formula, which is reasonable for thin elements.
* __`LoadCase`__
   * Steady Load: Uses an elliptic lift distribution, which is physically realistic for wings.
   * Gust Load: Applies a 1-cos gust profile with a 1.5× factor, simulating vertical gusts.
   * Maneuver Load: Scales loads by 2.5g with a root-weighted distribution, simulating high-g maneuvers.
   * The load vectors are correctly sized (2 DOFs per node), and the distributions are applied in the z-direction (y-index in `F`). The elliptic distribution and gust/maneuver factors align with aerospace standards.
* __`AntColonyOptimizer`__
   * Initialization: Sets uniform initial pheromones and creates a colony of ants.
   * Solution Construction: Uses a hybrid strategy (`q0 = 0.7` for exploitation vs. exploration) to select thicknesses based on pheromone (`alpha`) and heuristic (`beta`) weights. The heuristic (`1/thickness`) favors thinner elements to minimize weight.
   * Evaluation: Applies thicknesses to the wing, runs FEA, and computes fitness as weight if feasible (converged and no constraint violations).
   * Pheromone Update: Evaporates pheromones (`evaporationRate = 0.3`) and deposits based on fitness (`1/weight`) for feasible ants, with an elite deposit for the best solution.
   * The ACO implementation follows the standard algorithm, with proper roulette wheel selection and pheromone updates. The fitness function (weight) and constraint handling are appropriate for the problem.
* __`main`__
   * Creates a wing (10m span, 2m root chord, 1m tip chord, 25° sweep, 5° dihedral) with Aluminum 7075-T6, fixes 10 nodes at the root, and runs ACO with reasonable parameters (30 ants, 50 iterations, 3mm–15mm thickness range).
   * Boundary conditions (fixing x and y DOFs at the root) are appropriate for a cantilevered wing.
 
 ---

 Check this [webpage](https://aco-optimized-aircraft-wing-design.netlify.app/) about the solution 
