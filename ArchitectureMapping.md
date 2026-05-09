# MATLAB to C# FEA Mapping

## Core Logic (FeaSolver.Core)

| MATLAB Function | C# Class/Method | Responsibility |
|-----------------|-----------------|----------------|
| `feisot4.m` | `Elements.Tet4.GetShapeFunctions` | Linear Tet shape functions & derivatives |
| `feisot10.m` | `Elements.Tet10.GetShapeFunctions` | Quadratic Tet shape functions & derivatives |
| `fejacob3.m` | `Math.Jacobian.Compute` | 3D Jacobian matrix and determinant |
| `fekine3d.m` | `Math.Kinematics.ComputeBMatrix` | B matrix (strain-displacement) |
| `feglqd3t.m` | `Math.GaussQuadrature.GetTetPoints` | Gauss points for tetrahedra |
| `feasmbl1.m` | `Solver.Assembler.AssembleStiffness` | Global stiffness matrix assembly |
| `feasmbl2.m` | `Solver.Assembler.AssembleForce` | Global force vector assembly |
| `feaplyc2.m` | `Solver.BoundaryConditions.ApplyDirichlet` | Applying fixed BCs (penalty or substitution) |
| `VonMisesCPA.m` | `Materials.VonMisesModel.UpdateStress` | Von Mises cutting plane algorithm |
| `TrescaReturnMap.m`| `Materials.TrescaModel.UpdateStress` | Tresca stress update |
| `eqStressCal.m` | `Math.TensorMath.CalculateEquivalentStress`| Mises/Tresca equivalent stress |
| `main_tet.m` | `Solver.NewtonRaphsonSolver` | Main NR loop, increment control |

## Data Structures

- `Node`: Id, X, Y, Z, Dofs
- `Element`: Id, NodeIndices, Material, IntegrationPoints
- `Mesh`: Nodes, Elements, BoundaryConditions
- `Material`: E, nu, YieldStress, HardeningParameters

## Desktop UI (FeaSolver.Desktop)

- `MainWindow`: Main container, Ribbon/Toolbar
- `Viewport3D`: HelixToolkit implementation for STL/Mesh display
- `SelectionHandler`: Raycasting for node/face selection
- `PostProcessor`: Nodal averaging and color mapping

## Math Utilities (FeaSolver.Core.Math)

| MATLAB Function | C# Class/Method | Responsibility |
|-----------------|-----------------|----------------|
| `mkmatrixS.m` | `TensorConverter.ToMatrix` | 6x1 stress vector to 3x3 matrix |
| `mkvectorS.m` | `TensorConverter.ToVector` | 3x3 stress matrix to 6x1 vector |
| `stressinv.m` | `TensorMath.GetInvariants` | Calculate I1, J2, J3 |
| `plfun.m` | `Hardening.HardeningFunction` | Hardening law (yield stress vs plastic strain) |

## Selection and BCs

Selection in HelixToolkit can be achieved using `Viewport3D.FindNearestPoint` or `FindHits`.
Once a face/node is selected, its ID is stored in a `BoundaryCondition` collection.
Loads (Force/Pressure) are mapped to nodal equivalent forces or element surface tractions.

## Post-Processing

Post-processing involves mapping results from Gauss points back to nodes.
`PostProcessor.RecoverNodalStresses` implements the recovery matrix logic from MATLAB's `main_tet.m`.
The HelixToolkit viewport uses `MeshGeometry3D` and `VertexColor` mapping to display stress distributions.
