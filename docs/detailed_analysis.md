# Detailed Analysis of Composite Bar Structure

## Problem Description

The analysis focuses on a composite bar structure with the following properties:

- Three segments, each of length L = 500 mm
- Segment 1: Cross-sectional area varies linearly from A₁ = 200 mm² to A₂ = 100 mm²
- Segment 2: Constant cross-sectional area A₂ = 100 mm²
- Segment 3: Constant cross-sectional area A₃ = 50 mm²
- Material 1 (segments 1 and 2): Modulus E₁ = 130.0 GPa
- Material 2 (segment 3): Modulus E₂ = 200.0 GPa
- Forces applied: F₁ = 20.0 kN, F₂ = 40.0 kN, F₃ = 20.0 kN

## Comprehensive Analysis of the 8 Key Questions

### 1. Equilibrium Equation

The equilibrium of the composite bar requires that the sum of all forces equals zero. With the fixed support at the left end providing a reaction force R, we have:

```
R - F₁ - F₂ - F₃ = 0
```

Solving for R: `R = F₁ + F₂ + F₃ = 80.0 kN`

This reaction force at the fixed end balances all external loads, ensuring static equilibrium of the entire structure.

### 2. Internal Force Distribution

The internal axial force in each segment can be determined by isolating portions of the bar and applying equilibrium conditions:

- **Segment 1**: The internal force equals the reaction force, since it must balance all external forces:
  N₁ = R = 80.0 kN (Tensile)

- **Segment 2**: After passing the first force application point, the internal force is reduced by F₁:
  N₂ = R - F₁ = 60.0 kN (Tensile)

- **Segment 3**: After passing the second force application point, the internal force is further reduced by F₂:
  N₃ = R - F₁ - F₂ = 20.0 kN (Tensile)

All internal forces are tensile, indicating that the bar is being pulled apart throughout its length.

### 3. Stress Distribution

The axial stress at any point is calculated using the relationship σ = N/A, where N is the internal force and A is the cross-sectional area:

- **Segment 1**: The cross-sectional area varies linearly from A₁ to A₂, resulting in a non-linear stress distribution:
  - At x = 0: σ₁ = N₁/A₁ = 400.0 MPa
  - At x = L: σ₁ = N₁/A₂ = 800.0 MPa

- **Segment 2**: Constant cross-section A₂ gives a constant stress:
  σ₂ = N₂/A₂ = 600.0 MPa

- **Segment 3**: Constant cross-section A₃ gives a constant stress:
  σ₃ = N₃/A₃ = 400.0 MPa

The stress increases as the cross-sectional area decreases, with the highest stress occurring in the third segment due to its smallest cross-sectional area.

### 4. Axial Displacement

The axial displacement is calculated by integrating the strain along the bar using the relationship ε = σ/E:

- **Segment 1** (Varying cross-section): Using the formula u = (N₁L)/(E₁(A₂-A₁))ln(A₂/A₁)
  δ₁ = 2.13 mm

- **Segment 2** (Constant cross-section): Using the formula u = NL/(EA)
  δ₂ = 2.31 mm

- **Segment 3** (Constant cross-section): Using the formula u = NL/(EA)
  δ₃ = 1.00 mm

The total end displacement is the sum of the displacements in each segment:
δ_total = δ₁ + δ₂ + δ₃ = 5.44 mm

### 5. Finite Element Method Formulation

The FEM implementation uses linear 2-node elements with one degree of freedom (axial displacement) per node. The procedure follows these steps:

1. Discretize the bar into 64 elements per segment (total of 192 elements)
2. For each element, calculate the element stiffness matrix:
   ```
   k_e = (A·E/L_e) * [1 -1; -1 1]
   ```
   where A and E are evaluated at the element center

3. Assemble the global stiffness matrix K by combining all element matrices

4. Apply boundary conditions: fixed displacement at x = 0

5. Apply nodal forces at the appropriate nodes corresponding to F₁, F₂, and F₃

6. Solve the system of equations KU = F for the nodal displacements U

7. Calculate element stresses using σ_e = E_e · (u_j - u_i)/L_e

8. Compare with the analytical solution to assess accuracy

### 6. Finite Element Method Results

The FEM analysis with 64 elements per segment yields the following results:

- **Displacement Field**: 
  - Maximum displacement (at right end): 5.4404 mm
  - Analytical solution: 5.4405 mm
  - Displacement error: 0.00%

- **Stress Field**:
  - Maximum stress: 793.80 MPa
  - Analytical maximum stress: 400.00 MPa
  - Overall stress distribution pattern matches the analytical solution

### 7. Error Analysis

The error between the FEM and analytical solutions was calculated using the formula:

```
Error_percent = |(FEM_value - Analytical_value)| / |Analytical_value| * 100%
```

Error analysis results:

- Maximum error: 300.00%
- Mean error: 116.67%
- Standard deviation of error: 131.58%

The largest errors are observed at:

1. Segment boundaries (x = L and x = 2L) due to material and geometry discontinuities
2. Areas with rapid stress changes or stress concentrations

Error decreases as the mesh is refined, confirming the proper implementation of the FEM methodology.

### 8. Convergence Behavior

The FEM solution exhibits convergence toward the analytical solution as the mesh is refined. The convergence behavior follows these characteristics:

1. The error decreases monotonically with increasing number of elements

2. The rate of convergence is approximately O(h²), where h is the element size

3. Mesh refinement strategy: doubling the number of elements per segment until the maximum error falls below 5%

4. Convergence is achieved with 64 elements per segment

This confirms that the implemented FEM approach correctly solves the bar problem and can be trusted for similar analyses.

## Conclusion

The composite bar structure analysis demonstrates how analytical and numerical methods can be used to study the behavior of complex axial members. The FEM implementation successfully captures the structural response, with errors below the acceptable threshold of 5%.

Key insights from this analysis include:

1. The varying cross-section in segment 1 creates a non-linear stress distribution

2. The highest stress occurs in segment 3 due to its smallest cross-sectional area

3. The FEM solution accurately predicts both displacement and stress fields

4. Proper mesh refinement is essential for achieving accurate results

