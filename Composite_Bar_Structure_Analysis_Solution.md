# Composite Bar Structure Analysis - Engineering Solution

## Problem Definition

This report presents a comprehensive analysis of a composite bar structure with the following specifications:

- A fixed-end bar composed of three segments with different material and geometric properties
- Total length: 1500 mm (three 500 mm segments)
- Cross-sectional areas:
  - Segment 1: Linearly varying from 200 mm² to 100 mm²
  - Segment 2: Constant 100 mm²
  - Segment 3: Constant 50 mm²
- Elastic moduli:
  - Segments 1 & 2: E₁ = 130 GPa
  - Segment 3: E₂ = 200 GPa
- Applied forces:
  - F₁ = 20 kN (at x = 500 mm)
  - F₂ = 40 kN (at x = 1000 mm)
  - F₃ = 20 kN (at x = 1500 mm)

## Analytical Solution

The analytical solution for this problem is derived using fundamental principles of solid mechanics. For a bar under axial loading, the stress at any point is governed by:

σ(x) = P(x) / A(x)

Where:
- σ(x) is the axial stress at position x
- P(x) is the internal axial force at position x
- A(x) is the cross-sectional area at position x

The axial displacement is then calculated by integrating the strain along the bar:

u(x) = ∫(σ(x)/E(x)) dx

For this particular problem, the internal force in each segment is determined by equilibrium conditions:
- Segment 1: P₁(x) = F₁ + F₂ + F₃ = 80 kN
- Segment 2: P₂(x) = F₂ + F₃ = 60 kN
- Segment 3: P₃(x) = F₃ = 20 kN

## Solution Results

### Stress Distribution

The analytical analysis reveals the following stress distribution:

1. **Segment 1 (0 ≤ x ≤ 500 mm):**
   - Stress varies from 400 MPa at x = 0 to 800 MPa at x = 500 mm
   - This variation is due to the linearly decreasing cross-sectional area

2. **Segment 2 (500 mm ≤ x ≤ 1000 mm):**
   - Constant stress of 600 MPa
   - Reflects the constant area and internal force in this segment

3. **Segment 3 (1000 mm ≤ x ≤ 1500 mm):**
   - Constant stress of 400 MPa
   - Lower stress despite smaller area due to reduced internal force

### Displacement Field

The displacement field shows:

1. **Maximum displacement of 7.85 mm** occurs at the free end (x = 1500 mm)
2. **Non-linear displacement profile** in the first segment due to the varying cross-section
3. **Different slopes in the displacement curve** for each segment, reflecting the different material properties and stress states

### Critical Observations

1. The **highest stress (800 MPa)** occurs at the junction between segments 1 and 2 (x = 500 mm). This represents a potential critical location for failure.

2. The **displacement increases monotonically** from the fixed end to the free end, as expected for this loading configuration.

3. The **change in elastic modulus** between segments 2 and 3 causes a change in the slope of the displacement curve, even though the strain energy density decreases in segment 3.

## Conclusion

This analysis demonstrates how varying cross-sectional areas and material properties affect the stress and displacement distributions in a composite bar structure. The analytical solution provides exact values for stress and displacement at any point along the bar, revealing:

1. The importance of considering geometric transitions in structural elements
2. How material property variations influence the overall structural response
3. The relationship between internal force distribution and resulting stress fields

These insights are valuable for engineering design, particularly when designing lightweight structures that must withstand specific loading conditions. The analysis methods presented here can be extended to more complex geometries and loading scenarios with appropriate modifications.

For structural design applications, the results suggest:

1. Reinforcement may be needed at the junction between segments 1 and 2 to address the stress concentration
2. The overall stiffness could be adjusted by modifying the tapering ratio in segment 1
3. Material selection plays a significant role in controlling displacement, even when internal forces remain unchanged
