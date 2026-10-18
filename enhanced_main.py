"""
Enhanced main script for bar structure analysis
"""
import numpy as np
import matplotlib.pyplot as plt
from bar_analysis import BarAnalysis
import time
import sys

# Define parameters from the problem statement
A1 = 200  # mm²
A2 = 100  # mm²
A3 = 50   # mm²
E1 = 130  # GPa
E2 = 200  # GPa
L = 500   # mm
F1 = 20   # kN
F2 = 40   # kN
F3 = 20   # kN

# Convert units to consistent system (N and mm)
E1 = E1 * 1000  # Convert from GPa to N/mm²
E2 = E2 * 1000  # Convert from GPa to N/mm²
F1 = F1 * 1000  # Convert from kN to N
F2 = F2 * 1000  # Convert from kN to N
F3 = F3 * 1000  # Convert from kN to N

def main():
    try:
        print("========== BAR STRUCTURE ANALYSIS ==========")
        print("Composite Bar with Variable Cross-Section")
        print("===========================================")
        print(f"Parameters:")
        print(f"A1 = {A1} mm², A2 = {A2} mm², A3 = {A3} mm²")
        print(f"E1 = {E1/1000} GPa, E2 = {E2/1000} GPa")
        print(f"L = {L} mm")
        print(f"F1 = {F1/1000} kN, F2 = {F2/1000} kN, F3 = {F3/1000} kN")
        print()

        # 1. ANALYTICAL SOLUTION
        print("1. ANALYTICAL (CLOSED-FORM) SOLUTION:")
        print("-" * 40)
        
        # Calculate total reaction at fixed end
        R = F1 + F2 + F3
        print(f"Total reaction at fixed end: R = {R/1000:.1f} kN")
        
        # Internal forces in each segment
        N1 = R
        N2 = R - F1
        N3 = R - F1 - F2
        
        print(f"Internal forces:")
        print(f"  Segment 1: N1 = {N1/1000:.1f} kN")
        print(f"  Segment 2: N2 = {N2/1000:.1f} kN")
        print(f"  Segment 3: N3 = {N3/1000:.1f} kN")
        
        # Stress calculations
        sigma1_start = N1 / A1
        sigma1_end = N1 / A2
        sigma2 = N2 / A2
        sigma3 = N3 / A3
        
        print(f"Stress distribution:")
        print(f"  Segment 1: σ₁ varies from {sigma1_start:.1f} MPa to {sigma1_end:.1f} MPa")
        print(f"  Segment 2: σ₂ = {sigma2:.1f} MPa (constant)")
        print(f"  Segment 3: σ₃ = {sigma3:.1f} MPa (constant)")
        
        # Displacement calculations
        # For segment 1 with linearly varying cross-section
        # The formula is: δ₁ = (N₁*L)/(E₁(A₂-A₁))*ln(A₂/A₁)
        delta1 = (N1 * L) / (E1 * (A2 - A1)) * np.log(A2 / A1)
        
        # For segments 2 and 3 with constant cross-section
        # The formula is: δᵢ = (Nᵢ*L)/(Eᵢ*Aᵢ)
        delta2 = (N2 * L) / (E1 * A2)
        delta3 = (N3 * L) / (E2 * A3)
        
        # Total end displacement
        total_displacement = delta1 + delta2 + delta3
        
        print(f"Axial displacements:")
        print(f"  Segment 1: δ₁ = {delta1:.2f} mm")
        print(f"  Segment 2: δ₂ = {delta2:.2f} mm")
        print(f"  Segment 3: δ₃ = {delta3:.2f} mm")
        print(f"  Total end displacement: {total_displacement:.2f} mm")
        print()

        # 2. FINITE ELEMENT SOLUTION
        print("2. FINITE ELEMENT ANALYSIS:")
        print("-" * 40)
        
        # Initial number of elements per segment
        num_elements_per_segment = 4
        
        # Create the bar analysis object
        print(f"Creating initial mesh with {num_elements_per_segment} elements per segment...")
        bar = BarAnalysis(A1, A2, A3, E1, E2, L, F1, F2, F3, num_elements_per_segment)
        
        # Calculate analytical solution for comparison
        start_time = time.time()
        x_analytical, stress_analytical, displacement_analytical = bar.solve_analytical()
        analytical_time = time.time() - start_time
        print(f"Analytical solution computed in {analytical_time:.2f} seconds")
        
        # Use finite element method with different element counts until error is below 5%
        max_error = float('inf')
        max_iterations = 5
        iterations = 0
        
        while max_error > 0.05 and iterations < max_iterations:
            iterations += 1
            start_iter_time = time.time()
            print(f"FEM Iteration {iterations} with {num_elements_per_segment * 3} total elements")
            
            # Solve using FEM
            x_fem, nodal_displacements, element_stresses, error = bar.solve_fem(num_elements_per_segment)
            max_error = np.max(np.abs(error))
            
            print(f"  Maximum error: {max_error * 100:.2f}%")
            print(f"  FEM end displacement: {nodal_displacements[-1]:.4f} mm")
            print(f"  Iteration completed in {time.time() - start_iter_time:.2f} seconds")
            
            if max_error > 0.05 and iterations < max_iterations:
                num_elements_per_segment *= 2
                bar.update_num_elements(num_elements_per_segment)
                print(f"  Refining mesh to {num_elements_per_segment} elements per segment...")
        
        print(f"\nFinal FEM results:")
        print(f"  Number of elements: {num_elements_per_segment * 3}")
        print(f"  Maximum error: {max_error * 100:.2f}%")
        print(f"  End displacement: {nodal_displacements[-1]:.4f} mm")
        print(f"  Analytical end displacement: {total_displacement:.4f} mm")
        print(f"  Displacement error: {abs(nodal_displacements[-1] - total_displacement)/total_displacement * 100:.2f}%")
        print()
        
        # 3. VISUALIZATION
        print("3. GENERATING VISUALIZATIONS:")
        print("-" * 40)
        
        # Plot results
        try:
            bar.plot_results(x_analytical, stress_analytical, displacement_analytical, 
                           x_fem, nodal_displacements, element_stresses, error)
            print("Plots generated and saved successfully.")
        except Exception as e:
            print(f"Error generating plots: {str(e)}")
        
        # 4. REPORT GENERATION
        print("\n4. GENERATING REPORT:")
        print("-" * 40)
        
        # Generate report
        try:
            report_file = bar.generate_report()
            print(f"Basic report generated: {report_file}")
        except Exception as e:
            print(f"Error generating report: {str(e)}")
            
        print("\nAnalysis complete.")
            
    except KeyboardInterrupt:
        print("\nExecution interrupted by user. Terminating gracefully.")
    except Exception as e:
        print(f"\nAn error occurred: {str(e)}")
        print("Program terminating gracefully.")
    finally:
        print("\nExecution completed.")

if __name__ == "__main__":
    main()
