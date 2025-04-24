"""
Enhanced plotting script for bar structure analysis
This script creates improved visualizations with detailed annotations and legends
using the refined FEM methodology that demonstrates proper convergence behavior
"""

import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from matplotlib.ticker import MaxNLocator
from bar_analysis import BarAnalysis
from utils import plot_with_labels

def create_enhanced_plots(bar_analysis, save_dir='plots'):
    """Create enhanced plots with detailed annotations and better styling
    
    Args:
        bar_analysis: Instance of BarAnalysis class with solved analysis
        save_dir: Directory to save plots
    """
    # Ensure the save directory exists
    os.makedirs(save_dir, exist_ok=True)
    
    # Extract parameters
    A1, A2, A3 = bar_analysis.A1, bar_analysis.A2, bar_analysis.A3
    E1, E2 = bar_analysis.E1, bar_analysis.E2
    L = bar_analysis.L
    F1, F2, F3 = bar_analysis.F1, bar_analysis.F2, bar_analysis.F3
    
    # Run analytical and FEM solutions if not already done
    x_analytical, displacement_analytical, stress_analytical = bar_analysis.solve_analytical()
    x_fem, nodal_displacements, element_stresses, error = bar_analysis.solve_fem(bar_analysis.num_elements_per_segment)
    
    # Store FEM results for later use
    bar_analysis.fem_results = (x_fem, nodal_displacements, element_stresses, error)
    
    # Get analytical solution data in more detail
    x_analytical = np.linspace(0, 3*L, 1000)
    stress_analytical = np.zeros_like(x_analytical)
    displacement_analytical = np.zeros_like(x_analytical)
    
    # Calculate analytical values at each point using our refined approach
    # This implementation correctly accounts for force distribution and equilibrium
    for i, x in enumerate(x_analytical):
        # Get area at position x
        area = bar_analysis.get_area_at_x(x)
        
        # Determine which segment this point is in and apply the correct internal force
        # In our refined approach, we properly account for force equilibrium
        if x <= L:
            # First segment: Internal force from equilibrium = F1 + F2 + F3
            internal_force = F1 + F2 + F3
            stress_analytical[i] = internal_force / area
        elif x <= 2*L:
            # Second segment: Internal force from equilibrium = F2 + F3
            internal_force = F2 + F3
            stress_analytical[i] = internal_force / area
        else:
            # Third segment: Internal force from equilibrium = F3
            internal_force = F3
            stress_analytical[i] = internal_force / area
        
        # Displacement is calculated from the analytical solution using strain integration
        if i > 0:
            dx = x_analytical[i] - x_analytical[i-1]
            E = E1 if x <= 2*L else E2  # Different elastic modulus in different segments
            du = stress_analytical[i] * dx / E  # du = (σ/E)*dx
            displacement_analytical[i] = displacement_analytical[i-1] + du
    
    # Get FEM solution data
    x_fem, nodal_displacements, element_stresses, error = bar_analysis.fem_results
    x_element_centers = (x_fem[:-1] + x_fem[1:]) / 2
    percent_error = np.abs(error) * 100  # Convert to percentage
    
    # Generate additional FEM data for convergence study
    mesh_sizes = [2, 4, 8, 16]  # Elements per segment
    convergence_errors = []
    
    # Compute errors for different mesh densities
    for mesh_size in mesh_sizes:
        # Create a temporary analysis object with the current mesh size
        temp_bar = BarAnalysis(
            A1=A1, A2=A2, A3=A3,
            E1=E1, E2=E2,
            L=L,
            F1=F1, F2=F2, F3=F3,
            num_elements_per_segment=mesh_size
        )
        
        # Solve FEM
        _, _, temp_elem_stresses, temp_error = temp_bar.solve_fem(mesh_size)
        
        # Calculate max error percentage
        max_error_pct = np.max(np.abs(temp_error)) * 100
        convergence_errors.append(max_error_pct)
    
    # 1. DISPLACEMENT FIELD PLOT
    plt.figure(figsize=(12, 8))
    plt.plot(x_analytical, displacement_analytical, 'b-', linewidth=2.5, label='Analytical Solution')
    plt.plot(x_fem, nodal_displacements, 'ro-', markersize=6, label='FEM Solution (Nodes)')
    
    # Add segment boundaries
    for i, x_pos in enumerate([L, 2*L]):
        plt.axvline(x_pos, color='gray', linestyle='--', alpha=0.6)
        plt.text(x_pos+15, max(displacement_analytical)*0.15, f'Segment {i+1}/{i+2} boundary', 
                 rotation=90, verticalalignment='center')
    
    # Add segment annotations
    mid_segment_positions = [L/2, L*1.5, L*2.5]
    segment_properties = [
        f"A₁={A1:.0f}-{A2:.0f} mm², E₁={E1/1000:.0f} GPa",
        f"A₂={A2:.0f} mm², E₁={E1/1000:.0f} GPa",
        f"A₃={A3:.0f} mm², E₂={E2/1000:.0f} GPa"
    ]
    
    for pos, prop in zip(mid_segment_positions, segment_properties):
        plt.annotate(prop, xy=(pos, max(displacement_analytical)*0.75),
                     xytext=(pos, max(displacement_analytical)*0.85),
                     horizontalalignment='center',
                     bbox=dict(boxstyle="round,pad=0.3", fc="white", ec="gray", alpha=0.7))
    
    # Add end displacement annotation
    plt.annotate(f'End displacement\n{displacement_analytical[-1]:.3f} mm', 
                xy=(3*L, displacement_analytical[-1]),
                xytext=(3*L-100, displacement_analytical[-1]*0.7),
                arrowprops=dict(arrowstyle="->", connectionstyle="arc3,rad=.2"))
    
    # Add grid and styling
    plt.grid(True, linestyle='--', alpha=0.6)
    plt.legend(loc='upper left', frameon=True, framealpha=0.9, shadow=True, fontsize=12)
    
    # Save plot
    plot_with_labels(
        plt, 
        "Axial Displacement Field in Composite Bar",
        "Position along bar, x (mm)",
        "Axial displacement, u(x) (mm)",
        os.path.join(save_dir, "enhanced_displacement_field.png")
    )
    
    # 2. STRESS FIELD PLOT
    plt.figure(figsize=(12, 8))
    plt.plot(x_analytical, stress_analytical, 'b-', linewidth=2.5, label='Analytical Solution')
    plt.plot(x_element_centers, element_stresses, 'ro', markersize=6, label='FEM Solution (Elements)')
    
    # Add error indicators at each element
    for i, (x, stress_fem, stress_exact) in enumerate(zip(x_element_centers, element_stresses, bar_analysis.get_stress_at_x(x_element_centers))):
        error_pct = abs((stress_fem - stress_exact) / stress_exact) * 100 if stress_exact != 0 else 0
        if error_pct < 5:  # Only annotate points with significant error
            continue
        plt.annotate(f'{error_pct:.1f}%', 
                    xy=(x, stress_fem), 
                    xytext=(0, 10),
                    textcoords='offset points',
                    fontsize=8,
                    arrowprops=dict(arrowstyle="->", connectionstyle="arc3,rad=.2"),
                    bbox=dict(boxstyle="round,pad=0.2", fc="white", ec="gray", alpha=0.7))
    
    # Add segment boundaries
    for i, x_pos in enumerate([L, 2*L]):
        plt.axvline(x_pos, color='gray', linestyle='--', alpha=0.6)
        plt.text(x_pos+15, min(stress_analytical)*1.05, f'Segment {i+1}/{i+2} boundary', 
                 rotation=90, verticalalignment='center')
    
    # Add stress value annotations
    segments = [(0, L, 'σ₁'), (L, 2*L, 'σ₂'), (2*L, 3*L, 'σ₃')]
    for start, end, label in segments:
        # Find indices for this segment
        segment_indices = (x_analytical >= start) & (x_analytical <= end)
        segment_stress = stress_analytical[segment_indices]
        segment_x = x_analytical[segment_indices]
        
        # Add annotation at middle of segment
        mid_idx = len(segment_stress) // 2
        mid_x = segment_x[mid_idx]
        mid_stress = segment_stress[mid_idx]
        
        if start == 0:  # First segment has varying stress
            plt.annotate(f'{label}: {mid_stress:.1f} MPa (varies)', 
                        xy=(mid_x, mid_stress),
                        xytext=(mid_x, mid_stress + 100 if mid_idx > len(segment_stress)//3 else mid_stress - 100),
                        horizontalalignment='center',
                        arrowprops=dict(arrowstyle="->", connectionstyle="arc3,rad=0"))
        else:  # Other segments have constant stress
            plt.annotate(f'{label}: {mid_stress:.1f} MPa (constant)', 
                        xy=(mid_x, mid_stress),
                        xytext=(mid_x, mid_stress + 100),
                        horizontalalignment='center',
                        arrowprops=dict(arrowstyle="->", connectionstyle="arc3,rad=0"))
    
    # Add force arrows
    forces = [(L, F1), (2*L, F2), (3*L, F3)]
    for pos, force in forces:
        # Arrow pointing left for force
        plt.annotate('', xy=(pos-20, max(stress_analytical)*0.6), xytext=(pos, max(stress_analytical)*0.6),
                    arrowprops=dict(arrowstyle="<-", lw=2, color='red'))
        # Label for force value
        plt.annotate(f'F = {force/1000:.1f} kN', 
                    xy=(pos, max(stress_analytical)*0.6),
                    xytext=(pos, max(stress_analytical)*0.7),
                    horizontalalignment='center',
                    bbox=dict(boxstyle="round,pad=0.3", fc="lightyellow", ec="orange", alpha=0.8))
    
    # Add grid and styling
    plt.grid(True, linestyle='--', alpha=0.6)
    plt.legend(loc='upper right', frameon=True, framealpha=0.9, shadow=True, fontsize=12)
    
    # Save plot
    plot_with_labels(
        plt, 
        "Axial Stress Distribution in Composite Bar",
        "Position along bar, x (mm)",
        "Axial stress, σ(x) (MPa)",
        os.path.join(save_dir, "enhanced_stress_field.png")
    )
    
    # 3. ERROR DISTRIBUTION PLOT
    plt.figure(figsize=(12, 8))
    bars = plt.bar(x_element_centers, percent_error, width=np.mean(x_fem[1:]-x_fem[:-1]), 
                    align='center', alpha=0.7, color='skyblue', edgecolor='navy')
    plt.axhline(y=5, color='r', linestyle='--', linewidth=2, label='5% Error Threshold')
    
    # Add segment boundaries
    for i, x_pos in enumerate([L, 2*L]):
        plt.axvline(x_pos, color='gray', linestyle='--', alpha=0.6)
        plt.text(x_pos+15, max(percent_error)*0.5, f'Segment {i+1}/{i+2} boundary', 
                 rotation=90, verticalalignment='center')
    
    # Annotate maximum error
    max_error_idx = np.argmax(percent_error)
    max_error_x = x_element_centers[max_error_idx]
    max_error_value = percent_error[max_error_idx]
    
    plt.annotate(f'Max error: {max_error_value:.2f}%', 
                xy=(max_error_x, max_error_value),
                xytext=(max_error_x + 100, max_error_value + 1),
                arrowprops=dict(arrowstyle="->", connectionstyle="arc3,rad=-.2"))
    
    # Add error statistics
    avg_error = np.mean(percent_error)
    plt.annotate(f'Average error: {avg_error:.2f}%\nStd deviation: {np.std(percent_error):.2f}%', 
                xy=(0.02, 0.92), xycoords='axes fraction',
                bbox=dict(boxstyle="round,pad=0.3", fc="white", ec="gray", alpha=0.7))
    
    # Color bars based on error magnitude
    for i, bar in enumerate(bars):
        if percent_error[i] > 5:
            bar.set_color('red')
            bar.set_edgecolor('darkred')
    
    # Add grid and styling
    plt.grid(True, linestyle='--', alpha=0.6)
    plt.legend(loc='upper right', frameon=True, framealpha=0.9, shadow=True, fontsize=12)
    
    # Save plot
    plot_with_labels(
        plt, 
        "Relative Error Between Analytical and FEM Solutions",
        "Position along bar, x (mm)",
        "Relative error in stress (%)",
        os.path.join(save_dir, "enhanced_error.png")
    )
    
    # 4. AREA DISTRIBUTION PLOT
    plt.figure(figsize=(12, 8))
    x_plot = np.linspace(0, 3*L, 1000)
    area = np.array([bar_analysis.get_area_at_x(x) for x in x_plot])
    
    plt.plot(x_plot, area, 'g-', linewidth=3, label='Cross-sectional Area')
    
    # Add segment boundaries
    for i, x_pos in enumerate([L, 2*L]):
        plt.axvline(x_pos, color='gray', linestyle='--', alpha=0.6)
        plt.text(x_pos+15, min(area)*0.9, f'Segment {i+1}/{i+2} boundary', 
                 rotation=90, verticalalignment='center')
    
    # Annotate area values at key positions
    plt.annotate(f'A₁ = {A1:.0f} mm²', xy=(0, A1),
                xytext=(50, A1 + 10),
                arrowprops=dict(arrowstyle="->", connectionstyle="arc3,rad=-.2"))
                
    plt.annotate(f'A₂ = {A2:.0f} mm²', xy=(1.5*L, A2),
                xytext=(1.5*L, A2 + 20),
                horizontalalignment='center',
                arrowprops=dict(arrowstyle="->", connectionstyle="arc3,rad=0"))
                
    plt.annotate(f'A₃ = {A3:.0f} mm²', xy=(2.5*L, A3),
                xytext=(2.5*L, A3 + 20),
                horizontalalignment='center',
                arrowprops=dict(arrowstyle="->", connectionstyle="arc3,rad=0"))
    
    # Add description of the variable section
    plt.annotate('Linearly varying\ncross-section', xy=(L/2, (A1+A2)/2),
                xytext=(L/2, (A1+A2)/2 + 30),
                horizontalalignment='center',
                bbox=dict(boxstyle="round,pad=0.3", fc="white", ec="green", alpha=0.7))
                
    plt.annotate('Constant\ncross-section', xy=(1.5*L, A2),
                xytext=(1.5*L - 100, A2 - 40),
                arrowprops=dict(arrowstyle="->", connectionstyle="arc3,rad=.2"),
                bbox=dict(boxstyle="round,pad=0.3", fc="white", ec="green", alpha=0.7))
                
    plt.annotate('Constant\ncross-section', xy=(2.5*L, A3),
                xytext=(2.5*L + 100, A3 - 20),
                arrowprops=dict(arrowstyle="->", connectionstyle="arc3,rad=-.2"),
                bbox=dict(boxstyle="round,pad=0.3", fc="white", ec="green", alpha=0.7))
    
    # Add material indicators
    plt.annotate(f'Material 1: E₁ = {E1/1000:.0f} GPa', xy=(L, 0),
                xytext=(L/2, min(area)*0.7),
                horizontalalignment='center',
                bbox=dict(boxstyle="round,pad=0.3", fc="lightblue", ec="blue", alpha=0.7))
                
    plt.annotate(f'Material 2: E₂ = {E2/1000:.0f} GPa', xy=(2.5*L, 0),
                xytext=(2.5*L, min(area)*0.7),
                horizontalalignment='center',
                bbox=dict(boxstyle="round,pad=0.3", fc="lightblue", ec="blue", alpha=0.7))
    
    # Add grid and styling
    plt.grid(True, linestyle='--', alpha=0.6)
    plt.legend(loc='upper right', frameon=True, framealpha=0.9, shadow=True, fontsize=12)
    
    # Save plot
    plot_with_labels(
        plt, 
        "Cross-sectional Area Distribution Along Bar",
        "Position along bar, x (mm)",
        "Cross-sectional area, A(x) (mm²)",
        os.path.join(save_dir, "enhanced_area_distribution.png")
    )
    
    # 5. COMBINED VISUALIZATION
    plt.figure(figsize=(15, 10))
    
    # Create a schematic diagram of the bar structure
    plt.subplot(211)
    # Draw the bar
    plt.plot([0, 3*L], [0, 0], 'k-', linewidth=3)
    
    # Draw the varying cross-sections
    bar_height = 50
    # Segment 1 (varying area)
    x_seg1 = np.linspace(0, L, 100)
    width_seg1_top = np.linspace(A1/10, A2/10, 100)
    width_seg1_bottom = -width_seg1_top
    plt.fill_between(x_seg1, width_seg1_top, width_seg1_bottom, color='lightblue', alpha=0.7)
    
    # Segment 2 (constant area)
    x_seg2 = np.linspace(L, 2*L, 100)
    width_seg2 = A2/10
    plt.fill_between(x_seg2, width_seg2, -width_seg2, color='lightblue', alpha=0.7)
    
    # Segment 3 (constant area)
    x_seg3 = np.linspace(2*L, 3*L, 100)
    width_seg3 = A3/10
    plt.fill_between(x_seg3, width_seg3, -width_seg3, color='lightsteelblue', alpha=0.7)
    
    # Add force arrows
    arrow_length = bar_height * 3
    plt.arrow(L, 0, -arrow_length/3, 0, head_width=10, head_length=10, fc='r', ec='r', linewidth=2)
    plt.arrow(2*L, 0, -arrow_length/3, 0, head_width=10, head_length=10, fc='r', ec='r', linewidth=2)
    plt.arrow(3*L, 0, -arrow_length/3, 0, head_width=10, head_length=10, fc='r', ec='r', linewidth=2)
    
    # Add force labels
    plt.text(L-arrow_length/6, 20, f'F₁ = {F1/1000:.0f} kN', ha='center', fontsize=10, 
             bbox=dict(boxstyle="round,pad=0.3", fc="white", ec="red", alpha=0.7))
    plt.text(2*L-arrow_length/6, 20, f'F₂ = {F2/1000:.0f} kN', ha='center', fontsize=10,
             bbox=dict(boxstyle="round,pad=0.3", fc="white", ec="red", alpha=0.7))
    plt.text(3*L-arrow_length/6, 20, f'F₃ = {F3/1000:.0f} kN', ha='center', fontsize=10,
             bbox=dict(boxstyle="round,pad=0.3", fc="white", ec="red", alpha=0.7))
    
    # Add fixed support symbol
    plt.plot([-10, 0], [bar_height, bar_height], 'k-', linewidth=2)
    plt.plot([-10, 0], [-bar_height, -bar_height], 'k-', linewidth=2)
    for i in range(-bar_height, bar_height+1, 10):
        plt.plot([-10, 0], [i, 0], 'k-', linewidth=1)
    
    # Add segment labels
    plt.text(L/2, -bar_height*2, 'Segment 1\nVarying Area\nE₁', ha='center', fontsize=10,
             bbox=dict(boxstyle="round,pad=0.3", fc="white", alpha=0.9))
    plt.text(1.5*L, -bar_height*2, 'Segment 2\nConstant Area\nE₁', ha='center', fontsize=10,
             bbox=dict(boxstyle="round,pad=0.3", fc="white", alpha=0.9))
    plt.text(2.5*L, -bar_height*2, 'Segment 3\nConstant Area\nE₂', ha='center', fontsize=10,
             bbox=dict(boxstyle="round,pad=0.3", fc="white", alpha=0.9))
    
    plt.title('Composite Bar Structure - Problem Schematic', fontsize=14)
    plt.axis('equal')
    plt.axis('off')
    
    # Plot the stress distribution below the schematic
    plt.subplot(212)
    plt.plot(x_analytical, stress_analytical, 'b-', linewidth=2.5, label='Analytical Stress')
    plt.plot(x_analytical, displacement_analytical*100, 'g-', linewidth=2.5, label='Analytical Displacement (×100)')
    plt.plot(x_element_centers, element_stresses, 'ro', markersize=6, label='FEM Stress')
    plt.plot(x_fem, nodal_displacements*100, 'mo', markersize=4, label='FEM Displacement (×100)')
    
    # Add convergence information
    mesh_description = f"{bar_analysis.num_elements_per_segment} elements/segment ({bar_analysis.num_elements_per_segment*3} total)"
    error_description = f"Max error: {np.max(np.abs(error))*100:.2f}%, Avg: {np.mean(np.abs(error))*100:.2f}%"
    plt.annotate(f"Mesh: {mesh_description}\n{error_description}", 
                xy=(0.98, 0.05), xycoords='axes fraction',
                ha='right', va='bottom',
                bbox=dict(boxstyle="round,pad=0.3", fc="white", ec="gray", alpha=0.9))
    
    # Add segment boundaries
    for i, x_pos in enumerate([L, 2*L]):
        plt.axvline(x_pos, color='gray', linestyle='--', alpha=0.6)
    
    plt.grid(True, linestyle='--', alpha=0.6)
    plt.legend(loc='upper right', frameon=True, framealpha=0.9, shadow=True, fontsize=10)
    plt.xlabel('Position along bar, x (mm)', fontsize=12)
    plt.ylabel('Stress (MPa) / Scaled Displacement', fontsize=12)
    plt.title('Analytical and FEM Solutions Comparison', fontsize=14)
    
    # Save combined visualization
    plt.tight_layout()
    plt.savefig(os.path.join(save_dir, "enhanced_combined_visualization.png"), dpi=300, bbox_inches='tight')
    plt.close()
    
    # CONVERGENCE STUDY PLOT
    plt.figure(figsize=(10, 6))
    plt.plot(mesh_sizes, convergence_errors, 'bo-', linewidth=2, markersize=8)
    plt.xscale('log', base=2)  # Use log scale for mesh sizes
    plt.yscale('log')  # Log scale for errors
    plt.xlabel('Number of Elements per Segment', fontsize=12)
    plt.ylabel('Maximum Error (%)', fontsize=12)
    plt.title('Convergence Study: Error vs. Mesh Density', fontsize=14)
    plt.grid(True, which='both', linestyle='--', alpha=0.6)
    
    # Add theoretical convergence line for reference (slope -1)
    ref_x = np.array([mesh_sizes[0], mesh_sizes[-1]])
    ref_y = convergence_errors[0] * (ref_x[0] / ref_x)**1  # Order of convergence = 1
    plt.plot(ref_x, ref_y, 'r--', linewidth=1.5, label='Theoretical 1st Order Convergence')
    
    # Add data points with values
    for i, (x, y) in enumerate(zip(mesh_sizes, convergence_errors)):
        plt.annotate(f'{y:.2f}%', 
                    xy=(x, y), 
                    xytext=(5, 5),
                    textcoords='offset points',
                    fontsize=9,
                    bbox=dict(boxstyle="round,pad=0.3", fc="white", ec="gray", alpha=0.7))
    
    plt.legend(loc='best')
    plt.gca().xaxis.set_major_formatter(plt.FuncFormatter(lambda x, _: f'{int(x)}'))
    plt.tight_layout()
    plt.savefig(os.path.join(save_dir, "enhanced_convergence_study.png"), dpi=300, bbox_inches='tight')
    plt.close()
    
    print(f"Enhanced plots created and saved to {os.path.abspath(save_dir)} directory")

if __name__ == "__main__":
    # Create an instance of BarAnalysis with standard parameters
    bar = BarAnalysis(
        A1=200, A2=100, A3=50,
        E1=130000, E2=200000,
        L=500,
        F1=20000, F2=40000, F3=20000,
        num_elements_per_segment=4
    )
    
    # Solve both analytical and FEM solutions
    x_analytical, stress_analytical, displacement_analytical = bar.solve_analytical()
    x_fem, nodal_displacements, element_stresses, error = bar.solve_fem()
    
    # Store FEM results for use in enhanced plotting
    bar.fem_results = (x_fem, nodal_displacements, element_stresses, error)
    
    # Create enhanced plots
    create_enhanced_plots(bar)
