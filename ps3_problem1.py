"""
FEA Solver for Problem Set 3 - Problem 1: Three 2-node Elements
"""
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
import os

# Output directory for figures
output_dir = "figures"
if not os.path.exists(output_dir):
    os.makedirs(output_dir)

def solve_problem1():
    """
    Problem 1: Three 2-node elements
    - Node positions: x = 0, 2, 4, 6
    - Displacements at nodes 1, 2, 3, 4 are: 1, 3, 4, 6
    """
    print("Solving Problem 1: Three 2-node Elements")
    
    # Problem parameters
    node_positions = np.array([0, 2, 4, 6])  # x-coordinates of nodes
    displacements = np.array([1, 3, 4, 6])   # Given displacements at nodes
    
    # 1. Derive shape functions for all three elements
    print("\n1. Shape Functions for 2-node elements:")
    print("N₁(ξ) = (1-ξ)/2")
    print("N₂(ξ) = (1+ξ)/2")
    
    # For reference, ξ goes from -1 to 1 in each element
    # Element 1: nodes 1 and 2 (x = 0 to 2)
    # Element 2: nodes 2 and 3 (x = 2 to 4)
    # Element 3: nodes 3 and 4 (x = 4 to 6)
    
    # Element lengths
    L1 = node_positions[1] - node_positions[0]  # 2
    L2 = node_positions[2] - node_positions[1]  # 2
    L3 = node_positions[3] - node_positions[2]  # 2
    
    # 2. Find displacements in each element at ξ = 0.5
    xi = 0.5  # Given ξ value
    
    # Shape function values at ξ = 0.5
    N1_at_xi = (1 - xi) / 2  # 0.25
    N2_at_xi = (1 + xi) / 2  # 0.75
    
    # Displacements at ξ = 0.5 in each element
    u_elem1_at_xi = N1_at_xi * displacements[0] + N2_at_xi * displacements[1]
    u_elem2_at_xi = N1_at_xi * displacements[1] + N2_at_xi * displacements[2]
    u_elem3_at_xi = N1_at_xi * displacements[2] + N2_at_xi * displacements[3]
    
    print(f"\n2. Displacements at ξ = 0.5:")
    print(f"Element 1: u(ξ=0.5) = {u_elem1_at_xi:.4f}")
    print(f"Element 2: u(ξ=0.5) = {u_elem2_at_xi:.4f}")
    print(f"Element 3: u(ξ=0.5) = {u_elem3_at_xi:.4f}")
    
    # 3. Find strains in each element at ξ = 0.5
    # For a 2-node element, strain is constant: ε = (u₂-u₁)/L
    # dN₁/dξ = -0.5, dN₂/dξ = 0.5 (derivatives of shape functions)
    # Jacobian for each element: J = L/2
    # dN₁/dx = dN₁/dξ * dξ/dx = -0.5 * 2/L = -1/L
    # dN₂/dx = dN₂/dξ * dξ/dx = 0.5 * 2/L = 1/L
    
    strain_elem1 = (displacements[1] - displacements[0]) / L1
    strain_elem2 = (displacements[2] - displacements[1]) / L2
    strain_elem3 = (displacements[3] - displacements[2]) / L3
    
    print(f"\n3. Strains at ξ = 0.5:")
    print(f"Element 1: ε(ξ=0.5) = {strain_elem1:.4f}")
    print(f"Element 2: ε(ξ=0.5) = {strain_elem2:.4f}")
    print(f"Element 3: ε(ξ=0.5) = {strain_elem3:.4f}")
    
    # 4. Calculate the stiffness matrix of all elements using Gauss integration
    # For simplicity, we'll assume E (Young's modulus) = 1 and A (Cross-sectional area) = 1
    E = 1.0  # Young's modulus (assumed)
    A = 1.0  # Cross-sectional area (assumed)
    
    # Stiffness matrix for each element: k = ∫(B^T·E·B·A·det(J))dξ
    # For 2-node elements, B = [-1/L, 1/L], det(J) = L/2
    # Gauss integration points and weights for 2-point Gauss quadrature
    gauss_points = np.array([-1/np.sqrt(3), 1/np.sqrt(3)])
    gauss_weights = np.array([1.0, 1.0])
    
    # Calculate stiffness matrices
    k1 = calculate_stiffness_matrix(L1, E, A, gauss_points, gauss_weights)
    k2 = calculate_stiffness_matrix(L2, E, A, gauss_points, gauss_weights)
    k3 = calculate_stiffness_matrix(L3, E, A, gauss_points, gauss_weights)
    
    print(f"\n4. Element Stiffness Matrices:")
    print(f"Element 1 Stiffness Matrix:")
    print(k1)
    print(f"Element 2 Stiffness Matrix:")
    print(k2)
    print(f"Element 3 Stiffness Matrix:")
    print(k3)
    
    # Generate plots
    generate_plots_problem1(node_positions, displacements, u_elem1_at_xi, u_elem2_at_xi, u_elem3_at_xi, 
                           strain_elem1, strain_elem2, strain_elem3, k1, k2, k3)
    
    return {
        'node_positions': node_positions,
        'displacements': displacements,
        'element_strains': [strain_elem1, strain_elem2, strain_elem3],
        'element_stiffness': [k1, k2, k3]
    }

def calculate_stiffness_matrix(L, E, A, gauss_points, gauss_weights):
    """Calculate stiffness matrix for a 2-node element using Gauss integration"""
    k = np.zeros((2, 2))
    
    # Jacobian: J = L/2
    J = L / 2
    
    for i, xi in enumerate(gauss_points):
        # B matrix at this Gauss point: B = [-1/L, 1/L]
        B = np.array([-1/L, 1/L])
        
        # Contribution to stiffness matrix: w * B^T * E * B * A * det(J)
        k += gauss_weights[i] * np.outer(B, B) * E * A * J
    
    # Final stiffness matrix: k = [E*A/L, -E*A/L; -E*A/L, E*A/L]
    k = E * A / L * np.array([[1, -1], [-1, 1]])
    
    return k

def generate_plots_problem1(node_positions, displacements, u_elem1_at_xi, u_elem2_at_xi, u_elem3_at_xi,
                          strain_elem1, strain_elem2, strain_elem3, k1, k2, k3):
    """Generate plots for Problem 1"""
    
    # Plot 1: Bar configuration and displacements
    fig, ax = plt.subplots(figsize=(10, 6))
    
    # x-coordinates for plotting
    x_plot = np.linspace(0, node_positions[-1], 300)
    
    # Displacement function
    u_plot = np.zeros_like(x_plot)
    for i, x in enumerate(x_plot):
        if x <= node_positions[1]:  # Element 1
            xi = 2 * (x - node_positions[0]) / (node_positions[1] - node_positions[0]) - 1
            N1 = (1 - xi) / 2
            N2 = (1 + xi) / 2
            u_plot[i] = N1 * displacements[0] + N2 * displacements[1]
        elif x <= node_positions[2]:  # Element 2
            xi = 2 * (x - node_positions[1]) / (node_positions[2] - node_positions[1]) - 1
            N1 = (1 - xi) / 2
            N2 = (1 + xi) / 2
            u_plot[i] = N1 * displacements[1] + N2 * displacements[2]
        else:  # Element 3
            xi = 2 * (x - node_positions[2]) / (node_positions[3] - node_positions[2]) - 1
            N1 = (1 - xi) / 2
            N2 = (1 + xi) / 2
            u_plot[i] = N1 * displacements[2] + N2 * displacements[3]
    
    # Plot displacement field
    ax.plot(x_plot, u_plot, 'b-', linewidth=2, label='Displacement Field')
    ax.plot(node_positions, displacements, 'ro', markersize=8, label='Nodal Displacements')
    
    # Mark ξ = 0.5 points on each element
    x_elem1_mid = node_positions[0] + 0.75 * (node_positions[1] - node_positions[0])
    x_elem2_mid = node_positions[1] + 0.75 * (node_positions[2] - node_positions[1])
    x_elem3_mid = node_positions[2] + 0.75 * (node_positions[3] - node_positions[2])
    
    ax.plot([x_elem1_mid], [u_elem1_at_xi], 'g*', markersize=10)
    ax.plot([x_elem2_mid], [u_elem2_at_xi], 'g*', markersize=10)
    ax.plot([x_elem3_mid], [u_elem3_at_xi], 'g*', markersize=10)
    
    ax.text(x_elem1_mid, u_elem1_at_xi+0.2, f'ξ=0.5, u={u_elem1_at_xi:.4f}', ha='center', fontsize=9)
    ax.text(x_elem2_mid, u_elem2_at_xi+0.2, f'ξ=0.5, u={u_elem2_at_xi:.4f}', ha='center', fontsize=9)
    ax.text(x_elem3_mid, u_elem3_at_xi+0.2, f'ξ=0.5, u={u_elem3_at_xi:.4f}', ha='center', fontsize=9)
    
    # Add node annotations
    for i, (x, u) in enumerate(zip(node_positions, displacements)):
        ax.annotate(f'Node {i+1}\nu={u}', xy=(x, u), xytext=(x, u-0.5),
                   ha='center', va='top',
                   bbox=dict(boxstyle='round,pad=0.3', fc='white', ec='gray'))
    
    # Add element divisions
    for i in range(1, len(node_positions)-1):
        ax.axvline(x=node_positions[i], color='gray', linestyle='--', alpha=0.7)
    
    # Add element labels
    for i in range(len(node_positions)-1):
        mid_x = (node_positions[i] + node_positions[i+1]) / 2
        ax.text(mid_x, max(displacements)+0.7, f'Element {i+1}', ha='center',
               bbox=dict(boxstyle='round,pad=0.3', fc='lightyellow', ec='gray'))
    
    ax.set_xlabel('Position x')
    ax.set_ylabel('Displacement u')
    ax.set_title('Problem 1: Displacement Field for Three 2-node Elements')
    ax.grid(True, alpha=0.3)
    ax.legend()
    
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, 'ps3_problem1_displacement.png'), dpi=300)
    
    # Plot 2: Element strains
    fig, ax = plt.subplots(figsize=(8, 5))
    
    elements = [f'Element {i+1}' for i in range(3)]
    strains = [strain_elem1, strain_elem2, strain_elem3]
    
    ax.bar(elements, strains, color=['skyblue', 'lightgreen', 'salmon'], edgecolor='black', alpha=0.7)
    
    # Add strain values
    for i, strain in enumerate(strains):
        ax.text(i, strain+0.02, f'{strain:.4f}', ha='center', va='bottom')
    
    ax.set_ylabel('Strain')
    ax.set_title('Problem 1: Element Strains at ξ = 0.5')
    ax.grid(True, axis='y', alpha=0.3)
    
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, 'ps3_problem1_strains.png'), dpi=300)
    
    # Close all figures
    plt.close('all')

if __name__ == "__main__":
    solve_problem1()
