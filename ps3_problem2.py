"""
FEA Solver for Problem Set 3 - Problem 2: Two 3-node Elements
"""
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
import os

# Output directory for figures
output_dir = "figures"
if not os.path.exists(output_dir):
    os.makedirs(output_dir)

def solve_problem2():
    """
    Problem 2: Two 3-node elements
    - Node positions: x = 0, 2, 4, 6, 8
    - Displacements at nodes 1, 2, 3, 4, 5 are: 0, -1, 2, -1, 4
    - Element 1: nodes 1, 2, 3 (x = 0, 2, 4)
    - Element 2: nodes 3, 4, 5 (x = 4, 6, 8)
    """
    print("Solving Problem 2: Two 3-node Elements")
    
    # Problem parameters
    node_positions = np.array([0, 2, 4, 6, 8])  # x-coordinates of nodes
    displacements = np.array([0, -1, 2, -1, 4])   # Given displacements at nodes
    
    # 1. Derive shape functions for 3-node elements
    print("\n1. Shape Functions for 3-node elements:")
    print("N₁(ξ) = ξ(ξ-1)/2")
    print("N₂(ξ) = (1+ξ)(1-ξ)")
    print("N₃(ξ) = ξ(ξ+1)/2")
    
    # For reference, ξ goes from -1 to 1 in each element
    # Element 1: nodes 1, 2, 3 (x = 0, 2, 4) → ξ = -1, 0, 1
    # Element 2: nodes 3, 4, 5 (x = 4, 6, 8) → ξ = -1, 0, 1
    
    # Element lengths
    L1 = node_positions[2] - node_positions[0]  # 4
    L2 = node_positions[4] - node_positions[2]  # 4
    
    # 2. Find displacements in each element at ξ = -0.5
    xi = -0.5  # Given ξ value
    
    # Shape function values at ξ = -0.5
    N1_at_xi = xi * (xi - 1) / 2  # 0.125
    N2_at_xi = (1 + xi) * (1 - xi)  # 0.75
    N3_at_xi = xi * (xi + 1) / 2  # -0.125
    
    # Displacements at ξ = -0.5 in each element
    u_elem1_at_xi = N1_at_xi * displacements[0] + N2_at_xi * displacements[1] + N3_at_xi * displacements[2]
    u_elem2_at_xi = N1_at_xi * displacements[2] + N2_at_xi * displacements[3] + N3_at_xi * displacements[4]
    
    print(f"\n2. Displacements at ξ = -0.5:")
    print(f"Element 1: u(ξ=-0.5) = {u_elem1_at_xi:.4f}")
    print(f"Element 2: u(ξ=-0.5) = {u_elem2_at_xi:.4f}")
    
    # 3. Find strains in each element at ξ = -0.5
    # For a 3-node element:
    # dN₁/dξ = ξ - 0.5
    # dN₂/dξ = -2ξ
    # dN₃/dξ = ξ + 0.5
    
    # Derivatives of shape functions at ξ = -0.5
    dN1_dxi = xi - 0.5  # -1.0
    dN2_dxi = -2 * xi  # 1.0
    dN3_dxi = xi + 0.5  # 0.0
    
    # Jacobian for element 1: J = L1/2 = 2
    # Jacobian for element 2: J = L2/2 = 2
    J1 = L1 / 2
    J2 = L2 / 2
    
    # Strain = B * u where B = [dN₁/dx dN₂/dx dN₃/dx]
    # dN/dx = dN/dξ * dξ/dx = dN/dξ * 1/J
    
    # Element 1 strain at ξ = -0.5
    B1 = np.array([dN1_dxi, dN2_dxi, dN3_dxi]) / J1
    strain_elem1 = np.dot(B1, displacements[0:3])
    
    # Element 2 strain at ξ = -0.5
    B2 = np.array([dN1_dxi, dN2_dxi, dN3_dxi]) / J2
    strain_elem2 = np.dot(B2, displacements[2:5])
    
    print(f"\n3. Strains at ξ = -0.5:")
    print(f"Element 1: ε(ξ=-0.5) = {strain_elem1:.4f}")
    print(f"Element 2: ε(ξ=-0.5) = {strain_elem2:.4f}")
    
    # 4. Calculate the stiffness matrix of all elements using integration points ξ = -0.5 and ξ = 0.5
    # For simplicity, we'll assume E (Young's modulus) = 1 and A (Cross-sectional area) = 1
    E = 1.0  # Young's modulus (assumed)
    A = 1.0  # Cross-sectional area (assumed)
    
    # Integration points
    xi_points = np.array([-0.5, 0.5])
    weights = np.array([1.0, 1.0])
    
    # Calculate stiffness matrices
    k1 = calculate_stiffness_matrix(L1, E, A, xi_points, weights)
    k2 = calculate_stiffness_matrix(L2, E, A, xi_points, weights)
    
    print(f"\n4. Element Stiffness Matrices:")
    print(f"Element 1 Stiffness Matrix:")
    print(k1)
    print(f"Element 2 Stiffness Matrix:")
    print(k2)
    
    # Generate plots
    generate_plots_problem2(node_positions, displacements, u_elem1_at_xi, u_elem2_at_xi, 
                           strain_elem1, strain_elem2, k1, k2)
    
    return {
        'node_positions': node_positions,
        'displacements': displacements,
        'element_displacements': [u_elem1_at_xi, u_elem2_at_xi],
        'element_strains': [strain_elem1, strain_elem2],
        'element_stiffness': [k1, k2]
    }

def calculate_stiffness_matrix(L, E, A, xi_points, weights):
    """Calculate stiffness matrix for a 3-node element using numerical integration"""
    k = np.zeros((3, 3))
    
    # Jacobian: J = L/2
    J = L / 2
    
    for i, xi in enumerate(xi_points):
        # Derivatives of shape functions at this integration point
        dN1_dxi = xi - 0.5
        dN2_dxi = -2 * xi
        dN3_dxi = xi + 0.5
        
        # B matrix at this point: B = [dN₁/dx dN₂/dx dN₃/dx]
        B = np.array([dN1_dxi, dN2_dxi, dN3_dxi]) / J
        
        # Contribution to stiffness matrix: w * B^T * E * B * A * det(J)
        k += weights[i] * np.outer(B, B) * E * A * J
    
    return k

def generate_plots_problem2(node_positions, displacements, u_elem1_at_xi, u_elem2_at_xi,
                          strain_elem1, strain_elem2, k1, k2):
    """Generate plots for Problem 2"""
    
    # Plot 1: Bar configuration and displacements
    fig, ax = plt.subplots(figsize=(12, 7))
    
    # x-coordinates for plotting
    x_plot = np.linspace(0, node_positions[-1], 300)
    
    # Displacement function
    u_plot = np.zeros_like(x_plot)
    for i, x in enumerate(x_plot):
        if x <= node_positions[2]:  # Element 1
            # Map x to ξ in the element
            xi = 2 * (x - node_positions[0]) / (node_positions[2] - node_positions[0]) - 1
            
            # Compute shape functions
            N1 = xi * (xi - 1) / 2
            N2 = (1 + xi) * (1 - xi)
            N3 = xi * (xi + 1) / 2
            
            u_plot[i] = N1 * displacements[0] + N2 * displacements[1] + N3 * displacements[2]
            
        else:  # Element 2
            # Map x to ξ in the element
            xi = 2 * (x - node_positions[2]) / (node_positions[4] - node_positions[2]) - 1
            
            # Compute shape functions
            N1 = xi * (xi - 1) / 2
            N2 = (1 + xi) * (1 - xi)
            N3 = xi * (xi + 1) / 2
            
            u_plot[i] = N1 * displacements[2] + N2 * displacements[3] + N3 * displacements[4]
    
    # Plot displacement field
    ax.plot(x_plot, u_plot, 'b-', linewidth=2, label='Displacement Field')
    ax.plot(node_positions, displacements, 'ro', markersize=8, label='Nodal Displacements')
    
    # Mark ξ = -0.5 points on each element
    x_elem1_mid = node_positions[0] + 0.25 * (node_positions[2] - node_positions[0])  # x at ξ = -0.5 in element 1
    x_elem2_mid = node_positions[2] + 0.25 * (node_positions[4] - node_positions[2])  # x at ξ = -0.5 in element 2
    
    ax.plot([x_elem1_mid], [u_elem1_at_xi], 'g*', markersize=10)
    ax.plot([x_elem2_mid], [u_elem2_at_xi], 'g*', markersize=10)
    
    ax.text(x_elem1_mid, u_elem1_at_xi+0.3, f'ξ=-0.5, u={u_elem1_at_xi:.4f}', ha='center', fontsize=9)
    ax.text(x_elem2_mid, u_elem2_at_xi+0.3, f'ξ=-0.5, u={u_elem2_at_xi:.4f}', ha='center', fontsize=9)
    
    # Add node annotations
    for i, (x, u) in enumerate(zip(node_positions, displacements)):
        ax.annotate(f'Node {i+1}\nu={u}', xy=(x, u), xytext=(x, u+0.5*(-1)**(i+1)),
                   ha='center', va='center',
                   bbox=dict(boxstyle='round,pad=0.3', fc='white', ec='gray'))
    
    # Add element divisions
    ax.axvline(x=node_positions[2], color='gray', linestyle='--', alpha=0.7)
    
    # Add element labels
    mid_x1 = (node_positions[0] + node_positions[2]) / 2
    mid_x2 = (node_positions[2] + node_positions[4]) / 2
    
    ax.text(mid_x1, max(displacements)+1.0, f'Element 1 (3-node)', ha='center',
           bbox=dict(boxstyle='round,pad=0.3', fc='lightyellow', ec='gray'))
    ax.text(mid_x2, max(displacements)+1.0, f'Element 2 (3-node)', ha='center',
           bbox=dict(boxstyle='round,pad=0.3', fc='lightyellow', ec='gray'))
    
    ax.set_xlabel('Position x')
    ax.set_ylabel('Displacement u')
    ax.set_title('Problem 2: Displacement Field for Two 3-node Elements')
    ax.grid(True, alpha=0.3)
    ax.legend()
    
    # Set y-axis limits with some padding
    y_min = min(displacements) - 1.0
    y_max = max(displacements) + 1.5
    ax.set_ylim(y_min, y_max)
    
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, 'ps3_problem2_displacement.png'), dpi=300)
    
    # Plot 2: Element strains
    fig, ax = plt.subplots(figsize=(8, 5))
    
    elements = [f'Element {i+1}' for i in range(2)]
    strains = [strain_elem1, strain_elem2]
    
    ax.bar(elements, strains, color=['skyblue', 'lightgreen'], edgecolor='black', alpha=0.7)
    
    # Add strain values
    for i, strain in enumerate(strains):
        ax.text(i, strain+0.05*np.sign(strain), f'{strain:.4f}', ha='center', va='bottom' if strain > 0 else 'top')
    
    ax.set_ylabel('Strain')
    ax.set_title('Problem 2: Element Strains at ξ = -0.5')
    ax.grid(True, axis='y', alpha=0.3)
    
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, 'ps3_problem2_strains.png'), dpi=300)
    
    # Plot 3: Element stiffness matrices visualization as heatmaps
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))
    
    # Display stiffness matrices as heatmaps
    im1 = ax1.imshow(k1, cmap='viridis')
    ax1.set_title('Element 1 Stiffness Matrix')
    ax1.set_xticks(np.arange(3))
    ax1.set_yticks(np.arange(3))
    ax1.set_xticklabels([f'Node {i+1}' for i in range(3)])
    ax1.set_yticklabels([f'Node {i+1}' for i in range(3)])
    
    # Add text annotations to show matrix values
    for i in range(3):
        for j in range(3):
            ax1.text(j, i, f'{k1[i, j]:.3f}', ha='center', va='center', 
                    color='white' if k1[i, j] > 0.5*np.max(k1) else 'black')
    
    im2 = ax2.imshow(k2, cmap='viridis')
    ax2.set_title('Element 2 Stiffness Matrix')
    ax2.set_xticks(np.arange(3))
    ax2.set_yticks(np.arange(3))
    ax2.set_xticklabels([f'Node {i+3}' for i in range(3)])
    ax2.set_yticklabels([f'Node {i+3}' for i in range(3)])
    
    # Add text annotations to show matrix values
    for i in range(3):
        for j in range(3):
            ax2.text(j, i, f'{k2[i, j]:.3f}', ha='center', va='center', 
                    color='white' if k2[i, j] > 0.5*np.max(k2) else 'black')
    
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, 'ps3_problem2_stiffness.png'), dpi=300)
    
    # Close all figures
    plt.close('all')

if __name__ == "__main__":
    solve_problem2()
