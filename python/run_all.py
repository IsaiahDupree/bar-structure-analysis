#!/usr/bin/env python3
"""
Main entry script for Bar Structure Analysis
This script runs both the enhanced plotting and report generation
and automatically handles dependency installation if needed.
"""

import os
import sys
import time
import subprocess
import importlib.util

def check_and_install_dependencies():
    """
    Check if all required dependencies are installed and install them if not
    Returns True if all dependencies are available (either pre-installed or after installation)
    """
    required_packages = {
        'numpy': 'numpy',
        'scipy': 'scipy',
        'matplotlib': 'matplotlib',
        'docx': 'python-docx'
    }
    
    missing_packages = []
    
    print("Checking dependencies...")
    for module_name, package_name in required_packages.items():
        if importlib.util.find_spec(module_name) is None:
            missing_packages.append(package_name)
    
    if missing_packages:
        print(f"Missing dependencies: {', '.join(missing_packages)}")
        should_install = input("Would you like to install these packages now? (y/n): ").strip().lower()
        
        if should_install in ['y', 'yes']:
            print("\nInstalling dependencies...")
            requirements_path = os.path.join(os.path.dirname(__file__), 'requirements.txt')
            
            try:
                subprocess.check_call([sys.executable, '-m', 'pip', 'install', '-r', requirements_path])
                print("\nDependencies installed successfully!")
                return True
            except subprocess.CalledProcessError as e:
                print(f"\nError installing dependencies: {e}")
                print("Please install them manually using the command:")
                print(f"pip install -r {requirements_path}")
                return False
        else:
            print("\nPlease install the required dependencies manually using:")
            print(f"pip install {' '.join(missing_packages)}")
            return False
    else:
        print("All dependencies are already installed!")
        return True

# Import dependencies after checking they're installed
def import_dependencies():
    """
    Import dependencies after ensuring they're installed
    Returns True if all imports succeed
    """
    try:
        global plt, BarAnalysis, create_enhanced_plots, create_comprehensive_document
        import matplotlib.pyplot as plt
        from bar_analysis import BarAnalysis
        from enhanced_plotting import create_enhanced_plots
        from create_comprehensive_document import create_comprehensive_document
        return True
    except ImportError as e:
        print(f"Error importing dependencies: {e}")
        return False

def main():
    print("Bar Structure Analysis - Complete Solution")
    print("==========================================")
    
    # Check dependencies first
    if not check_and_install_dependencies():
        print("\nExiting due to missing dependencies.")
        return
    
    # Import dependencies after installation
    if not import_dependencies():
        print("\nExiting due to import errors.")
        return
    
    # Parameters
    A1, A2, A3 = 200, 100, 50  # mm²
    E1, E2 = 130 * 1000, 200 * 1000  # Convert GPa to MPa
    L = 500  # mm
    F1, F2, F3 = 20 * 1000, 40 * 1000, 20 * 1000  # Convert kN to N
    num_elements_per_segment = 8
    
    try:
        # Create bar analysis instance
        bar = BarAnalysis(
            A1=A1, A2=A2, A3=A3,
            E1=E1, E2=E2,
            L=L,
            F1=F1, F2=F2, F3=F3,
            num_elements_per_segment=num_elements_per_segment
        )
        
        # Step 1: Generate enhanced plots
        print("\nGenerating enhanced plots...")
        plots_dir = "../plots"
        os.makedirs(plots_dir, exist_ok=True)
        create_enhanced_plots(bar, save_dir=plots_dir)
        print(f"Plots saved to {os.path.abspath(plots_dir)}")
        
        # Step 2: Create comprehensive document
        print("\nGenerating comprehensive report...")
        create_comprehensive_document()
        print("Report generated successfully!")
        
        print("\nAll tasks completed successfully.")
        print("For more details, check the generated plots and report.")
    
    except Exception as e:
        import traceback
        print(f"\nError during execution: {e}")
        print("\nDetailed traceback:")
        traceback.print_exc()
        print("\nIf this is a dependency error, make sure all packages are properly installed.")
        print("For help, refer to the README.md file or open an issue on GitHub.")

if __name__ == "__main__":
    main()
