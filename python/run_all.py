#!/usr/bin/env python3
"""
Main entry script for Bar Structure Analysis
This script runs both the enhanced plotting and report generation
"""

import os
import sys
import time
import matplotlib.pyplot as plt

from bar_analysis import BarAnalysis
from enhanced_plotting import create_enhanced_plots
from create_comprehensive_document import create_comprehensive_document

def main():
    print("Bar Structure Analysis - Complete Solution")
    print("==========================================")
    
    # Parameters
    A1, A2, A3 = 200, 100, 50  # mm²
    E1, E2 = 130 * 1000, 200 * 1000  # Convert GPa to MPa
    L = 500  # mm
    F1, F2, F3 = 20 * 1000, 40 * 1000, 20 * 1000  # Convert kN to N
    num_elements_per_segment = 8
    
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

if __name__ == "__main__":
    main()
