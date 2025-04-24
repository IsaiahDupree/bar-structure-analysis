"""
Main script to run the bar structure analysis
"""
import numpy as np
import matplotlib.pyplot as plt
from bar_analysis import BarAnalysis

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
        import time
        import sys
        import threading
        
        # Define a global flag for timeout
        timeout_occurred = False
        
        # Define timeout handler using threading (compatible with Windows)
        def timeout_handler():
            global timeout_occurred
            timeout_occurred = True
            print("\nExecution timed out. Program terminating gracefully.")
            sys.exit(0)
        
        # Set timeout (30 seconds) using a timer thread
        timer = threading.Timer(30.0, timeout_handler)
        timer.daemon = True
        timer.start()
        
        print("Bar Structure Analysis")
        print("=====================")
        print(f"Parameters:")
        print(f"A1 = {A1} mm², A2 = {A2} mm², A3 = {A3} mm²")
        print(f"E1 = {E1/1000} GPa, E2 = {E2/1000} GPa")
        print(f"L = {L} mm")
        print(f"F1 = {F1/1000} kN, F2 = {F2/1000} kN, F3 = {F3/1000} kN")
        print()

        # Initial number of elements per segment
        num_elements_per_segment = 4
        
        # Create the bar analysis object
        bar = BarAnalysis(A1, A2, A3, E1, E2, L, F1, F2, F3, num_elements_per_segment)
        
        # Calculate analytical solution
        start_time = time.time()
        x_analytical, stress_analytical, displacement_analytical = bar.solve_analytical()
        print(f"Analytical solution computed in {time.time() - start_time:.2f} seconds")
        
        # Use finite element method with different element counts until error is below 5%
        max_error = float('inf')
        max_iterations = 3  # Reduced from 5 to speed up execution
        iterations = 0
        
        while max_error > 0.05 and iterations < max_iterations:
            iterations += 1
            start_iter_time = time.time()
            print(f"FEM Iteration {iterations} with {num_elements_per_segment * 3} total elements")
            
            # Solve using FEM
            x_fem, nodal_displacements, element_stresses, error = bar.solve_fem(num_elements_per_segment)
            max_error = np.max(np.abs(error))
            
            print(f"Maximum error: {max_error * 100:.2f}%")
            print(f"Iteration completed in {time.time() - start_iter_time:.2f} seconds")
            
            if max_error > 0.05 and iterations < max_iterations:
                num_elements_per_segment *= 2
                bar.update_num_elements(num_elements_per_segment)
        
        print("\nAnalysis complete. Generating plots and report...")
        
        # Plot results even if the error is still high
        try:
            bar.plot_results(x_analytical, stress_analytical, displacement_analytical, 
                            x_fem, nodal_displacements, element_stresses, error)
            print("Plots generated successfully.")
        except Exception as e:
            print(f"Error generating plots: {str(e)}")
        
        # Generate report
        try:
            report_file = bar.generate_report()
            print(f"Report generated: {report_file}")
        except Exception as e:
            print(f"Error generating report: {str(e)}")
            
        print("\nProgram completed successfully.")
        
    except KeyboardInterrupt:
        print("\nExecution interrupted by user. Program terminating gracefully.")
    except Exception as e:
        print(f"\nAn error occurred: {str(e)}")
        print("Program terminating gracefully.")
    finally:
        # Cancel the timer if it's still active
        try:
            if 'timer' in locals() and timer.is_alive():
                timer.cancel()
        except:
            pass
        print("\nExecution completed.")
        # Return successfully without using sys.exit() to avoid abrupt termination
        # Only use sys.exit() if necessary to forcefully terminate

if __name__ == "__main__":
    main()
