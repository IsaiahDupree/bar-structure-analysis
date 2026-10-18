# Bar Structure Analysis - How to Run

This document provides step-by-step instructions for running the Bar Structure Analysis project.

## Quick Start Guide

### Python Implementation (Recommended)

1. **Prerequisites**: 
   - Python 3.6 or later
   - All dependencies will be installed automatically when you run the main script

2. **Run the complete analysis**:
   ```
   cd python
   python run_all.py
   ```
   
   This will:
   - Check and install any required dependencies
   - Generate all enhanced plots
   - Create a comprehensive analysis report

3. **Output**:
   - Enhanced plots will be saved to the `plots` directory
   - A comprehensive Word document report will be generated in the project root

### MATLAB Implementation (Alternative)

1. **Prerequisites**: 
   - MATLAB R2018b or later
   - No additional toolboxes required

2. **Run the complete analysis**:
   ```
   cd matlab
   run run_all.m
   ```

3. **Output**:
   - Plots will be saved to the `plots` directory
   - Results will be displayed in the MATLAB console

## Directory Structure

- `python/` - Python implementation with automatic dependency management
- `matlab/` - MATLAB implementation 
- `docs/` - Documentation files including detailed analysis and tutorial
- `plots/` - Generated plots (created when you run the analysis)
- `working/` - Auxiliary development files (not needed for running the analysis)

## GitHub Repository

This project is also available on GitHub:
https://github.com/IsaiahDupree/bar-structure-analysis

## Troubleshooting

### Python Implementation

- If you encounter any issues with the Python implementation, ensure you have Python 3.6+ installed
- The script will automatically install dependencies, but you may manually install them:
  ```
  pip install -r python/requirements.txt
  ```
- If plots are not generated, check that the `plots` directory exists and is writable

### MATLAB Implementation

- If you encounter errors in MATLAB, ensure you have a compatible MATLAB version (R2018b+)
- Make sure all files are in their correct directories according to the repository structure

## Contact

For any questions or issues, please open an issue on the GitHub repository.
