# Bar Structure Analysis

This repository contains a comprehensive analysis of a composite bar structure using both analytical and finite element methods (FEM), implemented in both Python and MATLAB.

## Project Overview

The project analyzes a composite bar fixed at one end, consisting of three segments with different cross-sectional areas and elastic moduli. Forces are applied at the junctions between segments and at the free end. The analysis includes:

- Analytical solution for stress and displacement
- Finite Element Method (FEM) solution using 2-node linear elements
- Comparison between analytical and FEM solutions
- Convergence study showing proper FEM behavior
- Visualization of results with detailed annotations
- Comprehensive reports in text and Word format

## Directory Structure

- `/python/` - Main Python implementation files
- `/matlab/` - Main MATLAB implementation files
- `/plots/` - Generated plots and visualizations
- `/docs/` - Documentation files
- `/working/` - Development and auxiliary files

## Quick Start Guide

### Python Implementation

#### Generate Plots

To generate all plots using the Python implementation:

```bash
cd python
python enhanced_plotting.py
```

#### Generate Reports

To generate a comprehensive report in Word format:

```bash
cd python
python create_comprehensive_document.py
```

### MATLAB Implementation

#### Generate Plots

To generate all plots using the MATLAB implementation:

```matlab
cd matlab
run('enhanced_plotting_matlab.m')
```

#### Generate Reports

To generate a comprehensive report in MATLAB:

```matlab
cd matlab
run('generate_enhanced_report.m')
```

## Main Scripts

### Python

- `python/bar_analysis.py` - Core implementation of both analytical and FEM methods
- `python/enhanced_plotting.py` - Creates enhanced plots with detailed annotations
- `python/create_comprehensive_document.py` - Generates comprehensive Word document report

### MATLAB

- `matlab/main_bar_analysis_enhanced.m` - Main script for enhanced MATLAB analysis
- `matlab/enhanced_plotting_matlab.m` - Creates enhanced plots with detailed annotations
- `matlab/generate_enhanced_report.m` - Generates comprehensive report

## Problem Definition

- **Composite Bar**: Three segments with different properties
- **Cross-sectional Areas**: A₁ = 200 mm², A₂ = 100 mm², A₃ = 50 mm²
- **Elastic Moduli**: E₁ = 130 GPa, E₂ = 200 GPa
- **Segment Length**: L = 500 mm (each segment)
- **Applied Forces**: F₁ = 20 kN, F₂ = 40 kN, F₃ = 20 kN
- **Boundary Condition**: Fixed at the left end (x = 0)

## Output Examples

The analysis generates several visualizations, including:

1. Displacement Field: Comparison of analytical and FEM solutions
2. Stress Field: Comparison of analytical and FEM solutions
3. Convergence Study: Error reduction with mesh refinement
4. Combined Visualization: Bar schematic with stress and displacement

## Requirements

### Python

- Python 3.6+
- NumPy
- SciPy
- Matplotlib
- python-docx (for report generation)

### MATLAB

- MATLAB R2019b or later
- No additional toolboxes required

## Documentation

For more detailed information, refer to:

- `docs/detailed_analysis.md` - In-depth analysis of results
- `docs/tutorial.md` - Step-by-step tutorial on using the code
