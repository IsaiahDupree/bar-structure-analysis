# Python Implementation

This directory contains the Python implementation of the Bar Structure Analysis.

## Main Files

- `bar_analysis.py` - Core implementation of both analytical and FEM methods
- `enhanced_plotting.py` - Creates enhanced plots with detailed annotations
- `create_comprehensive_document.py` - Generates comprehensive Word document report
- `utils.py` - Utility functions used by other modules
- `run_all.py` - Entry script to run all functionality with automatic dependency management
- `requirements.txt` - List of required Python packages

## Usage

To run the complete analysis (generate plots and report):

```bash
python run_all.py
```

The `run_all.py` script will automatically:

1. Check if all required dependencies are installed
2. Offer to install any missing dependencies
3. Generate all plots
4. Create a comprehensive report

To generate only the plots:

```bash
python enhanced_plotting.py
```

To generate only the report:

```bash
python create_comprehensive_document.py
```

## Requirements

- Python 3.6+
- NumPy
- SciPy
- Matplotlib
- python-docx (for report generation)

All required packages are listed in `requirements.txt` and will be automatically installed by the `run_all.py` script if they are missing.
