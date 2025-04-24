"""
Utility functions for bar structure analysis
"""
import numpy as np
import matplotlib.pyplot as plt
import os

def linear_interpolation(x, x0, x1, y0, y1):
    """
    Linear interpolation function
    
    Parameters:
    -----------
    x : float
        Point at which to evaluate
    x0, x1 : float
        Range of x values
    y0, y1 : float
        Range of y values
        
    Returns:
    --------
    float
        Interpolated value at x
    """
    if x1 == x0:
        return y0
    return y0 + (x - x0) * (y1 - y0) / (x1 - x0)

def plot_with_labels(plt_obj, title, xlabel, ylabel, save_path=None):
    """
    Add labels to plot and optionally save
    
    Parameters:
    -----------
    plt_obj : matplotlib.pyplot
        The pyplot object
    title : str
        Plot title
    xlabel : str
        x-axis label
    ylabel : str
        y-axis label
    save_path : str, optional
        Path to save the figure
    """
    plt_obj.title(title)
    plt_obj.xlabel(xlabel)
    plt_obj.ylabel(ylabel)
    plt_obj.grid(True)
    plt_obj.legend()
    
    if save_path:
        directory = os.path.dirname(save_path)
        if directory and not os.path.exists(directory):
            os.makedirs(directory)
        plt_obj.savefig(save_path, dpi=300, bbox_inches='tight')
    
    # Save figure to file and close without showing interactively
    if save_path:
        plt.savefig(save_path, dpi=300, bbox_inches='tight')
    plt.close()

def create_report_with_plots(filename, plots, content):
    """
    Create a simple HTML report with embedded plots
    
    Parameters:
    -----------
    filename : str
        Output filename
    plots : list of str
        List of plot filenames to include
    content : str
        Text content for the report
    """
    html = "<html><head><title>Bar Structure Analysis Report</title>"
    html += "<style>body{font-family:Arial;max-width:800px;margin:auto;padding:20px;line-height:1.6}"
    html += "h1,h2,h3{color:#2c3e50}img{max-width:100%;display:block;margin:20px auto;}"
    html += "pre{background:#f8f8f8;padding:10px;border-radius:5px;overflow-x:auto}</style></head>"
    html += "<body><h1>Bar Structure Analysis Report</h1>"
    
    # Add content
    html += f"<div>{content.replace(chr(10), '<br>')}</div>"
    
    # Add plots
    html += "<h2>Results and Visualizations</h2>"
    for plot in plots:
        if os.path.exists(plot):
            html += f'<figure><img src="{plot}" alt="Plot"><figcaption>{os.path.basename(plot)}</figcaption></figure>'
    
    html += "</body></html>"
    
    with open(filename, 'w') as f:
        f.write(html)
    
    return filename
