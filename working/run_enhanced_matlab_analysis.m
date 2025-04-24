%% Run Enhanced MATLAB Analysis
% This is the main wrapper script to run the enhanced MATLAB implementation
% that matches the Python version with improved visualizations and comprehensive
% documentation including answers to the 8 main questions

% Clear workspace and figures
clear all;
close all;
clc;

disp('=======================================================');
disp('  BAR STRUCTURE ANALYSIS - ENHANCED MATLAB VERSION');
disp('=======================================================');
disp('This script runs the MATLAB implementation with features');
disp('matching the Python version, including enhanced plots');
disp('and comprehensive documentation.');
disp('=======================================================');

% Run the enhanced analysis script
main_bar_analysis_enhanced

% Display a summary of what was generated
disp('=======================================================');
disp('SUMMARY OF GENERATED FILES:');
disp('=======================================================');
disp('1. Enhanced plots in the "matlab" directory:');
disp('   - enhanced_displacement_field.png');
disp('   - enhanced_stress_field.png');
disp('   - enhanced_error.png');
disp('   - enhanced_area_distribution.png');
disp('   - enhanced_combined_visualization.png');
disp(' ');
disp('2. Documentation files:');
disp('   - enhanced_bar_analysis_report.txt');
disp('   - detailed_analysis.md (answers to the 8 main questions)');
disp('   - Composite_Bar_Structure_Analysis_Report.docx (if Report Generator Toolbox is available)');
disp(' ');
disp('=======================================================');
disp('Analysis complete!');
