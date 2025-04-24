%% Bar Structure Analysis - MODULAR Solution
% This script runs a comprehensive bar structure analysis with enhanced plotting
% and detailed report generation.
%
% This is the MODULAR version that calls external function files.
% For a standalone version with all functions included, use run_standalone.m

% Clear workspace and command window
clc;
clear all;
close all;

%% Setup Environment
disp('Bar Structure Analysis - Modular Solution');
disp('=======================================');

% Add parent directory to the MATLAB path
parent_dir = fullfile(pwd, '..');
addpath(parent_dir);
fprintf('Added parent directory to path: %s\n', parent_dir);

%% Problem Parameters
A1 = 200;  % mm² - Cross-sectional area 1
A2 = 100;  % mm² - Cross-sectional area 2
A3 = 50;   % mm² - Cross-sectional area 3
E1 = 130;  % GPa - Elastic modulus 1
E2 = 200;  % GPa - Elastic modulus 2
L = 500;   % mm - Segment length
F1 = 20;   % kN - Force 1
F2 = 40;   % kN - Force 2
F3 = 20;   % kN - Force 3
num_elements_per_segment = 8; % Number of elements per segment for FEM

% Display parameters
disp('Problem parameters:');
disp(['- A1 = ' num2str(A1) ' mm², A2 = ' num2str(A2) ' mm², A3 = ' num2str(A3) ' mm²']);
disp(['- E1 = ' num2str(E1) ' GPa, E2 = ' num2str(E2) ' GPa']);
disp(['- L = ' num2str(L) ' mm']);
disp(['- F1 = ' num2str(F1) ' kN, F2 = ' num2str(F2) ' kN, F3 = ' num2str(F3) ' kN']);
disp('');

%% Prepare Output Directory
% Create plots directory if it doesn't exist
plots_dir = fullfile(pwd, '..', 'plots');
if ~exist(plots_dir, 'dir')
    mkdir(plots_dir);
end
fprintf('Using plots directory: %s\n', plots_dir);

%% Step 1: Generate Enhanced Plots
disp('Generating enhanced plots...');

% Call the external plotting function to create enhanced visualizations
enhanced_plotting_matlab(A1, A2, A3, E1, E2, L, F1, F2, F3, num_elements_per_segment, plots_dir);
disp(['Plots generated and saved to ' plots_dir]);
disp('');

%% Step 2: Generate Comprehensive Report
disp('Generating comprehensive report...');
matlab_dir = fullfile(pwd);

% Call the external report generation function
report_file = generate_enhanced_report(A1, A2, A3, E1, E2, L, F1, F2, F3, num_elements_per_segment, matlab_dir);
disp(['Report generated: ' report_file]);
disp('');

%% Execution Complete
disp('All tasks completed successfully.');
disp('For more details, check the generated plots and report.');
