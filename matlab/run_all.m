%% Bar Structure Analysis - Complete MATLAB Solution
% This script runs both the enhanced plotting and report generation

clc;
clear all;
close all;

% Add parent directory to path to find all required functions
parent_dir = fullfile(pwd, '..');
addpath(parent_dir);
disp(['Added parent directory to path: ' parent_dir]);

disp('Bar Structure Analysis - Complete MATLAB Solution');
disp('===============================================');

%% Define parameters for the analysis
% Problem parameters
A1 = 200;  % mm²
A2 = 100;  % mm²
A3 = 50;   % mm²
E1 = 130;  % GPa
E2 = 200;  % GPa
L = 500;   % mm
F1 = 20;   % kN
F2 = 40;   % kN
F3 = 20;   % kN
num_elements_per_segment = 8; % Number of elements per segment

disp('Problem parameters:');
disp(['- A1 = ' num2str(A1) ' mm², A2 = ' num2str(A2) ' mm², A3 = ' num2str(A3) ' mm²']);
disp(['- E1 = ' num2str(E1) ' GPa, E2 = ' num2str(E2) ' GPa']);
disp(['- L = ' num2str(L) ' mm']);
disp(['- F1 = ' num2str(F1) ' kN, F2 = ' num2str(F2) ' kN, F3 = ' num2str(F3) ' kN']);
disp('');

%% Step 1: Generate enhanced plots
disp('Generating enhanced plots...');

% Create plots directory if it doesn't exist
plots_dir = '../plots';
if ~exist(plots_dir, 'dir')
    mkdir(plots_dir);
end

% Run the enhanced plotting script with parameters
enhanced_plotting_matlab(A1, A2, A3, E1, E2, L, F1, F2, F3, num_elements_per_segment, plots_dir);

disp(['Plots saved to ' fullfile(pwd, '../plots')]);
disp('');

%% Step 2: Create enhanced report
disp('Generating comprehensive report...');

% Run the report generation with parameters
report_file = generate_enhanced_report(A1, A2, A3, E1, E2, L, F1, F2, F3, num_elements_per_segment, 'matlab');

disp(['Report generated successfully: ' report_file]);
disp('');
disp('All tasks completed successfully.');
disp('For more details, check the generated plots and report.');
