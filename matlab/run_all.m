%% Bar Structure Analysis - Complete MATLAB Solution
% This script runs both the enhanced plotting and report generation

clc;
clear all;
close all;

disp('Bar Structure Analysis - Complete MATLAB Solution');
disp('===============================================');

%% Step 1: Generate enhanced plots
disp('');
disp('Generating enhanced plots...');

% Create plots directory if it doesn't exist
plots_dir = '../plots';
if ~exist(plots_dir, 'dir')
    mkdir(plots_dir);
end

% Run the enhanced plotting script
enhanced_plotting_matlab;

disp(['Plots saved to ' fullfile(pwd, '../plots')]);

%% Step 2: Create enhanced report
disp('');
disp('Generating comprehensive report...');

% Run the report generation
generate_enhanced_report;

disp('Report generated successfully!');
disp('');
disp('All tasks completed successfully.');
disp('For more details, check the generated plots and report.');
