%% main_axial_bar.m - Main script for bar structure analysis
% This script analyzes a composite bar with variable cross-section
% using both analytical and finite element methods.

clear; clc;

% Get input parameters
params = inputData();

% Generate initial mesh (4 elements per segment)
[mesh, areas, Es] = generateMesh(params);

% Solve the bar structure using FEM
[U, elemStress] = solveBar(mesh, areas, Es, params);

% Display analytical solution
displayAnalyticalSolution(params);

% Generate visualization plots
postPlots(mesh, U, elemStress, params);

% Error analysis
errorPlot(mesh, elemStress, @analyticStress, params);
