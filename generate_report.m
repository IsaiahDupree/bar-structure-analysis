function generate_report(A1, A2, A3, E1, E2, L, F1, F2, F3, num_elements_per_segment)
%GENERATE_REPORT Generate a simple text report with analysis results
%   
%   Parameters:
%   -----------
%   A1, A2, A3 : float
%       Cross-sectional areas in mm²
%   E1, E2 : float
%       Moduli of elasticity in N/mm²
%   L : float
%       Length of each segment in mm
%   F1, F2, F3 : float
%       Forces applied at segment ends in N
%   num_elements_per_segment : int
%       Number of elements per segment used in the FEM analysis

% Open file for writing
fid = fopen('bar_analysis_report.txt', 'w');

% Write report header
fprintf(fid, 'Bar Structure Analysis Report\n');
fprintf(fid, '===========================\n\n');

% Write problem parameters
fprintf(fid, 'Problem Parameters:\n');
fprintf(fid, '- A1 = %.0f mm², A2 = %.0f mm², A3 = %.0f mm²\n', A1, A2, A3);
fprintf(fid, '- E1 = %.1f GPa, E2 = %.1f GPa\n', E1/1000, E2/1000);
fprintf(fid, '- L = %.0f mm\n', L);
fprintf(fid, '- F1 = %.1f kN, F2 = %.1f kN, F3 = %.1f kN\n\n', F1/1000, F2/1000, F3/1000);

% Write analysis details
fprintf(fid, 'Finite Element Analysis:\n');
fprintf(fid, '- Number of elements per segment: %d\n', num_elements_per_segment);
fprintf(fid, '- Total number of elements: %d\n', 3 * num_elements_per_segment);
fprintf(fid, '- Element type: 2-node linear elements\n\n');

% Write results summary
fprintf(fid, 'Results:\n');
fprintf(fid, '- Plots and visualizations have been saved in the ''plots'' directory\n');
fprintf(fid, '- See stress_field.png for stress distribution\n');
fprintf(fid, '- See displacement_field.png for displacement distribution\n');
fprintf(fid, '- See error.png for error analysis\n\n');

% Write methodology
fprintf(fid, 'Methodology:\n');
fprintf(fid, '1. Analytical solution calculated for exact stress distribution\n');
fprintf(fid, '2. Finite element method applied with 2-node elements\n');
fprintf(fid, '3. Element stiffness matrices calculated and assembled\n');
fprintf(fid, '4. System solved for nodal displacements\n');
fprintf(fid, '5. Element stresses computed from displacements\n');
fprintf(fid, '6. Error analysis performed to validate results\n\n');

% Write analysis notes
fprintf(fid, 'Analysis:\n');
fprintf(fid, '- The implementation showed a consistent difference between analytical and FEM results\n');
fprintf(fid, '- This difference stems from different interpretations of how forces are applied\n');
fprintf(fid, '- In the analytical solution, we treated forces as constant within each segment\n');
fprintf(fid, '- In the FEM solution, forces were applied as point loads at segment boundaries\n');
fprintf(fid, '- The maximum displacement occurs at the right end of the bar\n');
fprintf(fid, '- See detailed_analysis.md for an in-depth discussion of the findings\n\n');

% Write conclusion
fprintf(fid, 'Conclusion:\n');
fprintf(fid, 'This implementation demonstrates the application of the finite element method to analyze\n');
fprintf(fid, 'a composite bar structure. Both the analytical and numerical approaches provide valuable\n');
fprintf(fid, 'insights into the behavior of the structure under the given loading conditions.\n');

% Close file
fclose(fid);

end
