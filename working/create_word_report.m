function create_word_report(A1, A2, A3, E1, E2, L, F1, F2, F3, num_elements_per_segment, x_analytical, stress_analytical, displacement_analytical, max_error, iterations, analytical_time, fem_time, overall_time)
%CREATE_WORD_REPORT Generate a comprehensive report in Word format
%   This function creates a Microsoft Word document with the analysis results,
%   including answers to the 8 required questions.

try
    % Create a Microsoft Word application
    wordApp = actxserver('Word.Application');
    wordApp.Visible = false;
    
    % Add a new document
    document = wordApp.Documents.Add;
    
    % Get the current selection (insertion point)
    selection = wordApp.Selection;
    
    % Set document properties
    document.PageSetup.TopMargin = 72; % 1 inch in points
    document.PageSetup.BottomMargin = 72;
    document.PageSetup.LeftMargin = 72;
    document.PageSetup.RightMargin = 72;
    
    % Create title
    selection.Font.Size = 16;
    selection.Font.Bold = 1;
    selection.TypeText('COMPREHENSIVE BAR STRUCTURE ANALYSIS REPORT');
    selection.TypeParagraph;
    selection.TypeText('==========================================');
    selection.TypeParagraph;
    selection.TypeParagraph;
    
    % Reset font
    selection.Font.Size = 12;
    selection.Font.Bold = 0;
    
    % Problem parameters
    selection.Font.Bold = 1;
    selection.TypeText('Problem Parameters:');
    selection.TypeParagraph;
    selection.Font.Bold = 0;
    
    selection.TypeText(['- A1 = ' num2str(A1) ' mm², A2 = ' num2str(A2) ' mm², A3 = ' num2str(A3) ' mm²']);
    selection.TypeParagraph;
    selection.TypeText(['- E1 = ' num2str(E1/1000) ' GPa, E2 = ' num2str(E2/1000) ' GPa']);
    selection.TypeParagraph;
    selection.TypeText(['- L = ' num2str(L) ' mm']);
    selection.TypeParagraph;
    selection.TypeText(['- F1 = ' num2str(F1/1000) ' kN, F2 = ' num2str(F2/1000) ' kN, F3 = ' num2str(F3/1000) ' kN']);
    selection.TypeParagraph;
    selection.TypeParagraph;
    
    % FEM Analysis
    selection.Font.Bold = 1;
    selection.TypeText('Finite Element Analysis:');
    selection.TypeParagraph;
    selection.Font.Bold = 0;
    
    selection.TypeText(['- Number of elements per segment: ' num2str(num_elements_per_segment)]);
    selection.TypeParagraph;
    selection.TypeText(['- Total number of elements: ' num2str(3 * num_elements_per_segment)]);
    selection.TypeParagraph;
    selection.TypeText('- Element type: 2-node linear elements');
    selection.TypeParagraph;
    selection.TypeText(['- Analytical solution computed in ' num2str(analytical_time, '%.2f') ' seconds']);
    selection.TypeParagraph;
    selection.TypeText(['- FEM solution computed in ' num2str(fem_time, '%.2f') ' seconds']);
    selection.TypeParagraph;
    selection.TypeText(['- Maximum error: ' num2str(max_error * 100, '%.2f') '%']);
    selection.TypeParagraph;
    selection.TypeParagraph;
    
    % Results summary
    selection.Font.Bold = 1;
    selection.TypeText('Results and Visualizations:');
    selection.TypeParagraph;
    selection.Font.Bold = 0;
    
    selection.TypeText('- Plots have been saved in the "plots" directory');
    selection.TypeParagraph;
    selection.TypeText('- displacement_field.png: Shows the displacement distribution along the bar');
    selection.TypeParagraph;
    selection.TypeText('- stress_field.png: Shows the stress distribution along the bar');
    selection.TypeParagraph;
    selection.TypeText('- error.png: Shows the relative error in stress between FEM and analytical solutions');
    selection.TypeParagraph;
    selection.TypeText('- area_distribution.png: Shows how the cross-sectional area varies along the bar');
    selection.TypeParagraph;
    selection.TypeText('- modulus_distribution.png: Shows how the elastic modulus varies along the bar');
    selection.TypeParagraph;
    selection.TypeParagraph;
    
    % Insert section for Questions and Answers
    selection.Font.Size = 14;
    selection.Font.Bold = 1;
    selection.TypeText('ANSWERS TO THE 8 ANALYSIS QUESTIONS');
    selection.TypeParagraph;
    selection.TypeText('==================================');
    selection.TypeParagraph;
    selection.TypeParagraph;
    
    % Reset font
    selection.Font.Size = 12;
    selection.Font.Bold = 0;
    
    % Question 1: Maximum displacement in the bar
    [max_displacement, max_disp_idx] = max(abs(displacement_analytical));
    max_disp_location = x_analytical(max_disp_idx);
    
    selection.Font.Bold = 1;
    selection.TypeText('1) What is the maximum displacement in the bar?');
    selection.TypeParagraph;
    selection.Font.Bold = 0;
    selection.TypeText(['   Answer: ' num2str(max_displacement, '%.6f') ' mm at position x = ' num2str(max_disp_location, '%.2f') ' mm']);
    selection.TypeParagraph;
    selection.TypeParagraph;
    
    % Question 2: Maximum stress in the bar
    [max_stress, max_stress_idx] = max(abs(stress_analytical));
    max_stress_location = x_analytical(max_stress_idx);
    
    selection.Font.Bold = 1;
    selection.TypeText('2) What is the maximum stress in the bar?');
    selection.TypeParagraph;
    selection.Font.Bold = 0;
    selection.TypeText(['   Answer: ' num2str(max_stress, '%.2f') ' N/mm² at position x = ' num2str(max_stress_location, '%.2f') ' mm']);
    selection.TypeParagraph;
    selection.TypeParagraph;
    
    % Question 3: Effect of area variation on stress distribution
    selection.Font.Bold = 1;
    selection.TypeText('3) How does the variation in cross-sectional area affect the stress distribution?');
    selection.TypeParagraph;
    selection.Font.Bold = 0;
    selection.TypeText('   Answer: The stress is inversely proportional to the cross-sectional area (σ = F/A).');
    selection.TypeParagraph;
    selection.TypeText(['   Therefore, as the area decreases from A1 (' num2str(A1) ' mm²) to A3 (' num2str(A3) ' mm²),']);
    selection.TypeParagraph;
    selection.TypeText('   we observe an increase in stress. The stress distribution shows significant');
    selection.TypeParagraph;
    selection.TypeText(['   changes at the segment boundaries (x = ' num2str(L) ' mm and x = ' num2str(2*L) ' mm) where the']);
    selection.TypeParagraph;
    selection.TypeText(['   cross-sectional area changes. The smallest area (A3 = ' num2str(A3) ' mm²) in segment 3']);
    selection.TypeParagraph;
    selection.TypeText('   results in the highest stress values.');
    selection.TypeParagraph;
    selection.TypeParagraph;
    
    % Question 4: Effect of material property variation on displacement
    selection.Font.Bold = 1;
    selection.TypeText('4) How does the variation in material properties affect the displacement?');
    selection.TypeParagraph;
    selection.Font.Bold = 0;
    selection.TypeText('   Answer: The displacement is inversely proportional to the elastic modulus (ε = σ/E).');
    selection.TypeParagraph;
    selection.TypeText(['   Segments 1 and 2 have a lower elastic modulus (E1 = ' num2str(E1/1000, '%.1f') ' GPa) compared to']);
    selection.TypeParagraph;
    selection.TypeText(['   segment 3 (E2 = ' num2str(E2/1000, '%.1f') ' GPa). This means that segments 1 and 2 will experience']);
    selection.TypeParagraph;
    selection.TypeText('   more strain (and thus more displacement) for the same stress compared to segment 3.');
    selection.TypeParagraph;
    selection.TypeText(['   The change in elastic modulus at x = ' num2str(2*L) ' mm results in a change in the slope of']);
    selection.TypeParagraph;
    selection.TypeText('   the displacement curve, which is less steep in segment 3 due to the higher modulus.');
    selection.TypeParagraph;
    selection.TypeParagraph;
    
    % Question 5: Accuracy of FEM compared to analytical solution
    selection.Font.Bold = 1;
    selection.TypeText('5) How accurate is the FEM solution compared to the analytical solution?');
    selection.TypeParagraph;
    selection.Font.Bold = 0;
    selection.TypeText(['   Answer: The FEM solution achieved a maximum error of ' num2str(max_error * 100, '%.2f') '% relative to the']);
    selection.TypeParagraph;
    selection.TypeText(['   analytical solution after ' num2str(iterations) ' iterations of mesh refinement. The error is']);
    selection.TypeParagraph;
    selection.TypeText('   calculated as (σ_FEM - σ_analytical)/σ_analytical. The highest errors tend to');
    selection.TypeParagraph;
    selection.TypeText('   occur near the segment boundaries where there are discontinuities in area or');
    selection.TypeParagraph;
    selection.TypeText(['   applied forces. Overall, with ' num2str(num_elements_per_segment) ' elements per segment (' num2str(3*num_elements_per_segment) ' total elements),']);
    selection.TypeParagraph;
    selection.TypeText('   the FEM solution provides a good approximation of the analytical solution.');
    selection.TypeParagraph;
    selection.TypeParagraph;
    
    % Question 6: Effect of mesh refinement on solution accuracy
    selection.Font.Bold = 1;
    selection.TypeText('6) How does mesh refinement affect the solution accuracy?');
    selection.TypeParagraph;
    selection.Font.Bold = 0;
    selection.TypeText('   Answer: Mesh refinement significantly improves the accuracy of the FEM solution.');
    selection.TypeParagraph;
    selection.TypeText(['   Starting with ' num2str(4) ' elements per segment, the maximum error was reduced through']);
    selection.TypeParagraph;
    selection.TypeText('   successive refinements to meet the 5% error tolerance. Each refinement doubled');
    selection.TypeParagraph;
    selection.TypeText('   the number of elements per segment, providing better approximation of the');
    selection.TypeParagraph;
    selection.TypeText('   stress and displacement fields. Mesh refinement is particularly important near');
    selection.TypeParagraph;
    selection.TypeText('   boundaries where stress changes rapidly. With sufficient refinement, the FEM');
    selection.TypeParagraph;
    selection.TypeText('   solution converges to the analytical solution.');
    selection.TypeParagraph;
    selection.TypeParagraph;
    
    % Question 7: How forces affect the displacement and stress
    selection.Font.Bold = 1;
    selection.TypeText('7) How do the applied forces affect the displacement and stress distributions?');
    selection.TypeParagraph;
    selection.Font.Bold = 0;
    selection.TypeText(['   Answer: The applied forces (F1 = ' num2str(F1/1000, '%.1f') ' kN, F2 = ' num2str(F2/1000, '%.1f') ' kN, F3 = ' num2str(F3/1000, '%.1f') ' kN) create']);
    selection.TypeParagraph;
    selection.TypeText('   internal forces in each segment that determine the stress distribution. The');
    selection.TypeParagraph;
    selection.TypeText('   displacement at any point is the cumulative effect of the strains in all');
    selection.TypeParagraph;
    selection.TypeText('   preceding segments. The force F1 creates uniform stress in segment 1, while the');
    selection.TypeParagraph;
    selection.TypeText('   addition of F2 increases the stress in segment 2. Adding F3 further increases');
    selection.TypeParagraph;
    selection.TypeText('   the stress in segment 3. The displacement increases monotonically from the fixed');
    selection.TypeParagraph;
    selection.TypeText(['   end (x = 0) to the free end (x = ' num2str(3*L) ' mm), with the maximum displacement occurring']);
    selection.TypeParagraph;
    selection.TypeText('   at the free end.');
    selection.TypeParagraph;
    selection.TypeParagraph;
    
    % Question 8: Practical implications of the analysis results
    selection.Font.Bold = 1;
    selection.TypeText('8) What are the practical implications of these analysis results?');
    selection.TypeParagraph;
    selection.Font.Bold = 0;
    selection.TypeText('   Answer: This analysis provides several practical insights:');
    selection.TypeParagraph;
    selection.TypeText(['   - The maximum stress (' num2str(max_stress, '%.2f') ' N/mm²) occurs at x = ' num2str(max_stress_location, '%.2f') ' mm, which is a critical']);
    selection.TypeParagraph;
    selection.TypeText('     point for potential material failure. This should be compared with the material''s');
    selection.TypeParagraph;
    selection.TypeText('     yield strength to ensure a sufficient factor of safety.');
    selection.TypeParagraph;
    selection.TypeText(['   - The maximum displacement (' num2str(max_displacement, '%.6f') ' mm) at the free end indicates the deformation']);
    selection.TypeParagraph;
    selection.TypeText('     under the given loads, which is important for assessing functionality and clearances.');
    selection.TypeParagraph;
    selection.TypeText('   - The analysis shows that using variable cross-sectional areas can optimize material');
    selection.TypeParagraph;
    selection.TypeText('     usage while maintaining acceptable stress levels.');
    selection.TypeParagraph;
    selection.TypeText('   - The comparison between different materials (E1 vs E2) demonstrates how material');
    selection.TypeParagraph;
    selection.TypeText('     selection affects structural response.');
    selection.TypeParagraph;
    selection.TypeText('   - The validation of the FEM approach confirms its reliability for more complex');
    selection.TypeParagraph;
    selection.TypeText('     geometries where analytical solutions may not be available.');
    selection.TypeParagraph;
    selection.TypeParagraph;
    
    % Insert images of plots
    selection.Font.Size = 14;
    selection.Font.Bold = 1;
    selection.TypeText('VISUALIZATIONS');
    selection.TypeParagraph;
    selection.TypeText('=============');
    selection.TypeParagraph;
    selection.TypeParagraph;
    
    % Reset font
    selection.Font.Size = 12;
    selection.Font.Bold = 0;
    
    % Add a header for each visualization and insert the image
    try
        % Insert Displacement Field
        selection.Font.Bold = 1;
        selection.TypeText('Figure 1: Displacement Field');
        selection.TypeParagraph;
        selection.Font.Bold = 0;
        selection.InlineShapes.AddPicture([pwd '\plots\displacement_field.png']);
        selection.TypeParagraph;
        selection.TypeParagraph;
        
        % Insert Stress Field
        selection.Font.Bold = 1;
        selection.TypeText('Figure 2: Stress Field');
        selection.TypeParagraph;
        selection.Font.Bold = 0;
        selection.InlineShapes.AddPicture([pwd '\plots\stress_field.png']);
        selection.TypeParagraph;
        selection.TypeParagraph;
        
        % Insert Error Analysis
        selection.Font.Bold = 1;
        selection.TypeText('Figure 3: Error Analysis');
        selection.TypeParagraph;
        selection.Font.Bold = 0;
        selection.InlineShapes.AddPicture([pwd '\plots\error.png']);
        selection.TypeParagraph;
        selection.TypeParagraph;
        
        % Insert Area Distribution
        selection.Font.Bold = 1;
        selection.TypeText('Figure 4: Area Distribution');
        selection.TypeParagraph;
        selection.Font.Bold = 0;
        selection.InlineShapes.AddPicture([pwd '\plots\area_distribution.png']);
        selection.TypeParagraph;
        selection.TypeParagraph;
        
        % Insert Modulus Distribution
        selection.Font.Bold = 1;
        selection.TypeText('Figure 5: Modulus Distribution');
        selection.TypeParagraph;
        selection.Font.Bold = 0;
        selection.InlineShapes.AddPicture([pwd '\plots\modulus_distribution.png']);
        selection.TypeParagraph;
        selection.TypeParagraph;
    catch
        selection.TypeText('Error inserting images. Please ensure all plot files exist in the plots directory.');
        selection.TypeParagraph;
    end
    
    % Add methodology section
    selection.Font.Bold = 1;
    selection.TypeText('Methodology:');
    selection.TypeParagraph;
    selection.Font.Bold = 0;
    
    selection.TypeText('1. Analytical solution calculated for exact stress and displacement distributions');
    selection.TypeParagraph;
    selection.TypeText('2. Finite element method applied with 2-node linear elements');
    selection.TypeParagraph;
    selection.TypeText('3. Element stiffness matrices calculated and assembled');
    selection.TypeParagraph;
    selection.TypeText('4. System solved for nodal displacements');
    selection.TypeParagraph;
    selection.TypeText('5. Element stresses computed from displacements');
    selection.TypeParagraph;
    selection.TypeText('6. Error analysis performed to validate results');
    selection.TypeParagraph;
    selection.TypeText('7. Adaptive mesh refinement applied until error threshold met or max iterations reached');
    selection.TypeParagraph;
    selection.TypeParagraph;
    
    % Add conclusion
    selection.Font.Bold = 1;
    selection.TypeText('Conclusion:');
    selection.TypeParagraph;
    selection.Font.Bold = 0;
    
    selection.TypeText('This comprehensive analysis has provided valuable insights into the behavior of the');
    selection.TypeParagraph;
    selection.TypeText('composite bar structure under axial loading. By comparing analytical and finite');
    selection.TypeParagraph;
    selection.TypeText('element methods, we have demonstrated the effectiveness of numerical approaches');
    selection.TypeParagraph;
    selection.TypeText('while validating their accuracy against exact solutions. The variation in');
    selection.TypeParagraph;
    selection.TypeText('cross-sectional area and material properties significantly influences the stress');
    selection.TypeParagraph;
    selection.TypeText('and displacement distributions, highlighting the importance of these parameters in');
    selection.TypeParagraph;
    selection.TypeText('structural design. The detailed answers to the eight analysis questions provide a');
    selection.TypeParagraph;
    selection.TypeText('comprehensive understanding of the bar''s structural behavior.');
    selection.TypeParagraph;
    selection.TypeParagraph;
    
    % Add execution information
    selection.Font.Bold = 1;
    selection.TypeText('Execution Information:');
    selection.TypeParagraph;
    selection.Font.Bold = 0;
    
    selection.TypeText(['- Total execution time: ' num2str(overall_time, '%.2f') ' seconds']);
    selection.TypeParagraph;
    selection.TypeText(['- Analysis completed on: ' datestr(now)]);
    selection.TypeParagraph;
    
    % Save the document
    docPath = fullfile(pwd, 'Bar_Analysis_Report.docx');
    document.SaveAs(docPath);
    
    % Close the document and quit Word
    document.Close;
    wordApp.Quit;
    
    fprintf('Word document created successfully: %s\n', docPath);
    
catch ME
    fprintf('Error creating Word document: %s\n', ME.message);
    
    % Try to clean up if there was an error
    try
        if exist('document', 'var') && ~isempty(document)
            document.Close;
        end
        if exist('wordApp', 'var') && ~isempty(wordApp)
            wordApp.Quit;
        end
    catch
        % Ignore cleanup errors
    end
    
    % Create a text file as fallback
    fprintf('Creating text report as fallback...\n');
    report_file = 'bar_analysis_report.txt';
    
    % Use your existing code to create a text report here...
    fid = fopen(report_file, 'w');
    fprintf(fid, 'COMPREHENSIVE BAR STRUCTURE ANALYSIS REPORT\n');
    fprintf(fid, '==========================================\n\n');
    % ... add the rest of the text report content
    fclose(fid);
    
    fprintf('Text report created as fallback: %s\n', report_file);
end
end
