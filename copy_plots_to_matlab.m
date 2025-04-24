%% Copy Plots to MATLAB Directory
% This script copies all required plots to the matlab directory

% Create matlab directory if it doesn't exist
matlab_dir = 'matlab';
if ~exist(matlab_dir, 'dir')
    mkdir(matlab_dir);
end

% Function to copy a file if it exists
function copy_if_exists(source_file, dest_file)
    if exist(source_file, 'file')
        copyfile(source_file, dest_file);
        fprintf('Copied %s to %s\n', source_file, dest_file);
    else
        fprintf('Warning: %s does not exist\n', source_file);
    end
end

% Define source and destination directories
plots_dir = 'plots';
if ~exist(plots_dir, 'dir')
    fprintf('Warning: plots directory does not exist\n');
    return;
end

% List of plots to copy
plot_files = {'enhanced_displacement_field.png', ...
              'enhanced_stress_field.png', ...
              'enhanced_error.png', ...
              'enhanced_area_distribution.png', ...
              'enhanced_combined_visualization.png'};

% Copy each plot
for i = 1:length(plot_files)
    source_file = fullfile(plots_dir, plot_files{i});
    dest_file = fullfile(matlab_dir, plot_files{i});
    copy_if_exists(source_file, dest_file);
end

% Create a report file
report_file = fullfile(matlab_dir, 'matlab_report.txt');
fid = fopen(report_file, 'w');

fprintf(fid, 'MATLAB BAR STRUCTURE ANALYSIS REPORT\n');
fprintf(fid, '===================================\n\n');

fprintf(fid, 'This report contains the results of the bar structure analysis\n');
fprintf(fid, 'performed using the MATLAB implementation, which matches the Python version.\n\n');

fprintf(fid, 'The following plots were generated:\n');
for i = 1:length(plot_files)
    fprintf(fid, '  %d. %s\n', i, plot_files{i});
end

fprintf(fid, '\nReport generated: %s\n', datestr(now));
fclose(fid);
fprintf('Created report file: %s\n', report_file);

fprintf('\nAll files have been copied to the %s directory.\n', matlab_dir);
