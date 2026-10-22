% Adding columns that provide the respective resultant vectors relating to trials in between the trial latency trial in question and the current trial
function [] = add_resultant_vector_columns()
scriptPath = mfilename('fullpath');
[currentDir,~,~] = fileparts(scriptPath);
analysisDir = fileparts(fileparts(currentDir));
data_dir = fullfile(analysisDir, 'data_n_back');
files = dir(fullfile(data_dir, '*.csv'));
X = 10;

for i = 1:numel(files)
    filename = fullfile(files(i).folder, files(i).name);

    fprintf('Reading: %s\n', files(i).name);
    study = readtable(filename, 'Delimiter', ';');
    study_trial_length = height(study);
    res_angle_table = nan(study_trial_length, X);
    res_mag_table = nan(study_trial_length, X);

    for x = 1:X
        for j = 1:study_trial_length
            trial_ang = study.theta(j-1:-1:max(1, j-x+1));  % start is i - 1, so one trial before the current, step size is -1, because we are going back and the stop is max(1, i-x+1) because we don't want to go back further than trial 1 and we don't want to include the trial for N - x
            if isempty(trial_ang)
                continue
            end
            trial_rad = circ_ang2rad(trial_ang);
            res_vector_rad = circ_mean(2 * trial_rad);
            res_angle = circ_rad2ang(res_vector_rad / 2);
            res_mag = circ_r(2 * trial_rad);
            res_angle_table(j,x) = res_angle;
            res_mag_table(j,x) = res_mag;

        end
        study.(sprintf("res_angle_%d", x)) = res_angle_table(:,x);
        study.(sprintf("res_mag_%d", x)) = res_mag_table(:,x);
    end
    writetable(study, filename, 'Delimiter', ';');
end

