function [] = paired_t_test_res_vector()
%% estimating error scatter at iso, mid and ortho using a weighted flexible-degree polynomial fitting
scriptPath = mfilename('fullpath');
[currentDir,~,~] = fileparts(scriptPath);
analysisDir = fileparts(fileparts(currentDir));
output_dir = fullfile(analysisDir, 'results');
table_dir = fullfile(output_dir, 'tables_first_half');
tbl = readtable(fullfile(analysisDir,'results','SD_ma_master_table_n_back_first_half.csv'),'Delimiter',',');

% Initializing some variables
X = 10;
error_variable = 'error_ori_deb';
delta_lower = 20;
delta_upper = 50;

table_dir_window_size = fullfile(table_dir, sprintf('%d_to_%d_degree_window', delta_lower, delta_upper));
if ~isfolder(table_dir_window_size)
    mkdir(table_dir_window_size);
end
for x = 2:(X-1)
    delta_variable = sprintf('delta_%d',x);
    res_angle_variable = sprintf('res_angle_%d', x + 1);

    mean_error_pairs = nan(0, 2);

    % Loop over studies
    studies = unique(tbl.studynum);

    for i = 1:numel(studies)
        tbl_i = tbl(tbl.studynum==studies(i),:);
        % Specify the window you want
        tbl_i = tbl_i(abs(tbl_i.(delta_variable)) > delta_lower & abs(tbl_i.(delta_variable)) < delta_upper, :);
        err = tbl_i.(error_variable);
        err(tbl_i.(delta_variable)<0 & tbl_i.(delta_variable)~=-90) = -err(tbl_i.(delta_variable)<0 & tbl_i.(delta_variable)~=-90);
        tbl_i.(error_variable) = err;

        % Loop over experiments and conditions
        experiments = unique(tbl_i.expnum);

        for k_index = 1:numel(experiments)
            k = experiments(k_index);
            tbl_i_k = tbl_i(tbl_i.expnum==k,:);
            cond = unique(tbl_i_k.cond);
            ncond = numel(cond);

            for j = 1:ncond
                tbl_i_k_j = tbl_i_k(tbl_i_k.cond==cond(j),:);
                obs = unique(tbl_i_k_j.obs);
                nobs = numel(obs);

                for o = 1:nobs
                    tbl_i_k_j_o = tbl_i_k_j(tbl_i_k_j.obs==obs(o),:);

                    stimulus = string(tbl_i_k_j_o.stimulus(1));
                    if strcmpi(strtrim(stimulus), "Motion")
                        period = 360;
                    else
                        period = 180;
                    end

                    theta = tbl_i_k_j_o.theta;
                    res_angles = tbl_i_k_j_o.(res_angle_variable);

                    angle_difference = mod(res_angles - theta + period/2, period) - period/2;
                    valid = isfinite(res_angles) & isfinite(theta);

                    same_direction = valid & abs(angle_difference) <= 20;
                    different_direction = valid & abs(angle_difference) > 20;

                    delta_values = tbl_i_k_j_o.(delta_variable);
                    error_values = tbl_i_k_j_o.(error_variable);
                    abs_delta_same = abs(delta_values(same_direction));
                    abs_delta_diff = abs(delta_values(different_direction));
                    folded_error_same = error_values(same_direction);
                    folded_error_diff = error_values(different_direction);

                    valid_same = ~isnan(abs_delta_same) & ~isnan(folded_error_same);
                    folded_error_same = folded_error_same(valid_same);

                    valid_diff = ~isnan(abs_delta_diff) & ~isnan(folded_error_diff);
                    folded_error_diff = folded_error_diff(valid_diff);

                    folded_error_same_avg = mean(folded_error_same);
                    folded_error_diff_avg = mean(folded_error_diff);
                    mean_error_pairs(end+1, :) = [folded_error_same_avg, folded_error_diff_avg];

                end
            end
        end
    end

    valid_pairs = all(isfinite(mean_error_pairs), 2);
    paired_means = mean_error_pairs(valid_pairs, :);

    if size(paired_means, 1) < 2
        warning('Not enough valid pairs for N-%d.', x);
        continue
    end

    [h, p, ci, stats] = ttest(paired_means(:,1), paired_means(:,2));
    mean_difference = mean(paired_means(:,1) - paired_means(:,2));

    stats_table = table(x, size(paired_means, 1), mean_difference, h, p, ci(1), ci(2), stats.tstat, stats.df, 'VariableNames', {'Lag', 'NPairs', 'MeanDifference', 'h', 'PValue', 'CILower', 'CIUpper', 'TValue', 'DF'});
    filename_stats_table = fullfile(table_dir_window_size, sprintf('paired_bias_stats_N_%d_folded.csv', x));
    writetable(stats_table, filename_stats_table, 'Delimiter', ';');
end
end