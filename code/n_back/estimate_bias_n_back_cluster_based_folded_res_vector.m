function [pos_clusters_same, neg_clusters_same, pos_cluster_p_same, neg_cluster_p_same, pos_clusters_diff, neg_clusters_diff, pos_cluster_p_diff, neg_cluster_p_diff] = estimate_bias_n_back_cluster_based_folded_res_vector()
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
alpha = 0.05;
bin_sizes = 11;        %%%%%%%%%%%%%%%%%
n_permutations = 500;       %%%%%%%%%%%%%%%%%

group_bias_same = cell(1, X);
group_bias_diff = cell(1, X);
t_values_same = nan(91, X);
t_values_diff = nan(91, X);
pos_mask_same = false(91, X);
neg_mask_same = false(91, X);
pos_mask_diff = false(91, X);
neg_mask_diff = false(91, X);

pos_clusters_same = cell(1, X);
neg_clusters_same = cell(1, X);
pos_clusters_diff = cell(1, X);
neg_clusters_diff = cell(1, X);

pos_counts_same = cell(1, X);
neg_counts_same = cell(1, X);
pos_counts_diff = cell(1, X);
neg_counts_diff = cell(1, X);

pos_cluster_p_same = cell(1, X);
neg_cluster_p_same = cell(1, X);
pos_cluster_p_diff = cell(1, X);
neg_cluster_p_diff = cell(1, X);

% % Correlation between resultant vector and error bias, taking into account the resultant vector magnitude
% for x = 2:(X-1)
%     delta_variable = sprintf('delta_%d',x);
%     res_angle_variable = sprintf('res_angle_%d', x + 1);
%     res_mag_variable = sprintf('res_mag_%d', x + 1);

%     correlations = [];
%     trial_counts = [];

%     % Loop over studies
%     studies = unique(tbl.studynum);

%     for i = 1:numel(studies)
%         tbl_i = tbl(tbl.studynum==studies(i),:);
%         tbl_i = tbl_i(tbl_i.(delta_variable)>-90 & tbl_i.(delta_variable)<=90,:);

%         err = tbl_i.(error_variable);
%         err(tbl_i.(delta_variable)<0 & tbl_i.(delta_variable)~=-90) = -err(tbl_i.(delta_variable)<0 & tbl_i.(delta_variable)~=-90);
%         tbl_i.(error_variable) = err;

%         % Loop over experiments and conditions
%         experiments = unique(tbl_i.expnum);

%         for k_index = 1:numel(experiments)
%             k = experiments(k_index);
%             tbl_i_k = tbl_i(tbl_i.expnum==k,:);
%             cond = unique(tbl_i_k.cond);
%             ncond = numel(cond);

%             for j = 1:ncond
%                 tbl_i_k_j = tbl_i_k(tbl_i_k.cond==cond(j),:);
%                 obs = unique(tbl_i_k_j.obs);
%                 nobs = numel(obs);

%                 for o = 1:nobs
%                     tbl_i_k_j_o = tbl_i_k_j(tbl_i_k_j.obs==obs(o),:);
%                     folded_error = tbl_i_k_j_o.(error_variable);

%                     stimulus = string(tbl_i_k_j_o.stimulus(1));
%                     if strcmpi(strtrim(stimulus), "Motion")
%                         period = 360;
%                     else
%                         period = 180;
%                     end
%                     angle_multiplier = 360 / period;

%                     theta = tbl_i_k_j_o.theta;
%                     res_angles = tbl_i_k_j_o.(res_angle_variable);
%                     res_magnitudes = tbl_i_k_j_o.(res_mag_variable);

%                     alignment = cosd(angle_multiplier * (res_angles - theta));
%                     % weighted_alignment = res_magnitudes .* alignment;
%                     weighted_alignment = alignment;

%                     valid = isfinite(weighted_alignment) & isfinite(folded_error);
%                     alignment_valid = weighted_alignment(valid);
%                     error_valid = folded_error(valid);

%                     r = NaN;

%                     if numel(alignment_valid) >= 3 && std(alignment_valid) > 0 && std(error_valid) > 0
%                         r = corr(alignment_valid(:), error_valid(:), 'Type', 'Pearson');
%                     end

%                     correlations(end+1, 1) = r;
%                     trial_counts(end+1, 1) = numel(alignment_valid);

%                 end
%             end
%         end
%     end

%     valid = isfinite(correlations) & abs(correlations) < 1 & trial_counts > 3;

%     r_values = correlations(valid);
%     weights = trial_counts(valid) - 3;

%     if isempty(r_values)
%         overall_r = NaN;
%     else
%         z_values = atanh(r_values);
%         mean_z = sum(weights .* z_values) / sum(weights);
%         overall_r = tanh(mean_z);
%     end

%     fprintf('The overall correlation for N-%d = %.4f\n', x, overall_r);
%     correlation_table = table(correlations, trial_counts, z_values, 'VariableNames', {'Correlation', 'NTrials', 'ZValues'});
%     writetable(correlation_table, fullfile(table_dir, sprintf('alignment_error_correlations_N_%d.csv', x)), 'Delimiter', ';');

% end
% end


% Computing bias for true delta values
for bin_size = bin_sizes
    table_dir_bin_size = fullfile(table_dir, sprintf('%d_degree_window', bin_size));
    if ~isfolder(table_dir_bin_size)
        mkdir(table_dir_bin_size);
    end
    for x = 2:X
        delta_variable = sprintf('delta_%d',x);
        res_angle_variable = sprintf('res_angle_%d', x);

        mvav_bias_aggregated_same = nan(91, 0);
        mvav_bias_aggregated_diff = nan(91, 0);

        % Loop over studies
        studies = unique(tbl.studynum);

        for i = 1:numel(studies)
            tbl_i = tbl(tbl.studynum==studies(i),:);
            tbl_i = tbl_i(tbl_i.(delta_variable)>-90 & tbl_i.(delta_variable)<=90,:);

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

                        % theta_n_back = tbl_i_k_j_o.theta + tbl_i_k_j_o.(delta_variable);
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
                        abs_delta_same = abs_delta_same(valid_same);
                        folded_error_same = folded_error_same(valid_same);

                        valid_diff = ~isnan(abs_delta_diff) & ~isnan(folded_error_diff);
                        abs_delta_diff = abs_delta_diff(valid_diff);
                        folded_error_diff = folded_error_diff(valid_diff);

                        % abs_delta = abs(tbl_i_k_j_o.(delta_variable));
                        % folded_error = tbl_i_k_j_o.(error_variable);

                        % valid = ~isnan(abs_delta) & ~isnan(folded_error);
                        % abs_delta = abs_delta(valid);
                        % folded_error = folded_error(valid);

                        if ~isempty(abs_delta_same)
                            biass_same = grpstats(folded_error_same,abs_delta_same,'mean');
                            [counts_same, ~] = groupcounts(abs_delta_same);
                            delta_fit_same = unique(abs_delta_same);

                            mvav_same = nan(91,1); mvav_same(delta_fit_same+1) = biass_same; mvav_same = [-mvav_same(end-1:-1:2);mvav_same]; mvav_same = repmat(mvav_same,3,1);
                            weights_same = nan(91,1); weights_same(delta_fit_same+1) = counts_same; weights_same = [weights_same(end-1:-1:2);weights_same]; weights_same = repmat(weights_same,3,1);

                            weighted_moving_avg_same = movsum(mvav_same.*weights_same,bin_size,'omitnan')./movsum(weights_same,bin_size,'omitnan');
                            weighted_moving_avg_same = weighted_moving_avg_same(180+[90:180]);
                            mvav_bias_aggregated_same = [mvav_bias_aggregated_same weighted_moving_avg_same];
                        end

                        if ~isempty(abs_delta_diff)
                            biass_diff = grpstats(folded_error_diff,abs_delta_diff,'mean');
                            [counts_diff, ~] = groupcounts(abs_delta_diff);
                            delta_fit_diff = unique(abs_delta_diff);

                            mvav_diff = nan(91,1); mvav_diff(delta_fit_diff+1) = biass_diff; mvav_diff = [-mvav_diff(end-1:-1:2);mvav_diff]; mvav_diff = repmat(mvav_diff,3,1);
                            weights_diff = nan(91,1); weights_diff(delta_fit_diff+1) = counts_diff; weights_diff = [weights_diff(end-1:-1:2);weights_diff]; weights_diff = repmat(weights_diff,3,1);

                            weighted_moving_avg_diff = movsum(mvav_diff.*weights_diff,bin_size,'omitnan')./movsum(weights_diff,bin_size,'omitnan');
                            weighted_moving_avg_diff = weighted_moving_avg_diff(180+[90:180]);
                            mvav_bias_aggregated_diff = [mvav_bias_aggregated_diff weighted_moving_avg_diff];
                        end

                        % % Calculate mean signed error and trial count per delta                   
                        % biass = grpstats(folded_error,abs_delta,'mean');
                        % [counts, ~] = groupcounts(abs_delta);
                        % delta_fit = unique(abs_delta);

                        % % Aggregated moving average bias
                        % mvav = nan(91,1); mvav(delta_fit+1) = biass; mvav = [-mvav(end-1:-1:2);mvav]; mvav=repmat(mvav,3,1);
                        % weights = nan(91,1); weights(delta_fit+1) = counts; weights = [weights(end-1:-1:2);weights]; weights=repmat(weights,3,1);
                        % weighted_moving_avg = movsum(mvav.*weights,bin_size,'omitnan')./movsum(weights,bin_size,'omitnan');
                        % weighted_moving_avg = weighted_moving_avg(180+[90:180]);
                        % mvav_bias_aggregated = [mvav_bias_aggregated weighted_moving_avg];

                    end
                end
            end
        end

        group_bias_same{x} = mean(mvav_bias_aggregated_same,2,'omitnan');
        group_bias_diff{x} = mean(mvav_bias_aggregated_diff,2,'omitnan');
        group_delta = (0:90)';

        % One test per row, across participant columns, against zero
        n_same = sum(isfinite(mvav_bias_aggregated_same), 2);
        n_diff = sum(isfinite(mvav_bias_aggregated_diff), 2);
        fprintf('N-%d: points with fewer than 2 observations: same=%d, different=%d\n', x, sum(n_same < 2), sum(n_diff < 2));

        [~, ~, ~, stats_same] = ttest(mvav_bias_aggregated_same, 0, 'Dim', 2);
        [~, ~, ~, stats_diff] = ttest(mvav_bias_aggregated_diff, 0, 'Dim', 2);

        t_values_same(:, x) = stats_same.tstat;
        t_values_diff(:, x) = stats_diff.tstat;
        t_threshold_same = tinv(1 - alpha/2, stats_same.df);
        t_threshold_diff = tinv(1 - alpha/2, stats_diff.df);

        pos_mask_same(:, x) = stats_same.tstat > t_threshold_same;
        neg_mask_same(:, x) = stats_same.tstat < -t_threshold_same;
        pos_mask_diff(:, x) = stats_diff.tstat > t_threshold_diff;
        neg_mask_diff(:, x) = stats_diff.tstat < -t_threshold_diff;

        pos_clusters_same{x} = find_clusters(pos_mask_same(:, x), t_values_same(:, x));
        neg_clusters_same{x} = find_clusters(neg_mask_same(:, x), t_values_same(:, x));
        pos_clusters_diff{x} = find_clusters(pos_mask_diff(:, x), t_values_diff(:, x));
        neg_clusters_diff{x} = find_clusters(neg_mask_diff(:, x), t_values_diff(:, x));

        % figure;
        % plot(group_delta, t_values(:, x), 'LineWidth', 1.5);
        % yline(0, '--');
        % xline(0, '--');
        % xlabel('Angular difference (degrees)');
        % ylabel('t-value');
        % title(sprintf('Group bias t-values: N-%d', x));

        % figure;
        % plot(group_delta, group_bias, 'LineWidth', 1.5);
        % yline(0, '--');
        % xline(0, '--');

        % xlabel('Angular difference (degrees)');
        % ylabel('Mean signed error (degrees)');
        % title(sprintf('Group bias: N-%d', x));
        % grid on;

        bias_table_same = table(group_delta, group_bias_same{x}, 'VariableNames', {'DeltaDegrees', 'MeanBias'});
        filename_same = fullfile(table_dir_bin_size, sprintf('group_bias_N_%d_folded_theta_20_same_direction_res_vector.csv', x));
        writetable(bias_table_same, filename_same, 'Delimiter', ';');
        bias_table_diff = table(group_delta, group_bias_diff{x}, 'VariableNames', {'DeltaDegrees', 'MeanBias'});
        filename_diff = fullfile(table_dir_bin_size, sprintf('group_bias_N_%d_folded_theta_20_diff_direction_res_vector.csv', x));
        writetable(bias_table_diff, filename_diff, 'Delimiter', ';');

    end

    % Computing bias values for shuffled delta values
    for x = 2:X
        delta_variable = sprintf('delta_%d',x);
        res_angle_variable = sprintf('res_angle_%d', x);
        lag_timer = tic;
        % These are the clusters from the original, unshuffled data
        observed_pos_mass_same = [pos_clusters_same{x}.mass];
        observed_neg_mass_same = [neg_clusters_same{x}.mass];
        observed_pos_mass_diff = [pos_clusters_diff{x}.mass];
        observed_neg_mass_diff = [neg_clusters_diff{x}.mass];

        pos_counts_same{x} = zeros(size(observed_pos_mass_same));
        neg_counts_same{x} = zeros(size(observed_neg_mass_same));
        pos_counts_diff{x} = zeros(size(observed_pos_mass_diff));
        neg_counts_diff{x} = zeros(size(observed_neg_mass_diff));

        fprintf('\nStarting %d permutations for N-%d...\n',n_permutations,x);

        % Permutation test on shuffled data
        for permutation = 1:n_permutations
            mvav_bias_aggregated_same = nan(91, 0);
            mvav_bias_aggregated_diff = nan(91, 0);

            t_values_perm_same = nan(91, X);
            pos_mask_perm_same = false(91, X);
            neg_mask_perm_same = false(91, X);
            t_values_perm_diff = nan(91, X);
            pos_mask_perm_diff = false(91, X);
            neg_mask_perm_diff = false(91, X);

            % Loop over studies
            studies = unique(tbl.studynum);

            for i = 1:numel(studies)
                tbl_i = tbl(tbl.studynum==studies(i),:);
                tbl_i = tbl_i(tbl_i.(delta_variable)>-90 & tbl_i.(delta_variable)<=90,:);

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

                            theta_n_back = tbl_i_k_j_o.theta + tbl_i_k_j_o.(delta_variable);
                            res_angles = tbl_i_k_j_o.(res_angle_variable);

                            angle_difference = mod(res_angles - theta_n_back + period/2, period) - period/2;
                            valid = isfinite(res_angles) & isfinite(theta_n_back);

                            same_direction = valid & abs(angle_difference) <= 45;
                            different_direction = valid & abs(angle_difference) > 45;

                            delta_values = tbl_i_k_j_o.(delta_variable);
                            error_values = tbl_i_k_j_o.(error_variable);
                            delta_same = delta_values(same_direction);
                            delta_diff = delta_values(different_direction);
                            folded_error_same = error_values(same_direction);
                            folded_error_diff = error_values(different_direction);

                            valid_same = ~isnan(delta_same) & ~isnan(folded_error_same);
                            delta_same = delta_same(valid_same);
                            folded_error_same = folded_error_same(valid_same);

                            valid_diff = ~isnan(delta_diff) & ~isnan(folded_error_diff);
                            delta_diff = delta_diff(valid_diff);
                            folded_error_diff = folded_error_diff(valid_diff);

                            %%%
                            shuffled_delta_same = delta_same(randperm(numel(delta_same)));
                            shuffled_delta_diff = delta_diff(randperm(numel(delta_diff)));

                            flip_same = shuffled_delta_same < 0 & shuffled_delta_same ~= -90;
                            folded_error_same(flip_same) = -folded_error_same(flip_same);
                            abs_shuffled_delta_same = abs(shuffled_delta_same);

                            flip_diff = shuffled_delta_diff < 0 & shuffled_delta_diff ~= -90;
                            folded_error_diff(flip_diff) = -folded_error_diff(flip_diff);
                            abs_shuffled_delta_diff = abs(shuffled_delta_diff);

                            if ~isempty(abs_shuffled_delta_same)
                                biass_same = grpstats(folded_error_same,abs_shuffled_delta_same,'mean');
                                [counts_same, ~] = groupcounts(abs_shuffled_delta_same);
                                delta_fit_same = unique(abs_shuffled_delta_same);

                                mvav_same = nan(91,1); mvav_same(delta_fit_same+1) = biass_same; mvav_same = [-mvav_same(end-1:-1:2);mvav_same]; mvav_same = repmat(mvav_same,3,1);
                                weights_same = nan(91,1); weights_same(delta_fit_same+1) = counts_same; weights_same = [weights_same(end-1:-1:2);weights_same]; weights_same = repmat(weights_same,3,1);

                                weighted_moving_avg_same = movsum(mvav_same.*weights_same,bin_size,'omitnan')./movsum(weights_same,bin_size,'omitnan');
                                weighted_moving_avg_same = weighted_moving_avg_same(180+[90:180]);
                                mvav_bias_aggregated_same = [mvav_bias_aggregated_same weighted_moving_avg_same];
                            end

                            if ~isempty(abs_shuffled_delta_diff)
                                biass_diff = grpstats(folded_error_diff,abs_shuffled_delta_diff,'mean');
                                [counts_diff, ~] = groupcounts(abs_shuffled_delta_diff);
                                delta_fit_diff = unique(abs_shuffled_delta_diff);

                                mvav_diff = nan(91,1); mvav_diff(delta_fit_diff+1) = biass_diff; mvav_diff = [-mvav_diff(end-1:-1:2);mvav_diff]; mvav_diff = repmat(mvav_diff,3,1);
                                weights_diff = nan(91,1); weights_diff(delta_fit_diff+1) = counts_diff; weights_diff = [weights_diff(end-1:-1:2);weights_diff]; weights_diff = repmat(weights_diff,3,1);

                                weighted_moving_avg_diff = movsum(mvav_diff.*weights_diff,bin_size,'omitnan')./movsum(weights_diff,bin_size,'omitnan');
                                weighted_moving_avg_diff = weighted_moving_avg_diff(180+[90:180]);
                                mvav_bias_aggregated_diff = [mvav_bias_aggregated_diff weighted_moving_avg_diff];
                            end

                        end
                    end
                end
            end

            % One test per row, across participant columns, against zero
            [~, ~, ~, stats_same] = ttest(mvav_bias_aggregated_same, 0, 'Dim', 2);
            [~, ~, ~, stats_diff] = ttest(mvav_bias_aggregated_diff, 0, 'Dim', 2);

            t_values_perm_same(:, x) = stats_same.tstat;
            t_values_perm_diff(:, x) = stats_diff.tstat;

            t_threshold_same = tinv(1 - alpha/2, stats_same.df);
            t_threshold_diff = tinv(1 - alpha/2, stats_diff.df);

            pos_mask_perm_same(:, x) = stats_same.tstat > t_threshold_same;
            neg_mask_perm_same(:, x) = stats_same.tstat < -t_threshold_same;
            pos_mask_perm_diff(:, x) = stats_diff.tstat > t_threshold_diff;
            neg_mask_perm_diff(:, x) = stats_diff.tstat < -t_threshold_diff;

            pos_clusters_perm_same{x} = find_clusters(pos_mask_perm_same(:, x), t_values_perm_same(:, x));
            neg_clusters_perm_same{x} = find_clusters(neg_mask_perm_same(:, x), t_values_perm_same(:, x));
            pos_clusters_perm_diff{x} = find_clusters(pos_mask_perm_diff(:, x), t_values_perm_diff(:, x));
            neg_clusters_perm_diff{x} = find_clusters(neg_mask_perm_diff(:, x), t_values_perm_diff(:, x));

            % Largest mass in this permutation; zero if no clusters exist
            max_mass_same = max([0, pos_clusters_perm_same{x}.mass, neg_clusters_perm_same{x}.mass]);
            max_mass_diff = max([0, pos_clusters_perm_diff{x}.mass, neg_clusters_perm_diff{x}.mass]);
            
            % Update one count for each observed cluster
            pos_counts_same{x} = pos_counts_same{x} + (max_mass_same >= observed_pos_mass_same);
            neg_counts_same{x} = neg_counts_same{x} + (max_mass_same >= observed_neg_mass_same);
            pos_counts_diff{x} = pos_counts_diff{x} + (max_mass_diff >= observed_pos_mass_diff);
            neg_counts_diff{x} = neg_counts_diff{x} + (max_mass_diff >= observed_neg_mass_diff);

            if mod(permutation,50) == 0
                fprintf('N-%d: %d/%d permutations completed (%.1f minutes)\n',x,permutation,n_permutations,toc(lag_timer)/60);
            end
        end

        pos_cluster_p_same{x} = (pos_counts_same{x} + 1) / (n_permutations + 1);
        neg_cluster_p_same{x} = (neg_counts_same{x} + 1) / (n_permutations + 1);
        pos_cluster_p_diff{x} = (pos_counts_diff{x} + 1) / (n_permutations + 1);
        neg_cluster_p_diff{x} = (neg_counts_diff{x} + 1) / (n_permutations + 1);

        % An attempt to have results saved every lag, in order to not lose progress when running fails halfway through
        save_cluster_csv( ...
            pos_clusters_same(1:x), pos_cluster_p_same(1:x), ...
            fullfile(table_dir_bin_size, 'positive_clusters_folded_theta_20_same_res_vector.csv'));

        save_cluster_csv( ...
            neg_clusters_same(1:x), neg_cluster_p_same(1:x), ...
            fullfile(table_dir_bin_size, 'negative_clusters_folded_theta_20_same_res_vector.csv'));

        save_cluster_csv( ...
            pos_clusters_diff(1:x), pos_cluster_p_diff(1:x), ...
            fullfile(table_dir_bin_size, 'positive_clusters_folded_theta_20_diff_res_vector.csv'));

        save_cluster_csv( ...
            neg_clusters_diff(1:x), neg_cluster_p_diff(1:x), ...
            fullfile(table_dir_bin_size, 'negative_clusters_folded_theta_20_diff_res_vector.csv'));

    end
    % save_cluster_csv(pos_clusters, pos_cluster_p, fullfile(table_dir_bin_size, 'positive_clusters_folded.csv'));
    % save_cluster_csv(neg_clusters, neg_cluster_p, fullfile(table_dir_bin_size, 'negative_clusters_folded.csv'));
end
end

function clusters = find_clusters(mask, t)
    % Each cluster will contain its indices and mass
    clusters = struct('indices', {}, 'mass', {});

    current_indices = [];

    for idx = 1:numel(mask)

        if mask(idx)
            % This point belongs to the current cluster
            current_indices(end + 1) = idx;
        end

        % Save when we reach a false point or the end of the array
        if (~mask(idx) || idx == numel(mask)) && ~isempty(current_indices)

            cluster_number = numel(clusters) + 1;

            clusters(cluster_number).indices = current_indices;
            clusters(cluster_number).mass = abs(sum(t(current_indices)));

            % Reset, ready for the next cluster
            current_indices = [];
        end
    end
end

function save_cluster_csv(clusters, cluster_p, filename)

    results = table('Size', [0, 7], ...
        'VariableTypes', ...
        {'double','double','string','string', ...
         'double','double','double'}, ...
        'VariableNames', ...
        {'Lag','Cluster','Indices','DeltaDegrees', ...
         'NumPoints','Mass','P_values'});

    delta_values = (0:90)';

    for lag = 1:numel(clusters)
        for c = 1:numel(clusters{lag})

            indices = clusters{lag}(c).indices;
            angles = delta_values(indices);
            p_tail = cluster_p{lag}(c);

            results(end + 1, :) = { ...
                lag, ...
                c, ...
                string(jsonencode(indices(:)')), ...
                string(jsonencode(angles(:)')), ...
                numel(indices), ...
                clusters{lag}(c).mass, ...
                p_tail
            };
        end
    end

    writetable(results, filename, 'Delimiter', ';');
end