function [pos_clusters, neg_clusters, pos_cluster_p, neg_cluster_p] = estimate_bias_n_back_cluster_based()
%% estimating error scatter at iso, mid and ortho using a weighted flexible-degree polynomial fitting
scriptPath = mfilename('fullpath');
[currentDir,~,~] = fileparts(scriptPath);
analysisDir = fileparts(fileparts(currentDir));
output_dir = fullfile(analysisDir, 'results');
table_dir = fullfile(output_dir, 'tables_first_half');
tbl = readtable(fullfile(analysisDir,'results','SD_ma_master_table_n_back_first_half.csv'),'Delimiter',';');

% Initializing some variables
X = 10;
error_variable = 'error_ori_deb';
alpha = 0.05;
bin_sizes = [3, 11];        %%%%%%%%%%%%%%%%%
n_permutations = 10000;     %%%%%%%%%%%%%%%%%

group_bias = cell(1, X);
t_values = nan(180, X);
pos_mask = false(180, X);
neg_mask = false(180, X);

pos_clusters = cell(1, X);
neg_clusters = cell(1, X);

pos_counts = cell(1, X);
neg_counts = cell(1, X);

pos_cluster_p = cell(1, X);
neg_cluster_p = cell(1, X);

% Computing bias for true delta values
for bin_size = bin_sizes
    table_dir_bin_size = fullfile(table_dir, sprintf('%d_degree_window', bin_size));
    if ~isfolder(table_dir_bin_size)
        mkdir(table_dir_bin_size);
    end
    for x = 1:X
        delta_variable = sprintf('delta_%d',x);
        mvav_bias_aggregated = [];

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

                        signed_delta = tbl_i_k_j_o.(delta_variable);
                        signed_error = tbl_i_k_j_o.(error_variable);

                        valid = ~isnan(signed_delta) & ~isnan(signed_error);
                        signed_delta = signed_delta(valid);
                        signed_error = signed_error(valid);

                        if isempty(signed_delta)
                            continue
                        end

                        % Calculate mean signed error and trial count per delta                   
                        biass = grpstats(signed_error,signed_delta,'mean');
                        [counts, ~] = groupcounts(signed_delta);
                        delta_fit = unique(signed_delta);

                        % Aggregated moving average bias
                        mvav = nan(180,1); mvav(delta_fit+90) = biass; mvav=repmat(mvav,3,1);
                        weights = nan(180,1); weights(delta_fit+90) = counts; weights=repmat(weights,3,1);
                        weighted_moving_avg = movsum(mvav.*weights,bin_size,'omitnan')./movsum(weights,bin_size,'omitnan');
                        weighted_moving_avg = weighted_moving_avg(180+[1:180]);
                        mvav_bias_aggregated = [mvav_bias_aggregated weighted_moving_avg];

                    end
                end
            end
        end

        group_bias{x} = mean(mvav_bias_aggregated,2,'omitnan');
        group_delta = (-89:90)';

        % One test per row, across participant columns, against zero
        [~, ~, ~, stats] = ttest(mvav_bias_aggregated, 0, 'Dim', 2);

        t_values(:, x) = stats.tstat;
        t_threshold = tinv(1 - alpha/2, stats.df);

        pos_mask(:, x) = stats.tstat > t_threshold;
        neg_mask(:, x) = stats.tstat < -t_threshold;

        pos_clusters{x} = find_clusters(pos_mask(:, x), t_values(:, x));
        neg_clusters{x} = find_clusters(neg_mask(:, x), t_values(:, x));

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

        bias_table = table(group_delta, group_bias{x}, 'VariableNames', {'DeltaDegrees', 'MeanBias'});
        filename = fullfile(table_dir_bin_size, sprintf('group_bias_N_%d.csv', x));
        writetable(bias_table, filename, 'Delimiter', ';');

    end

    % Computing bias values for shuffled delta values
    for x = 1:X
        delta_variable = sprintf('delta_%d',x);
        lag_timer = tic;
        % These are the clusters from the original, unshuffled data
        observed_pos_mass = [pos_clusters{x}.mass];
        observed_neg_mass = [neg_clusters{x}.mass];

        pos_counts{x} = zeros(size(observed_pos_mass));
        neg_counts{x} = zeros(size(observed_neg_mass));

        fprintf('\nStarting %d permutations for N-%d...\n',n_permutations,x);

        % Permutation test on shuffled data
        for permutation = 1:n_permutations
            mvav_bias_aggregated = [];
            t_values_perm = nan(180, X);
            pos_mask_perm = false(180, X);
            neg_mask_perm = false(180, X);

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

                            signed_delta = tbl_i_k_j_o.(delta_variable);
                            signed_error = tbl_i_k_j_o.(error_variable);

                            valid = ~isnan(signed_delta) & ~isnan(signed_error);
                            signed_delta = signed_delta(valid);
                            signed_error = signed_error(valid);

                            if isempty(signed_delta)
                                continue
                            end

                            shuffled_delta = signed_delta(randperm(numel(signed_delta)));

                            % Calculate mean signed error and count per delta for the shuffled deltas
                            biass = grpstats(signed_error,shuffled_delta,'mean');
                            [counts, ~] = groupcounts(shuffled_delta);
                            delta_fit = unique(shuffled_delta);

                            % Aggregated moving average bias
                            mvav = nan(180,1); mvav(delta_fit+90) = biass; mvav=repmat(mvav,3,1);
                            weights = nan(180,1); weights(delta_fit+90) = counts; weights=repmat(weights,3,1);
                            weighted_moving_avg = movsum(mvav.*weights,bin_size,'omitnan')./movsum(weights,bin_size,'omitnan');
                            weighted_moving_avg = weighted_moving_avg(180+[1:180]);
                            mvav_bias_aggregated = [mvav_bias_aggregated weighted_moving_avg];

                        end
                    end
                end
            end

            % One test per row, across participant columns, against zero
            [~, ~, ~, stats] = ttest(mvav_bias_aggregated, 0, 'Dim', 2);

            t_values_perm(:, x) = stats.tstat;

            t_threshold = tinv(1 - alpha/2, stats.df);

            pos_mask_perm(:, x) = stats.tstat > t_threshold;
            neg_mask_perm(:, x) = stats.tstat < -t_threshold;

            pos_clusters_perm{x} = find_clusters(pos_mask_perm(:, x), t_values_perm(:, x));
            neg_clusters_perm{x} = find_clusters(neg_mask_perm(:, x), t_values_perm(:, x));

            % Largest mass in this permutation; zero if no clusters exist
            max_mass = max([0, pos_clusters_perm{x}.mass, neg_clusters_perm{x}.mass]);
            
            % Update one count for each observed cluster
            pos_counts{x} = pos_counts{x} + (max_mass >= observed_pos_mass);
            neg_counts{x} = neg_counts{x} + (max_mass >= observed_neg_mass);

            if mod(permutation,1000) == 0
                fprintf('N-%d: %d/%d permutations completed (%.1f minutes)\n',x,permutation,n_permutations,toc(lag_timer)/60);
            end
        end

        pos_cluster_p{x} = (pos_counts{x} + 1) / (n_permutations + 1);
        neg_cluster_p{x} = (neg_counts{x} + 1) / (n_permutations + 1);

    end
    save_cluster_csv(pos_clusters, pos_cluster_p, fullfile(table_dir_bin_size, 'positive_clusters.csv'));
    save_cluster_csv(neg_clusters, neg_cluster_p, fullfile(table_dir_bin_size, 'negative_clusters.csv'));
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

    delta_values = (-89:90)';

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