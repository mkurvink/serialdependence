function [amplitudes,p_values] = estimate_bias_n_back()
%% estimating error scatter at iso, mid and ortho using a weighted flexible-degree polynomial fitting
scriptPath = mfilename('fullpath');
[currentDir,~,~] = fileparts(scriptPath);
analysisDir = fileparts(fileparts(currentDir));

tbl = readtable(fullfile(analysisDir,'results','SD_ma_master_table_n_back_first_half.csv'),'Delimiter',';');

% Initializing some variables
X = 4;
error_variable = 'error_ori_deb';
c = sqrt(2) / exp(-0.5);
amplitudes = nan(1,X);
p_values = nan(1,X);
bin_size = 11;
n_permutations = 10000;

% Derivative of Gaussian curve
DoG = @(parameters,x) x .* parameters(1) .* (1./(sqrt(2).*parameters(2))) .* c .* exp(-(x./(sqrt(2).*parameters(2))).^2);

% Arbitrary specifications for the curve fitting 
% Because of the DoG formula construction the second input for the bounds represents degrees directly.
starting_values = [2,25];
lower_bounds = [-20,3];
upper_bounds = [20,85];
options = optimoptions('lsqcurvefit','Display','off');

% Computing bias for true delta values 
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

    group_bias = mean(mvav_bias_aggregated,2,'omitnan');
    group_delta = (-89:90)';

    x_data = group_delta;
    y_data = group_bias;

    valid = ~isnan(x_data) & ~isnan(y_data);
    x_data = x_data(valid);
    y_data = y_data(valid);

    % Fit model
    fitted_parameters = lsqcurvefit(DoG,starting_values,x_data,y_data,lower_bounds,upper_bounds,options);

    % Save X-specific amplitudes
    amplitudes(x) = fitted_parameters(1);
    fprintf('Observed N-%d amplitude: %.4f\n',x,amplitudes(x));
end


% Computing bias values for shuffled delta values
for x = 1:X
    delta_variable = sprintf('delta_%d',x);
    permuted_amplitudes = nan(n_permutations,1);
    lag_timer = tic;

    fprintf('\nStarting %d permutations for N-%d...\n',n_permutations,x);

    % Permutation test on shuffled data
    for permutation = 1:n_permutations
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

        group_bias = mean(mvav_bias_aggregated,2,'omitnan');
        group_delta = (-89:90)';

        x_data = group_delta;
        y_data = group_bias;

        valid = ~isnan(x_data) & ~isnan(y_data);
        x_data = x_data(valid);
        y_data = y_data(valid);

        % Fit model
        fitted_parameters = lsqcurvefit(DoG,starting_values,x_data,y_data,lower_bounds,upper_bounds,options);

        % Save permuted amplitudes
        permuted_amplitudes(permutation) = fitted_parameters(1);
        if mod(permutation,500) == 0
            fprintf('N-%d: %d/%d permutations completed (%.1f minutes)\n',x,permutation,n_permutations,toc(lag_timer)/60);
        end
    end

    % Calculate X-specific p-values by taking the absolute amplitudes and comparing to absolute permuted amplitudes
    p_values(x) = (1 + sum(abs(permuted_amplitudes) >= abs(amplitudes(x)))) / (n_permutations + 1);
    fprintf('Finished N-%d: amplitude = %.4f, p = %.5f, time = %.1f minutes\n',x,amplitudes(x),p_values(x),toc(lag_timer)/60);
end

results = table((1:X)',amplitudes',p_values','VariableNames',{'Lag','Amplitude','PValue'});
fprintf('\nFinal results:\n');
disp(results)
end
