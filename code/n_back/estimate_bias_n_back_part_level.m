function [amplitudes,amp_std_delta] = estimate_bias_n_back_part_level()
%% estimating error scatter at iso, mid and ortho using a weighted flexible-degree polynomial fitting
scriptPath = mfilename('fullpath');
[currentDir,~,~] = fileparts(scriptPath);
analysisDir = fileparts(fileparts(currentDir));
resultsDir = fullfile(analysisDir,'results');
p_value_output_file = fullfile(resultsDir,'fit_p_values_n_back_first_half.csv');

tbl = readtable(fullfile(analysisDir,'results','SD_ma_master_table_n_back_first_half.csv'),'Delimiter',';');

% Initializing some variables
X = 6;
error_variable = 'error_ori_deb';
c = sqrt(2) / exp(-0.5);
amp_std_delta = nan(X,2);
amplitudes = cell(1,X);
fit_p_values = cell(1,X);
signed_p_values = nan(1,X);
absolute_p_values = nan(1,X);
n_permutations = 10000;  % Number of permutations for the permutation test

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
    lag_amplitudes = [];
    count_exc = 0;
    count_x_data = 0;
    studies = unique(tbl.studynum);

    for i = 1:numel(studies)
        study_number = studies(i);
        tbl_i = tbl(tbl.studynum==study_number,:);
        tbl_i = tbl_i(tbl_i.(delta_variable)>-90 & tbl_i.(delta_variable)<=90,:);

        % Loop over experiments and conditions
        experiments = unique(tbl_i.expnum);

        for k = 1:numel(experiments)
            experiment_number = experiments(k);
            tbl_i_k = tbl_i(tbl_i.expnum==experiment_number,:);
            cond = unique(tbl_i_k.cond);

            for j = 1:numel(cond)
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

                    count_x_data = count_x_data + 1;
                    if isempty(signed_delta)
                        count_exc = count_exc + 1;
                        continue
                    end

                    x_data = signed_delta;
                    y_data = signed_error;

                    if numel(unique(x_data)) < 5
                        count_exc = count_exc + 1;
                        continue
                    end

                    % Fit model
                    fitted_parameters = lsqcurvefit(DoG,starting_values,x_data,y_data,lower_bounds,upper_bounds,options);

                    % Save observer-specific amplitudes
                    lag_amplitudes(end+1) = fitted_parameters(1);

                end
            end
        end
    end
    
    fprintf('%d out of %d participant-level datasets were excluded for delta_%d (%.2f%%).\n',count_exc,count_x_data,x,100*count_exc/count_x_data);

    % Save X-specific amplitudes
    amplitudes{x} = lag_amplitudes;
    amp_std_delta(x,:) = [mean(lag_amplitudes,'omitnan'),std(lag_amplitudes,'omitnan')];

end


% Computing bias values for shuffled delta values
for x = 1:X
    delta_variable = sprintf('delta_%d',x);
    lag_timer = tic;

    observed_amplitudes = amplitudes{x}(:);
    n_observed_fits = numel(observed_amplitudes);

    fit_perm_higher_counts = zeros(n_observed_fits,1);
    permuted_signed_means = nan(n_permutations,1);
    permuted_absolute_means = nan(n_permutations,1);

    fprintf('\nStarting %d permutations for N-%d...\n',n_permutations,x);

    % Permutation test on shuffled data
    for permutation = 1:n_permutations
        lag_amplitudes = [];

        studies = unique(tbl.studynum);

        for i = 1:numel(studies)
            study_number = studies(i);
            tbl_i = tbl(tbl.studynum==study_number,:);
            tbl_i = tbl_i(tbl_i.(delta_variable)>-90 & tbl_i.(delta_variable)<=90,:);

            % Loop over experiments and conditions
            experiments = unique(tbl_i.expnum);

            for k = 1:numel(experiments)
                experiment_number = experiments(k);
                tbl_i_k = tbl_i(tbl_i.expnum==experiment_number,:);
                cond = unique(tbl_i_k.cond);

                for j = 1:numel(cond)
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

                        x_data = shuffled_delta;
                        y_data = signed_error;

                        if numel(unique(x_data)) < 5
                            continue
                        end

                        % Fit model
                        fitted_parameters = lsqcurvefit(DoG,starting_values,x_data,y_data,lower_bounds,upper_bounds,options);

                        % Save observer-specific amplitudes
                        lag_amplitudes(end+1) = fitted_parameters(1);

                    end
                end
            end
        end

        permuted_lag_amplitudes = lag_amplitudes(:);
        if numel(permuted_lag_amplitudes) ~= n_observed_fits
            error('Number of fitted amplitudes differs during permutation %d for N-%d.',permutation,x)
        end

        fit_perm_higher_counts = fit_perm_higher_counts + (abs(permuted_lag_amplitudes) >= abs(observed_amplitudes));   % Adds a 1 if the permuted amplitude is greater than or equal to the observed amplitude, resulting in f.e. [0, 534, 82, 25, 10, 5245,...]
        permuted_signed_means(permutation) = mean(permuted_lag_amplitudes,'omitnan');   % Takes the mean of the permuted amplitudes across all observers for this permutation
        permuted_absolute_means(permutation) = mean(abs(permuted_lag_amplitudes),'omitnan');    % Takes the mean of the absolute values of the permuted amplitudes across all observers for this permutation

        if mod(permutation,1000) == 0 % Change in relation to number of permutations for sueful progress reports
            fprintf('N-%d: %d/%d permutations completed (%.1f minutes)\n',x,permutation,n_permutations,toc(lag_timer)/60);
        end

    end

    % Calculate p-value per x
    fit_p_values{x} = (1 + fit_perm_higher_counts) ./ (n_permutations + 1);  % Essentially the count of when the permuted amplitude is greater than or equal to the observed amplitude, divided by the number of permutations + 1
    bonferroni_alpha = 0.05 / n_observed_fits;
    significant_bonferroni = fit_p_values{x} < bonferroni_alpha;

    positive = observed_amplitudes > 0;
    negative = observed_amplitudes < 0;

    positive_count = sum(positive);
    negative_count = sum(negative);
    significant_positive_count = sum(significant_bonferroni & positive);
    significant_negative_count = sum(significant_bonferroni & negative);

    fprintf('\nN-%d fit-level results:\n',x);
    fprintf('Positive amplitudes: %d/%d (%.2f%%)\n',positive_count,n_observed_fits,100*positive_count/n_observed_fits);
    fprintf('Negative amplitudes: %d/%d (%.2f%%)\n',negative_count,n_observed_fits,100*negative_count/n_observed_fits);
    fprintf('Significant positive amplitudes: %d/%d (%.2f%%)\n',significant_positive_count,n_observed_fits,100*significant_positive_count/n_observed_fits);
    fprintf('Significant negative amplitudes: %d/%d (%.2f%%)\n',significant_negative_count,n_observed_fits,100*significant_negative_count/n_observed_fits);

    observed_signed_mean = mean(observed_amplitudes,'omitnan');
    observed_absolute_mean = mean(abs(observed_amplitudes),'omitnan');

    signed_p_values(x) = (1 + sum(abs(permuted_signed_means) >= abs(observed_signed_mean))) / (n_permutations + 1);
    absolute_p_values(x) = (1 + sum(permuted_absolute_means >= observed_absolute_mean)) / (n_permutations + 1);

    fprintf('Signed mean amplitude p-value for N-%d: %.5f\n',x,signed_p_values(x));
    fprintf('Absolute mean amplitude p-value for N-%d: %.5f\n',x,absolute_p_values(x));

    writetable(fit_p_values,p_value_output_file,'Delimiter',';');   % Save the fit-level p-values to a CSV file to inspect

end
end


% addpath('/Users/Kurvi001/Documents/Serial_Dependence/Analysis/code/n_back')
% savepath
% [amplitudes,amp_std_delta] = estimate_bias_n_back_part_level();