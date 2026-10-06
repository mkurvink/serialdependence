function estimate_polynomial_mvav_n_back()
%% estimating error scatter at iso, mid and ortho using a weighted flexible-degree polynomial fitting
scriptPath = mfilename('fullpath');
[currentDir,~,~] = fileparts(scriptPath);
analysisDir = fileparts(fileparts(currentDir));

tbl = readtable(fullfile(analysisDir,'results','SD_ma_master_table_n_back_first_half.csv'),'Delimiter',';');

% Initializing some variables
X = 6;
degrees = 2:3;
bin_size = 11;
error_variable = 'error_ori_deb';
nanunique = @(values) unique(values(~isnan(values)));

for x = 1:X
    delta_variable = sprintf('delta_%d',x);
    too_few_delta_count = 0;

    tbl_scatter = [];
    best_fit_all = [];
    best_polynomial_degree = [];
    mvav_scatter_aggregated = [];
    mvav_bias_aggregated = [];
    stimulus = [];
    winning_degree = [];

    % Loop over studies
    for i = 1:length(unique(tbl.studynum))
        tbl_i = tbl(tbl.studynum==i,:);
        tbl_i = tbl_i(abs(tbl_i.(delta_variable)) <= 90, :);  % Restrict to |delta| <= 90

        % fold if sd_deb is not removed already...
        if ~contains(error_variable,'sd_deb')
            err = tbl_i.(error_variable);
            err(tbl_i.(delta_variable)<0 & tbl_i.(delta_variable)~=-90) = -err(tbl_i.(delta_variable)<0 & tbl_i.(delta_variable)~=-90);
            tbl_i.(error_variable) = err;
        end

        % Loop over experiments and conditions
        nexps = max(tbl_i.expnum);
        for k = 1:nexps
            tbl_i_k = tbl_i(tbl_i.expnum==k,:);
            cond = unique(tbl_i_k.cond);
            ncond = numel(cond);

            for j = 1:ncond
                tbl_i_k_j = tbl_i_k(tbl_i_k.cond==cond(j),:);
                obs = unique(tbl_i_k_j.obs);
                nobs = numel(obs);

                for o = 1:nobs
                    tbl_i_k_j_o = tbl_i_k_j(tbl_i_k_j.obs==obs(o),:);

                    % Polynomial fitting
                    scatt = grpstats(tbl_i_k_j_o.(error_variable), abs(tbl_i_k_j_o.(delta_variable)), 'std');
                    biass  = grpstats(tbl_i_k_j_o.(error_variable), abs(tbl_i_k_j_o.(delta_variable)), 'mean');
                    [counts, ~] = groupcounts(abs(tbl_i_k_j_o.(delta_variable)));
                    delta_fit = nanunique(abs(tbl_i_k_j_o.(delta_variable)));
                    idx   = scatt ~= 0 & ~isnan(scatt);
                    biass = biass(idx);
                    scatt_fit = scatt(idx);
                    delta_fit = delta_fit(idx);
                    counts = counts(idx);

                    if numel(delta_fit) < max(degrees) + 1
                        too_few_delta_count = too_few_delta_count + 1;
                        continue
                    end

                    weights = counts./max(counts) ;

                    % BIC calculation
                    BIC_values = [];
                    for degree = degrees
                        fitOptions = fitoptions('Method', 'LinearLeastSquares', 'Weights', weights);
                        model = fit(delta_fit, scatt_fit, ['poly' num2str(degree)], fitOptions);

                        % Calculate residuals and BIC
                        y_fit = feval(model, delta_fit);
                        residuals = scatt_fit - y_fit;
                        SSR = sum((sqrt(weights) .* residuals).^2./sum(weights));
                        numParams = degree + 1;  % Polynomial degree + 1 for intercept
                        n = length(delta_fit);
                        BIC = n * log(SSR / n) + numParams * log(n);
                        BIC_values = [BIC_values; BIC];
                    end

                    % Select best polynomial degree using BIC
                    [sortedBIC, sortedIndex] = sort(BIC_values);
                    bestDegree = degrees(min(sortedIndex(sortedBIC < sortedBIC(1) + 2)));
                    winning_degree = [winning_degree, bestDegree];

                    best_model = fit(delta_fit, scatt_fit, ['poly' num2str(bestDegree)], fitOptions);
                    % Predict over range 0:90
                    best_fit = feval(best_model, [0:90]');

                    % Skip if scatter becomes negative or NaN
                    if any(best_fit([1, 46, 91]) <= 0 | any(isnan(best_fit([1, 46, 91]))))
                        continue
                    end

                    % Aggregated moving average bias
                    mvav = nan(91,1); mvav(delta_fit+1) = biass; mvav = [-mvav(end-1:-1:2);mvav]; mvav=repmat(mvav,3,1);
                    weights = nan(91,1); weights(delta_fit+1) = counts; weights = [weights(end-1:-1:2);weights]; weights=repmat(weights,3,1);
                    weighted_moving_avg = movsum(mvav.*weights,bin_size,'omitnan')./movsum(weights,bin_size,'omitnan');
                    weighted_moving_avg = weighted_moving_avg(180+[90:180]);
                    mvav_bias_aggregated = [mvav_bias_aggregated weighted_moving_avg];

                    % Aggregated moving average scatter
                    mvav = nan(91,1); mvav(delta_fit+1) = scatt_fit; mvav = [mvav(end-1:-1:2);mvav]; mvav=repmat(mvav,3,1);
                    weights = nan(91,1); weights(delta_fit+1) = counts; weights = [weights(end-1:-1:2);weights]; weights=repmat(weights,3,1);
                    weighted_moving_avg = movsum(mvav.*weights,bin_size,'omitnan')./movsum(weights,bin_size,'omitnan');
                    weighted_moving_avg = weighted_moving_avg(180+[90:180]);
                    mvav_scatter_aggregated = [mvav_scatter_aggregated weighted_moving_avg];

                    % Save results for this level for plotting
                    best_polynomial_degree = [best_polynomial_degree; bestDegree];
                    best_fit_all = [best_fit_all, best_fit];
                    stimulus = [stimulus, strcmpi(tbl_i_k_j_o.stimulus{1},'Motion')];

                    % Store polynomial fit results
                    tmp = repmat(tbl_i_k_j_o(1,ismember(tbl_i_k_j_o.Properties.VariableNames,{'obsid','codenum'})),3,1);
                    tmp.SI = categorical({'iso' 'mid' 'ortho'})';

                    absolute_delta = abs(tbl_i_k_j_o.(delta_variable));
                    current_bin = nan(height(tbl_i_k_j_o),1);

                    current_bin(absolute_delta<=10) = 1;
                    current_bin(absolute_delta>=40 & absolute_delta<=50) = 2;
                    current_bin(absolute_delta>=80 & absolute_delta<=90) = 3;

                    bin_scatter = nan(3,1);

                    for bin_number = 1:3
                        bin_errors = tbl_i_k_j_o.(error_variable)(current_bin==bin_number);
                        bin_scatter(bin_number) = std(bin_errors,'omitnan');
                    end

                    tmp.ES = best_fit([1,46,91]);
                    tmp.bin_scatter = bin_scatter;
                    tmp.ntrials = repmat(length(counts),3,1);
                    tbl_scatter = [tbl_scatter; tmp];
                end
            end
        end
    end

    fprintf('%s: skipped %d participant-condition combinations with fewer than %d usable delta values.\n',delta_variable,too_few_delta_count,max(degrees)+1);
    output_name = sprintf('sd_estimates_first_half_delta_%d.mat',x);
    output_file = fullfile(analysisDir,'results',output_name);
    save(output_file,'tbl_scatter','best_fit_all','mvav_bias_aggregated','mvav_scatter_aggregated','stimulus','winning_degree');

end
end