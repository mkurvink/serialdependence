function figure3B_S7_n_back()

scriptPath = mfilename('fullpath');
[currentDir,~,~] = fileparts(scriptPath);
analysisDir = fileparts(fileparts(currentDir));

addpath(genpath(fullfile(analysisDir,'code','utils')))

resultsDir = fullfile(analysisDir,'results');
figuresDir = fullfile(resultsDir,'figures_first_half');

if ~exist(figuresDir,'dir')
    mkdir(figuresDir)
end

X = 6;  % Number of delta conditions
for x = 1:X
    %% Figure 3B: Mean bias moving average and mean error scatter + polynomial fit
    load(fullfile(resultsDir,sprintf('sd_estimates_first_half_delta_%d.mat',x)))    
    FontName = 'Arial';

    figure('Units','normalized','position',[.05 .4 .25 .7]);
    % prepare bias for plotting
    mvav_bias_aggregated_curve = mvav_bias_aggregated(:)';
    % prepare scatter for plotting
    mvav_scatter_aggregated_curve = mvav_scatter_aggregated(:)';
    aggregated_fit_curve = best_fit_all(:)';

    subplot(211);hold off;
    hold on;
    p=plot_mean_ci(repmat(0:90,1,length(aggregated_fit_curve)/91),mvav_bias_aggregated_curve,1,1,'k','mean',[]);
    p.LineWidth = 6;
    p.LineStyle = ':';

    subplot(212)
    hold on;
    p=plot_mean_ci(repmat(0:90,1,length(aggregated_fit_curve)/91),mvav_scatter_aggregated_curve,1,1,'k','mean',[]);
    p.LineWidth = 6;
    p.LineStyle = ':';
    p=plot_mean_ci(repmat(0:90,1,length(aggregated_fit_curve)/91),aggregated_fit_curve,1,1,'b','mean',[]);
    p.LineWidth = 6;
    legend(p,'Polynomial fit with 95% CI')

    subplot(211)
    % Configure plot
    xlabel('|Δ| (°)');
    ylabel('Bias (°)');
    title('Data');
    xticks([0 45 90]);
    set(gca, 'FontSize', 20, 'LineWidth', 1.9, 'FontName', FontName);
    xticklabels({'iso' 'mid' 'ortho'});

    subplot(212)
    % Configure plot
    xlabel('|Δ| (°)');
    ylabel('Error Scatter (°)');
    title('Data vs. Polynomial Fit');
    xticks([0 45 90]);
    set(gca, 'FontSize', 20, 'LineWidth', 1.9, 'FontName', FontName);
    xticklabels({'iso' 'mid' 'ortho'});
    
    set(gcf,'PaperOrientation','landscape')
    set(gcf, 'PaperUnits', 'normalized');
    exportgraphics(gcf,fullfile(figuresDir,sprintf('Figure3B_delta_%d.pdf',x)),'BackgroundColor','none')
    savefig(gcf,fullfile(figuresDir,sprintf('Figure3B_delta_%d.fig',x)));
    %% Figure S7: Mean bias moving average and mean error scatter + polynomial fit
    figure('Units','normalized','position',[.05 .4 .5 .7]);
    for i = 0:1
        % prepare bias for plotting
        mvav_bias_aggregated_curve = mvav_bias_aggregated(:,stimulus==i); mvav_bias_aggregated_curve=mvav_bias_aggregated_curve(:)';
        % prepare scatter for plotting
        mvav_scatter_aggregated_curve = mvav_scatter_aggregated(:,stimulus==i);mvav_scatter_aggregated_curve=mvav_scatter_aggregated_curve(:)';
        aggregated_fit_curve = best_fit_all(:,stimulus==i);aggregated_fit_curve = aggregated_fit_curve(:)';

        subplot(2,2,i+1);hold off;
        hold on;
        p=plot_mean_ci(repmat(0:90,1,length(aggregated_fit_curve)/91),mvav_bias_aggregated_curve,1,1,'k','mean',[]);
        p.LineWidth = 6;
        p.LineStyle = ':';

        subplot(2,2,i+3)
        hold on;
        p=plot_mean_ci(repmat(0:90,1,length(aggregated_fit_curve)/91),mvav_scatter_aggregated_curve,1,1,'k','mean',[]);
        p.LineWidth = 6;
        p.LineStyle = ':';
        p=plot_mean_ci(repmat(0:90,1,length(aggregated_fit_curve)/91),aggregated_fit_curve,1,1,'b','mean',[]);
        p.LineWidth = 6;


        subplot(2,2,i+1)
        % Configure plot
        xlabel('|Δ| (°)');
        ylabel('Bias (°)');
        title('Data');
        xticks([0 45 90]);
        set(gca, 'FontSize', 20, 'LineWidth', 1.9, 'FontName', FontName);
        xticklabels({'iso' 'mid' 'ortho'});
        if i == 0 
            title({'Orientation'})
        end
        if i==1
            title({'Motion'})
        end

        subplot(2,2,i+3)
        % Configure plot
        xlabel('|Δ| (°)');
        ylabel('Error Scatter (°)');
        xticks([0 45 90]);
        set(gca, 'FontSize', 20, 'LineWidth', 1.9, 'FontName', FontName);
        xticklabels({'iso' 'mid' 'ortho'});
        sgtitle('Empirical Data Split Based On Stimulus Type','FontName',FontName,'FontSize',22)

        set(gcf,'PaperOrientation','landscape')
        set(gcf, 'PaperUnits', 'normalized');
        exportgraphics(gcf,fullfile(figuresDir,sprintf('Figure3B_ori_vs_motion_delta_%d.pdf',x)),'BackgroundColor','none')
        savefig(gcf,fullfile(figuresDir,sprintf('Figure3B_ori_vs_motion_delta_%d.fig',x)));

    end
end
end