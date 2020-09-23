clear all;

% Path vars
PATH_EEGLAB           = 'add_a_path_here';
PATH_AUTOCLEANED      = 'add_a_path_here';
PATH_TFDECOMP         = 'add_a_path_here';
PATH_PLOT             = 'add_a_path_here';
PATH_CORTISOL         = 'add_a_path_here';
PATH_VEUSZ            = 'add_a_path_here';

% ======================= SUBJECTS =========================================================================================================

% Define subjects (N=32)
subject_list = [1, 2, 4, 5, 6, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 23, 24, 25, 26, 27, 28];

% Determine and drop subjects with too few trials (<100)
preprostats = dlmread([PATH_AUTOCLEANED 'preprostats.csv']);
todrop = unique(preprostats(preprostats(:, 6) < 100, 1));
subject_list = setdiff(subject_list, todrop);

% ======================= OPTIONS =========================================================================================================

% Switch parts of the script on/off
to_execute = {'part3'};

% ============================ Part 1: Calculate lateralization index ============================================================================
if ismember('part1', to_execute)

    % Init eeglab
    addpath(PATH_EEGLAB);
    eeglab;

    % Get tf params
    tf_times = dlmread([PATH_TFDECOMP 'tf_times.csv']); 
    tf_freqs = dlmread([PATH_TFDECOMP 'tf_freqs.csv']);

    % Define lateralized channel pairs
    chanpairs = {[1, 2],... 
                 [3, 7], [4, 6],... 
                 [28, 32], [8, 11], [9, 10],...
                 [12, 16], [13, 15],...
                 [17, 22], [18, 21], [19, 20],...
                 [23, 27], [24, 26],...
                 [29, 31]};
    
    % Init latidx matrices
    latidx_cold      = zeros(length(subject_list), numel(chanpairs), length(tf_freqs), length(tf_times));
    latidx_warm      = zeros(length(subject_list), numel(chanpairs), length(tf_freqs), length(tf_times));
    latidx_cold_fake = zeros(length(subject_list), numel(chanpairs), length(tf_freqs), length(tf_times));
    latidx_warm_fake = zeros(length(subject_list), numel(chanpairs), length(tf_freqs), length(tf_times));

    % Iterate subjects
    for s = 1 : length(subject_list)

        % Talk
        fprintf('\nCalculating lateralization indices subject %i/%i\n', s, length(subject_list));

        % Load meta
        meta = dlmread([PATH_TFDECOMP num2str(subject_list(s)) '_powcube_meta.csv']); 

        % Iterate chanpairs
        for cp = 1 : numel(chanpairs)

            % Select chanpair
            chanpair = chanpairs{cp};

            % Load left and right channel pows
            load([PATH_TFDECOMP num2str(subject_list(s)) '_powcube_chan_' num2str(chanpair(1))]);
            pow_l = powcube;
            load([PATH_TFDECOMP num2str(subject_list(s)) '_powcube_chan_' num2str(chanpair(2))]);
            pow_r = powcube;

            % Get tf data by trials sorted by condition
            cold_ipsi_left    = pow_l(:, :, meta(:, 4) == 1 & meta(:, 3) == 1);
            cold_ipsi_right   = pow_r(:, :, meta(:, 4) == 2 & meta(:, 3) == 1);
            cold_contra_left  = pow_l(:, :, meta(:, 4) == 2 & meta(:, 3) == 1);
            cold_contra_right = pow_r(:, :, meta(:, 4) == 1 & meta(:, 3) == 1);
            warm_ipsi_left    = pow_l(:, :, meta(:, 4) == 1 & meta(:, 3) == 2);
            warm_ipsi_right   = pow_r(:, :, meta(:, 4) == 2 & meta(:, 3) == 2);
            warm_contra_left  = pow_l(:, :, meta(:, 4) == 2 & meta(:, 3) == 2);
            warm_contra_right = pow_r(:, :, meta(:, 4) == 1 & meta(:, 3) == 2);

            % Permute cold left
            idx_ipsi_swap = randsample(size(cold_ipsi_left, 3), round(size(cold_ipsi_left, 3) / 2));
            idx_ipsi_stay = setdiff([1 : size(cold_ipsi_left, 3)], idx_ipsi_swap);
            idx_contra_swap = randsample(size(cold_contra_left, 3), round(size(cold_contra_left, 3) / 2));
            idx_contra_stay = setdiff([1 : size(cold_contra_left, 3)], idx_contra_swap);
            fake_cold_ipsi_left   = squeeze(mean(cat(3, cold_ipsi_left(:, :, idx_ipsi_stay), cold_contra_left(:, :, idx_contra_swap)), 3));
            fake_cold_contra_left = squeeze(mean(cat(3, cold_contra_left(:, :, idx_contra_stay), cold_ipsi_left(:, :, idx_ipsi_swap)), 3));
    
            % Permute cold right
            idx_ipsi_swap = randsample(size(cold_ipsi_right, 3), round(size(cold_ipsi_right, 3) / 2));
            idx_ipsi_stay = setdiff([1 : size(cold_ipsi_right, 3)], idx_ipsi_swap);
            idx_contra_swap = randsample(size(cold_contra_right, 3), round(size(cold_contra_right, 3) / 2));
            idx_contra_stay = setdiff([1 : size(cold_contra_right, 3)], idx_contra_swap);
            fake_cold_ipsi_right   = squeeze(mean(cat(3, cold_ipsi_right(:, :, idx_ipsi_stay), cold_contra_right(:, :, idx_contra_swap)), 3));
            fake_cold_contra_right = squeeze(mean(cat(3, cold_contra_right(:, :, idx_contra_stay), cold_ipsi_right(:, :, idx_ipsi_swap)), 3));
    
            % Permute warm left
            idx_ipsi_swap = randsample(size(warm_ipsi_left, 3), round(size(warm_ipsi_left, 3) / 2));
            idx_ipsi_stay = setdiff([1 : size(warm_ipsi_left, 3)], idx_ipsi_swap);
            idx_contra_swap = randsample(size(warm_contra_left, 3), round(size(warm_contra_left, 3) / 2));
            idx_contra_stay = setdiff([1 : size(warm_contra_left, 3)], idx_contra_swap);
            fake_warm_ipsi_left   = squeeze(mean(cat(3, warm_ipsi_left(:, :, idx_ipsi_stay), warm_contra_left(:, :, idx_contra_swap)), 3));
            fake_warm_contra_left = squeeze(mean(cat(3, warm_contra_left(:, :, idx_contra_stay), warm_ipsi_left(:, :, idx_ipsi_swap)), 3));
    
            % Permute warm right
            idx_ipsi_swap = randsample(size(warm_ipsi_right, 3), round(size(warm_ipsi_right, 3) / 2));
            idx_ipsi_stay = setdiff([1 : size(warm_ipsi_right, 3)], idx_ipsi_swap);
            idx_contra_swap = randsample(size(warm_contra_right, 3), round(size(warm_contra_right, 3) / 2));
            idx_contra_stay = setdiff([1 : size(warm_contra_right, 3)], idx_contra_swap);
            fake_warm_ipsi_right   = squeeze(mean(cat(3, warm_ipsi_right(:, :, idx_ipsi_stay), warm_contra_right(:, :, idx_contra_swap)), 3));
            fake_warm_contra_right = squeeze(mean(cat(3, warm_contra_right(:, :, idx_contra_stay), warm_ipsi_right(:, :, idx_ipsi_swap)), 3));

            % Average left and right channels
            cold_ipsi   = (squeeze(mean(cold_ipsi_left, 3)) + squeeze(mean(cold_ipsi_right, 3))) / 2;
            cold_contra = (squeeze(mean(cold_contra_left, 3)) + squeeze(mean(cold_contra_right, 3))) / 2;
            warm_ipsi   = (squeeze(mean(warm_ipsi_left, 3)) + squeeze(mean(warm_ipsi_right, 3))) / 2;
            warm_contra = (squeeze(mean(warm_contra_left, 3)) + squeeze(mean(warm_contra_right, 3))) / 2;

            cold_ipsi_fake   = (fake_cold_ipsi_left + fake_cold_ipsi_right) / 2;
            cold_contra_fake = (fake_cold_contra_left + fake_cold_contra_right) / 2;
            warm_ipsi_fake   = (fake_warm_ipsi_left + fake_warm_ipsi_right) / 2;
            warm_contra_fake = (fake_warm_contra_left + fake_warm_contra_right) / 2;

            % Calculate lateralization indices
            latidx_cold(s, cp, :, :)      = (cold_ipsi - cold_contra) ./ (cold_ipsi + cold_contra);
            latidx_warm(s, cp, :, :)      = (warm_ipsi - warm_contra) ./ (warm_ipsi + warm_contra);
            latidx_cold_fake(s, cp, :, :) = (cold_ipsi_fake - cold_contra_fake) ./ (cold_ipsi_fake + cold_contra_fake);
            latidx_warm_fake(s, cp, :, :) = (warm_ipsi_fake - warm_contra_fake) ./ (warm_ipsi_fake + warm_contra_fake);

        end % End chanpair iteration
    end % End subject iteration

    % Save latidx
    save([PATH_TFDECOMP 'latidx_cold'], 'latidx_cold');
    save([PATH_TFDECOMP 'latidx_warm'], 'latidx_warm');
    save([PATH_TFDECOMP 'latidx_cold_fake'], 'latidx_cold_fake');
    save([PATH_TFDECOMP 'latidx_warm_fake'], 'latidx_warm_fake');

end % End part1

% ============================ Part 2: Plot topographies of lateralization ============================================================================
if ismember('part2', to_execute)

    % Init eeglab
    addpath(PATH_EEGLAB);
    eeglab;

    % Get chanlocs
    EEG = pop_loadset('filename', [num2str(subject_list(1)) '_autocleaned.set'], 'filepath', PATH_AUTOCLEANED, 'loadmode', 'info');
    chanlocs = EEG.chanlocs;

    % Get tf params
    tf_times = dlmread([PATH_TFDECOMP 'tf_times.csv']); 
    tf_freqs = dlmread([PATH_TFDECOMP 'tf_freqs.csv']);

    % Define lateralized channel pairs
    chanpairs = {[1, 2],... 
                    [3, 7], [4, 6],... 
                    [28, 32], [8, 11], [9, 10],...
                    [12, 16], [13, 15],...
                    [17, 22], [18, 21], [19, 20],...
                    [23, 27], [24, 26],...
                    [29, 31]};

    % Load latidx
    load([PATH_TFDECOMP 'latidx_cold']);
    load([PATH_TFDECOMP 'latidx_warm']);

    % Iterate time windows
    counter = 0;
    for t = -400 : 400 : 2800
        counter = counter + 1;

        % Init topo
        topo_cold = zeros(numel(chanlocs), 1);
        topo_warm = zeros(numel(chanlocs), 1);

        % Iterate chanpairs
        for cp = 1 : numel(chanpairs)

            % Get average chanpair latidx
            lidx_cold = squeeze(mean(squeeze(latidx_cold(:, cp, :, :)), 1));
            lidx_warm = squeeze(mean(squeeze(latidx_warm(:, cp, :, :)), 1));

            % Get average for time window in alpha range for topo
            topo_cold([chanpairs{cp}(1), chanpairs{cp}(2)]) = mean(squeeze(mean(lidx_cold(tf_freqs >= 8 & tf_freqs <= 12, tf_times >= t & tf_times < t + 400), 1)), 2);
            topo_warm([chanpairs{cp}(1), chanpairs{cp}(2)]) = mean(squeeze(mean(lidx_warm(tf_freqs >= 8 & tf_freqs <= 12, tf_times >= t & tf_times < t + 400), 1)), 2);

        end

        clim = [-0.05, 0.05];


        figure('Visible', 'off'); clf;
        pd = topo_cold;
        topoplot(pd, chanlocs, 'plotrad', 0.7, 'intrad', 0.7, 'intsquare', 'on', 'conv', 'off', 'electrodes', 'on');
        colormap('jet');
        caxis(clim);
        saveas(gcf, [PATH_PLOT 'topo_latidx_cold' num2str(counter) '.png']);

        figure('Visible', 'off'); clf;
        pd = topo_warm;
        topoplot(pd, chanlocs, 'plotrad', 0.7, 'intrad', 0.7, 'intsquare', 'on', 'conv', 'off', 'electrodes', 'on');
        colormap('jet');
        caxis(clim);
        saveas(gcf, [PATH_PLOT 'topo_latidx_warm' num2str(counter) '.png']);
    end

end % End part2

% ============================ Part 3: Some statistics ============================================================================
if ismember('part3', to_execute)

    % Get tf params
    tf_times = dlmread([PATH_TFDECOMP 'tf_times.csv']); 
    tf_freqs = dlmread([PATH_TFDECOMP 'tf_freqs.csv']);

    % Define lateralized channel pairs
    chanpairs = {[1, 2],... 
                 [3, 7], [4, 6],... 
                 [28, 32], [8, 11], [9, 10],...
                 [12, 16], [13, 15],...
                 [17, 22], [18, 21], [19, 20],...
                 [23, 27], [24, 26],...
                 [29, 31]};

    % Load latidx
    load([PATH_TFDECOMP 'latidx_cold']);
    load([PATH_TFDECOMP 'latidx_warm']);
    load([PATH_TFDECOMP 'latidx_cold_fake']);
    load([PATH_TFDECOMP 'latidx_warm_fake']);

    % Define channelpair for cluster
    idx = [10, 12, 13]; 

    % Average across channel pairs
    latcold = squeeze(mean(latidx_cold(:, idx, :, :), 2));
    latwarm = squeeze(mean(latidx_warm(:, idx, :, :), 2));
    latcold_fake = squeeze(mean(latidx_cold_fake(:, idx, :, :), 2));
    latwarm_fake = squeeze(mean(latidx_warm_fake(:, idx, :, :), 2));

    % Smooth data in time
    winlength = 50;
    if winlength
        latcold = movmean(latcold, winlength, 3);
        latwarm = movmean(latwarm, winlength, 3);
        latcold_fake = movmean(latcold_fake, winlength, 3);
        latwarm_fake = movmean(latwarm_fake, winlength, 3);
    end

    % Prune
    prune_idx_time = tf_times >= -500 & tf_times <= 3200;
    prune_idx_freq = tf_freqs >= 0 & tf_freqs <= 20;
    tf_times = tf_times(prune_idx_time);
    tf_freqs = tf_freqs(prune_idx_freq);
    latcold = latcold(:, prune_idx_freq, prune_idx_time);
    latwarm = latwarm(:, prune_idx_freq, prune_idx_time);
    latcold_fake = latcold_fake(:, prune_idx_freq, prune_idx_time);
    latwarm_fake = latwarm_fake(:, prune_idx_freq, prune_idx_time);

    % Calculate average of warm & cold
    latgeneral = (latwarm + latcold) / 2;
    latgeneral_fake = (latwarm_fake + latcold_fake) / 2;

    % Cluster permutation tests
    [latgen.sig_flag, latgen.ave, latgen.ave_fake, latgen.outline, latgen.apes, latgen.clust_sumt, latgen.clust_pvals, latgen.clust_apes, latgen.time_limits, latgen.freq_limits, latgen.cluster_idx]...
     = cluststats_2d_data(latgeneral, latgeneral_fake, tf_times, tf_freqs, 'pval_voxel', 0.05, 'tail', 1);
    [cold_vs_warm.sig_flag, cold_vs_warm.ave_cold, cold_vs_warm.ave_warm, cold_vs_warm.outline, cold_vs_warm.apes,...
     cold_vs_warm.clust_sumt, cold_vs_warm.clust_pvals, cold_vs_warm.clust_apes, cold_vs_warm.time_limits, cold_vs_warm.freq_limits, cold_vs_warm.cluster_idx]...
     = cluststats_2d_data(latcold, latwarm, tf_times, tf_freqs, 'pval_voxel', 0.05);

    % Save matrices for plotting
    dlmwrite([PATH_VEUSZ, 'latidx_general_ave.csv'], latgen.ave);
    dlmwrite([PATH_VEUSZ, 'latidx_general_contour.csv'], latgen.outline);
    dlmwrite([PATH_VEUSZ, 'latidx_general_apes.csv'], latgen.apes);

    dlmwrite([PATH_VEUSZ, 'latidx_ave_cold.csv'], cold_vs_warm.ave_cold);
    dlmwrite([PATH_VEUSZ, 'latidx_ave_warm.csv'], cold_vs_warm.ave_warm);
    dlmwrite([PATH_VEUSZ, 'latidx_ave_contour.csv'], cold_vs_warm.outline);
    dlmwrite([PATH_VEUSZ, 'latidx_ave_apes.csv'], cold_vs_warm.apes);

    % Plot cluster effects
    figure('Visible', 'off'); clf;

    cmap = 'jet';
    clim = [-0.05, 0.05];

    subplot(2, 1, 1)
    pd = latgen.ave;
    contourf(tf_times, tf_freqs, pd, 40, 'linecolor','none')
    hold on
    contour(tf_times, tf_freqs, latgen.outline, 1, 'linecolor', 'k', 'LineWidth', 2)
    colormap(cmap)
    set(gca, 'clim', clim, 'YScale', 'lin', 'YTick', [4, 8, 12, 20])
    colorbar;
    title('latidx')

    subplot(2, 1, 2)
    pd = latgen.apes;
    contourf(tf_times, tf_freqs, pd, 40, 'linecolor','none')
    colormap(cmap)
    set(gca, 'clim', [-0.5, 0.5], 'YScale', 'lin', 'YTick', [4, 8, 12, 20])
    colorbar;
    title('apes')

    saveas(gcf, [PATH_PLOT 'main_effect.png']);

    % Plot cluster effects
    figure('Visible', 'off'); clf;

    cmap = 'jet';
    clim = [-0.05, 0.05];

    subplot(3, 1, 1)
    pd = cold_vs_warm.ave_warm;
    contourf(tf_times, tf_freqs, pd, 40, 'linecolor','none')
    hold on
    contour(tf_times, tf_freqs, cold_vs_warm.outline, 1, 'linecolor', 'k', 'LineWidth', 2)
    colormap(cmap)
    set(gca, 'clim', clim, 'YScale', 'lin', 'YTick', [4, 8, 12, 20])
    colorbar;
    title('latidx warm')

    subplot(3, 1, 2)
    pd = cold_vs_warm.ave_cold;
    contourf(tf_times, tf_freqs, pd, 40, 'linecolor','none')
    hold on
    contour(tf_times, tf_freqs, cold_vs_warm.outline, 1, 'linecolor', 'k', 'LineWidth', 2)
    colormap(cmap)
    set(gca, 'clim', clim, 'YScale', 'lin', 'YTick', [4, 8, 12, 20])
    colorbar;
    title('latidx cold')

    subplot(3, 1, 3)
    pd = cold_vs_warm.apes;
    contourf(tf_times, tf_freqs, pd, 40, 'linecolor','none')
    colormap(cmap)
    set(gca, 'clim', [-0.5, 0.5], 'YScale', 'lin', 'YTick', [4, 8, 12, 20])
    colorbar;
    title('apes')

    saveas(gcf, [PATH_PLOT 'cold_vs_warm.png']);

    % Read cortisol data: VP | condition | C_delta_max | AUCi | AUCg
    cortisol_data = dlmread([PATH_CORTISOL, 'Cort.csv'], ';');

    % Cortisol measure labels
    cort_labels = {'C delta max', 'AUCi', 'AUCg'};

    % Get session cortisol 
    session_corts = [];
    for cm = 1 : 3
        for s = 1 : length(subject_list)
            session_corts(cm, s, 1) = cortisol_data(cortisol_data(:, 1) == subject_list(s) & cortisol_data(:, 2) == 1, 2 + cm); % stress
            session_corts(cm, s, 2) = cortisol_data(cortisol_data(:, 1) == subject_list(s) & cortisol_data(:, 2) == 2, 2 + cm); % control
        end
    end
    
    % Choose a measure
    cm = 3;

    % Determine cortisol increase as stress minus control
    cortisol_increase = (squeeze(session_corts(cm, :, 1)) - squeeze(session_corts(cm, :, 2)))';
    cortisol_stress = squeeze(session_corts(cm, :, 1))';

    % Calculate correlation of latdiff in cluster and cortisol increase
    lat_clust_diffs = [];
    lat_clust_stress = [];
    alpha_clust_time = [];
    latdiff = latcold - latwarm;
    for s = 1 : size(latdiff, 1)
        tmp = squeeze(latdiff(s, :, :));
        lat_clust_diffs(end + 1) = mean(tmp(cold_vs_warm.cluster_idx{1}));
        tmp = squeeze(latcold(s, :, :));
        lat_clust_stress(end + 1) = mean(tmp(cold_vs_warm.cluster_idx{1}));

        tmp = squeeze(latdiff(s, :, :));
        tmp = tmp(tf_freqs >= 8 & tf_freqs <= 12, tf_times >= 1100 & tf_times <= 1300);
        alpha_clust_time(end + 1) = mean(tmp(:));
    end

    cortisol_stress = log10(cortisol_stress);

    figure
    [r,p] = corrcoef(cortisol_stress, lat_clust_diffs');
    plot(cortisol_stress, lat_clust_diffs', 'k*')
    lsline

    % Lateralizations to 2d
    latdiff_2d = reshape(latdiff, size(latdiff, 1), size(latdiff, 2) * size(latdiff, 3));

    % Correlate
    session_diff_corrs = reshape(corr(cortisol_increase, latdiff_2d), size(latdiff, 2), size(latdiff, 3));
    dlmwrite([PATH_VEUSZ, 'correlation_diffs.csv'], session_diff_corrs);

    % Plot
    figure('Visible', 'off'); clf;
    cmap = 'jet';
    clim = [-0.5, 0.5];
    pd = session_diff_corrs;
    contourf(tf_times, tf_freqs, pd, 40, 'linecolor','none')
    colormap(cmap)
    set(gca, 'clim', clim, 'YScale', 'lin', 'YTick', [4, 8, 12, 20])
    colorbar;
    title('correlations')
    saveas(gcf, [PATH_PLOT 'session_difference_correlation.png']);

end % End part3