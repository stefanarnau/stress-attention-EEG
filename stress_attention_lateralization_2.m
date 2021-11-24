clear all;

% Path vars
PATH_EEGLAB           = '/home/plkn/eeglab2021.1/';
PATH_FIELDTRIP        = '/home/plkn/fieldtrip-master/';
PATH_AUTOCLEANED      = '/mnt/data_heap/exp1027/eeg/2_autocleaned/';
PATH_TFDECOMP         = '/mnt/data_heap/exp1027/eeg/3_tfdecomp/';
PATH_PLOT             = '/mnt/data_heap/exp1027/plt_new/';

% ======================= SUBJECTS =========================================================================================================

% Define subjects (N=32)
subject_list = [1, 2, 4, 5, 6, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 23, 24, 25, 26, 27, 28];

% Determine and drop subjects with too few trials (<100)
preprostats = dlmread([PATH_AUTOCLEANED 'preprostats.csv']);
todrop = unique(preprostats(preprostats(:, 6) < 100, 1));
subject_list = setdiff(subject_list, todrop);

% ======================= OPTIONS =========================================================================================================

% Switch parts of the script on/off
to_execute = {'part2a'};

% ============================ Part 1: Calculate lateralization index ============================================================================
if ismember('part1', to_execute)

    % Init eeglab
    addpath(PATH_EEGLAB);
    eeglab;

    % Get tf params
    tf_times = dlmread([PATH_TFDECOMP 'tf_times.csv']); 
    tf_freqs = dlmread([PATH_TFDECOMP 'tf_freqs.csv']);

    % ERSP baseline latencies
    ersp_bl = [-500, -200];
    [~, blidx1] = min(abs(tf_times - ersp_bl(1)));
    [~, blidx2] = min(abs(tf_times - ersp_bl(2)));
    zersp_bl = [-200, 0];
    [~, zblidx1] = min(abs(tf_times - ersp_bl(1)));
    [~, zblidx2] = min(abs(tf_times - ersp_bl(2)));

    % Channel lateralizations
    EEG = pop_loadset('filename', [num2str(subject_list(1)) '_autocleaned.set'], 'filepath', PATH_AUTOCLEANED, 'loadmode', 'info');
    chanlats = cell2mat({EEG.chanlocs.Y});
    chanlats(abs(chanlats) < 1e-12) = 0; % center
    chanlats(chanlats > 0) = 1; % left
    chanlats(chanlats < 0) = 2; % right

    % ERSP matrices
    ersp_cold_ipsi = zeros(length(subject_list), 32, length(tf_freqs), length(tf_times));
    ersp_cold_contra = zeros(length(subject_list), 32, length(tf_freqs), length(tf_times));
    ersp_warm_ipsi = zeros(length(subject_list), 32, length(tf_freqs), length(tf_times));
    ersp_warm_contra = zeros(length(subject_list), 32, length(tf_freqs), length(tf_times));

    % Alpha matrices
    atra_cold_ipsi   = zeros(length(subject_list), 32, length(tf_times)); 
    atra_cold_contra = zeros(length(subject_list), 32, length(tf_times));
    atra_warm_ipsi   = zeros(length(subject_list), 32, length(tf_times)); 
    atra_warm_contra = zeros(length(subject_list), 32, length(tf_times));
    atra_db_cold_ipsi   = zeros(length(subject_list), 32, length(tf_times)); 
    atra_db_cold_contra = zeros(length(subject_list), 32, length(tf_times));
    atra_db_warm_ipsi   = zeros(length(subject_list), 32, length(tf_times)); 
    atra_db_warm_contra = zeros(length(subject_list), 32, length(tf_times));
    atra_z_cold_ipsi   = zeros(length(subject_list), 32, length(tf_times)); 
    atra_z_cold_contra = zeros(length(subject_list), 32, length(tf_times));
    atra_z_warm_ipsi   = zeros(length(subject_list), 32, length(tf_times)); 
    atra_z_warm_contra = zeros(length(subject_list), 32, length(tf_times));
    
    % Iterate subjects
    for s = 1 : length(subject_list)

        % Load meta
        meta = dlmread([PATH_TFDECOMP num2str(subject_list(s)) '_powcube_meta.csv']); 

        % Iterate channels
        for ch = 1 : 32

            fprintf('\n subject %i/%i, channel %i/32', s, length(subject_list), ch)

            % Load power
            load([PATH_TFDECOMP num2str(subject_list(s)) '_powcube_chan_' num2str(ch)]);

            % If left channel
            if chanlats(ch) == 1
                cold_ipsi   = powcube(:, :, meta(:, 4) == 1 & meta(:, 3) == 1);
                cold_contra = powcube(:, :, meta(:, 4) == 2 & meta(:, 3) == 1);
                warm_ipsi   = powcube(:, :, meta(:, 4) == 1 & meta(:, 3) == 2);
                warm_contra = powcube(:, :, meta(:, 4) == 2 & meta(:, 3) == 2);
            % If right channel
            elseif chanlats(ch) == 2
                cold_ipsi   = powcube(:, :, meta(:, 4) == 2 & meta(:, 3) == 1);
                cold_contra = powcube(:, :, meta(:, 4) == 1 & meta(:, 3) == 1);
                warm_ipsi   = powcube(:, :, meta(:, 4) == 2 & meta(:, 3) == 2);
                warm_contra = powcube(:, :, meta(:, 4) == 1 & meta(:, 3) == 2);
            elseif chanlats(ch) == 0
                cold_ipsi   = powcube(:, :, meta(:, 3) == 1);
                cold_contra = powcube(:, :, meta(:, 3) == 1);
                warm_ipsi   = powcube(:, :, meta(:, 3) == 2);
                warm_contra = powcube(:, :, meta(:, 3) == 2);
            end

            % Calculate ERSPs
            % ersp_cold_ipsi(s, ch, :, :)   = znorm(cold_ipsi, blidx1, blidx2);
            % ersp_cold_contra(s, ch, :, :) = znorm(cold_contra, blidx1, blidx2);
            % ersp_warm_ipsi(s, ch, :, :)   = znorm(warm_ipsi, blidx1, blidx2);
            % ersp_warm_contra(s, ch, :, :) = znorm(warm_contra, blidx1, blidx2);

            % Calculate ERSPs
            % ersp_cold_ipsi(s, ch, :, :)   = 10 * log10(bsxfun(@rdivide, squeeze(mean(cold_ipsi, 3)),   squeeze(mean(cold_ipsi(:,   blidx1 : blidx2, :), [2, 3]))));
            % ersp_cold_contra(s, ch, :, :) = 10 * log10(bsxfun(@rdivide, squeeze(mean(cold_contra, 3)), squeeze(mean(cold_contra(:, blidx1 : blidx2, :), [2, 3]))));
            % ersp_warm_ipsi(s, ch, :, :)   = 10 * log10(bsxfun(@rdivide, squeeze(mean(warm_ipsi, 3)),   squeeze(mean(warm_ipsi(:,   blidx1 : blidx2, :), [2, 3]))));
            % ersp_warm_contra(s, ch, :, :) = 10 * log10(bsxfun(@rdivide, squeeze(mean(warm_contra, 3)), squeeze(mean(warm_contra(:, blidx1 : blidx2, :), [2, 3]))));

            % Get alpha traces
            atra_cold_ipsi(s, ch, :)   = squeeze(mean(cold_ipsi(tf_freqs >= 8 & tf_freqs <= 12, :, :), [1, 3]));
            atra_cold_contra(s, ch, :) = squeeze(mean(cold_contra(tf_freqs >= 8 & tf_freqs <= 12, :, :), [1, 3]));
            atra_warm_ipsi(s, ch, :)   = squeeze(mean(warm_ipsi(tf_freqs >= 8 & tf_freqs <= 12, :, :), [1, 3]));
            atra_warm_contra(s, ch, :) = squeeze(mean(warm_contra(tf_freqs >= 8 & tf_freqs <= 12, :, :), [1, 3]));

            % Get alpha traces baselined
            tmp1 = 10 * log10(bsxfun(@rdivide, squeeze(mean(cold_ipsi, 3)),   squeeze(mean(cold_ipsi(:,   blidx1 : blidx2, :), [2, 3]))));
            tmp2 = 10 * log10(bsxfun(@rdivide, squeeze(mean(cold_contra, 3)), squeeze(mean(cold_contra(:, blidx1 : blidx2, :), [2, 3]))));
            tmp3 = 10 * log10(bsxfun(@rdivide, squeeze(mean(warm_ipsi, 3)),   squeeze(mean(warm_ipsi(:,   blidx1 : blidx2, :), [2, 3]))));
            tmp4 = 10 * log10(bsxfun(@rdivide, squeeze(mean(warm_contra, 3)), squeeze(mean(warm_contra(:, blidx1 : blidx2, :), [2, 3]))));

            atra_db_cold_ipsi(s, ch, :)   = squeeze(mean(tmp1(tf_freqs >= 8 & tf_freqs <= 12, :), 1));
            atra_db_cold_contra(s, ch, :) = squeeze(mean(tmp2(tf_freqs >= 8 & tf_freqs <= 12, :), 1));
            atra_db_warm_ipsi(s, ch, :)   = squeeze(mean(tmp3(tf_freqs >= 8 & tf_freqs <= 12, :), 1));
            atra_db_warm_contra(s, ch, :) = squeeze(mean(tmp4(tf_freqs >= 8 & tf_freqs <= 12, :), 1));

            tmp1 = znorm(cold_ipsi,   zblidx1, zblidx2);
            tmp2 = znorm(cold_contra, zblidx1, zblidx2);
            tmp3 = znorm(warm_ipsi,   zblidx1, zblidx2);
            tmp4 = znorm(warm_contra, zblidx1, zblidx2);

            atra_z_cold_ipsi(s, ch, :)   = squeeze(mean(tmp1(tf_freqs >= 8 & tf_freqs <= 12, :), 1));
            atra_z_cold_contra(s, ch, :) = squeeze(mean(tmp2(tf_freqs >= 8 & tf_freqs <= 12, :), 1));
            atra_z_warm_ipsi(s, ch, :)   = squeeze(mean(tmp3(tf_freqs >= 8 & tf_freqs <= 12, :), 1));
            atra_z_warm_contra(s, ch, :) = squeeze(mean(tmp4(tf_freqs >= 8 & tf_freqs <= 12, :), 1));

        end % End channel iteration
    end % End subject iteration

    % Save ERSPs
    % save([PATH_TFDECOMP 'ersp_cold_ipsi'], 'ersp_cold_ipsi');
    % save([PATH_TFDECOMP 'ersp_cold_contra'], 'ersp_cold_contra');
    % save([PATH_TFDECOMP 'ersp_warm_ipsi'], 'ersp_warm_ipsi');
    % save([PATH_TFDECOMP 'ersp_warm_contra'], 'ersp_warm_contra');

    save([PATH_TFDECOMP 'atra_cold_ipsi'],   'atra_cold_ipsi');
    save([PATH_TFDECOMP 'atra_cold_contra'], 'atra_cold_contra');
    save([PATH_TFDECOMP 'atra_warm_ipsi'],   'atra_warm_ipsi');
    save([PATH_TFDECOMP 'atra_warm_contra'], 'atra_warm_contra');

    save([PATH_TFDECOMP 'atra_db_cold_ipsi'],   'atra_db_cold_ipsi');
    save([PATH_TFDECOMP 'atra_db_cold_contra'], 'atra_db_cold_contra');
    save([PATH_TFDECOMP 'atra_db_warm_ipsi'],   'atra_db_warm_ipsi');
    save([PATH_TFDECOMP 'atra_db_warm_contra'], 'atra_db_warm_contra');

    save([PATH_TFDECOMP 'atra_z_cold_ipsi'],   'atra_z_cold_ipsi');
    save([PATH_TFDECOMP 'atra_z_cold_contra'], 'atra_z_cold_contra');
    save([PATH_TFDECOMP 'atra_z_warm_ipsi'],   'atra_z_warm_ipsi');
    save([PATH_TFDECOMP 'atra_z_warm_contra'], 'atra_z_warm_contra');

end % End part1

% ============================ Part 2: stats ============================================================================
if ismember('part2a', to_execute)

    % Load chanlocs
    EEG = pop_loadset('filename', [num2str(subject_list(1)) '_autocleaned.set'], 'filepath', PATH_AUTOCLEANED, 'loadmode', 'info');

    % Get tf params
    tf_times = dlmread([PATH_TFDECOMP 'tf_times.csv']); 
    tf_freqs = dlmread([PATH_TFDECOMP 'tf_freqs.csv']);

    % Load ERSPs
    load([PATH_TFDECOMP 'atra_z_cold_ipsi']);
    load([PATH_TFDECOMP 'atra_z_cold_contra']);
    load([PATH_TFDECOMP 'atra_z_warm_ipsi']);
    load([PATH_TFDECOMP 'atra_z_warm_contra']);


    cluster_channels = [18, 21, 23, 27, 24, 26]

    figure
    for chan = 1 : length(cluster_channels)

        ch = cluster_channels(chan);

        pd1 = squeeze(mean(atra_z_cold_ipsi(:, ch, :), 1));
        pd2 = squeeze(mean(atra_z_cold_contra(:, ch, :), 1));
        pd3 = squeeze(mean(atra_z_warm_ipsi(:, ch, :), 1));
        pd4 = squeeze(mean(atra_z_warm_contra(:, ch, :), 1));


        subplot(2, 3, chan)
        plot(tf_times, [pd1, pd2, pd3, pd4])
        title(EEG.chanlocs(ch).labels)

        if chan == 6
            legend({'coldips', 'coldcon', 'warmips', 'warmcon'})
        end
    end


end

% ============================ Part 2: stats ============================================================================
if ismember('part2', to_execute)

    % Init eeglab
    addpath(PATH_EEGLAB);
    eeglab;

    % Init fieldtrip
    addpath(PATH_FIELDTRIP);
    ft_defaults;

    % Get tf params
    tf_times = dlmread([PATH_TFDECOMP 'tf_times.csv']); 
    tf_freqs = dlmread([PATH_TFDECOMP 'tf_freqs.csv']);

    % Load ERSPs
    load([PATH_TFDECOMP 'ersp_cold_ipsi']);
    load([PATH_TFDECOMP 'ersp_cold_contra']);
    load([PATH_TFDECOMP 'ersp_warm_ipsi']);
    load([PATH_TFDECOMP 'ersp_warm_contra']);

    % Smooth data in time
    winlength = 50;
    if winlength
        ersp_cold_ipsi = movmean(ersp_cold_ipsi, winlength, 3);
        ersp_cold_contra = movmean(ersp_cold_contra, winlength, 3);
        ersp_warm_ipsi = movmean(ersp_warm_ipsi, winlength, 3);
        ersp_warm_contra = movmean(ersp_warm_contra, winlength, 3);
    end

    % Get channel labels and coordinates
    EEG = pop_loadset('filename', [num2str(subject_list(1)) '_autocleaned.set'], 'filepath', PATH_AUTOCLEANED, 'loadmode', 'info');
    chanlocs = EEG.chanlocs;
    chanlabs = {};
    coords = [];
    for c = 1 : numel(chanlocs)
        chanlabs{c} = chanlocs(c).labels;
        coords(c, :) = [chanlocs(c).X, chanlocs(c).Y, chanlocs(c).Z];
    end
 
    % A sensor struct
    sensors = struct();
    sensors.label = chanlabs;
    sensors.chanpos = coords;
    sensors.elecpos = coords;

    % Prepare neighbor struct
    cfg                 = [];
    cfg.elec            = sensors;
    cfg.feedback        = 'no';
    cfg.method          = 'triangulation'; 
    neighbours          = ft_prepare_neighbours(cfg);

    % A template for GA structs
    ga_template = [];
    ga_template.dimord = 'chan_freq_time';
    ga_template.label = chanlabs;
    ga_template.freq = tf_freqs;
    ga_template.time = tf_times;

    % Build condition ga
    cfg=[];
    cfg.keepindividual = 'yes';
    D = {};
    for s = 1 : length(subject_list)
        ga_template.powspctrm = squeeze(ersp_cold_ipsi(s, :, :, :));
        D{s} = ga_template;
    end 
    GA_cold_ipsi = ft_freqgrandaverage(cfg, D{1, :});
    D = {};
    for s = 1 : length(subject_list)
        ga_template.powspctrm = squeeze(ersp_cold_contra(s, :, :, :));
        D{s} = ga_template;
    end 
    GA_cold_contra = ft_freqgrandaverage(cfg, D{1, :});
    D = {};
    for s = 1 : length(subject_list)
        ga_template.powspctrm = squeeze(ersp_warm_ipsi(s, :, :, :));
        D{s} = ga_template;
    end 
    GA_warm_ipsi = ft_freqgrandaverage(cfg, D{1, :});
    D = {};
    for s = 1 : length(subject_list)
        ga_template.powspctrm = squeeze(ersp_warm_contra(s, :, :, :));
        D{s} = ga_template;
    end 
    GA_warm_contra = ft_freqgrandaverage(cfg, D{1, :});

    % Build difference ga
    cfg=[];
    cfg.keepindividual = 'yes';
    D = {};
    for s = 1 : length(subject_list)
        ga_template.powspctrm = squeeze(ersp_cold_ipsi(s, :, :, :)) - squeeze(ersp_cold_contra(s, :, :, :));
        D{s} = ga_template;
    end 
    GA_diff_in_cold = ft_freqgrandaverage(cfg, D{1, :});
    D = {};
    for s = 1 : length(subject_list)
        ga_template.powspctrm = squeeze(ersp_warm_ipsi(s, :, :, :)) - squeeze(ersp_warm_contra(s, :, :, :));
        D{s} = ga_template;
    end 
    GA_diff_in_warm = ft_freqgrandaverage(cfg, D{1, :});
    D = {};
    for s = 1 : length(subject_list)
        ga_template.powspctrm = squeeze(ersp_warm_ipsi(s, :, :, :)) - squeeze(ersp_cold_ipsi(s, :, :, :));
        D{s} = ga_template;
    end 
    GA_diff_in_ipsi = ft_freqgrandaverage(cfg, D{1, :});
    D = {};
    for s = 1 : length(subject_list)
        ga_template.powspctrm = squeeze(ersp_warm_contra(s, :, :, :)) - squeeze(ersp_cold_contra(s, :, :, :));
        D{s} = ga_template;
    end 
    GA_diff_in_contra = ft_freqgrandaverage(cfg, D{1, :});

    % Testparams
    testalpha  = 0.025;
    voxelalpha  = 0.01;
    nperm = 1000;

    % Set config. Same for all tests
    cfg = [];
    cfg.tail             = 0; % Two sided test!
    cfg.statistic        = 'depsamplesT';
    cfg.alpha            = testalpha;
    cfg.neighbours       = neighbours;
    cfg.minnbchan        = 2;
    cfg.method           = 'montecarlo';
    cfg.correctm         = 'cluster';
    cfg.clustertail      = 0;
    cfg.clusteralpha     = voxelalpha;
    cfg.clusterstatistic = 'maxsum';
    cfg.numrandomization = nperm;
    cfg.computecritval   = 'yes'; 
    cfg.ivar             = 1;
    cfg.uvar             = 2;
    cfg.design           = [ones(1, numel(subject_list)), 2 * ones(1, numel(subject_list)); 1 : numel(subject_list), 1 : numel(subject_list)];

    % The tests
    [stat_interaction] = ft_freqstatistics(cfg, GA_cold_contra, GA_cold_ipsi);  


    % Calculate effect sizes
    adjpetasq_cold = [];
    for ch = 1 : 32

        petasq = (squeeze(stat_interaction.stat(ch, :, :)) .^ 2) ./ ((squeeze(stat_interaction.stat(ch, :, :)) .^ 2) + (numel(subject_list) - 1));
        adj_petasq = petasq - (1 - petasq) .* (1 / (numel(subject_list) - 1));
        adjpetasq_cold(ch, :, :) = adj_petasq;


    end




    %pd = GA_warm_contra;
    %topoplot(pd, chanlocs, 'plotrad', 0.7, 'intrad', 0.7, 'intsquare', 'on', 'conv', 'off', 'electrodes', 'on');






    % Repeat testalpha
    testalpha  = 0.025;

    % Set colors
    cmap = 'jet';
    clinecol = 'k';

    % Identify significant clusters
    clusts = struct();
    cnt = 0;
    stat_names = {'stat_interaction'};
    for s = 1 : numel(stat_names)
        stat = eval(stat_names{s});
        if ~isempty(stat.negclusters)
            neg_idx = find([stat.negclusters(1, :).prob] < testalpha);
            for c = 1 : numel(neg_idx)
                cnt = cnt + 1;
                clusts(cnt).testlabel = stat_names{s};
                clusts(cnt).clustnum = cnt;
                clusts(cnt).time = stat.time;
                clusts(cnt).freq = stat.freq;
                clusts(cnt).polarity = -1;
                clusts(cnt).prob = stat.negclusters(1, neg_idx(c)).prob;
                clusts(cnt).idx = stat.negclusterslabelmat == neg_idx(c);
                clusts(cnt).stats = clusts(cnt).idx .* stat.stat * -1;
                clusts(cnt).chans_sig = find(logical(mean(clusts(cnt).idx, [2,3])));
            end
        end
        if ~isempty(stat.posclusters)
            pos_idx = find([stat.posclusters(1, :).prob] < testalpha);
            for c = 1 : numel(pos_idx)
                cnt = cnt + 1;
                clusts(cnt).testlabel = stat_names{s};
                clusts(cnt).clustnum = cnt;
                clusts(cnt).time = stat.time;
                clusts(cnt).freq = stat.freq;
                clusts(cnt).polarity = 1;
                clusts(cnt).prob = stat.posclusters(1, pos_idx(c)).prob;
                clusts(cnt).idx = stat.posclusterslabelmat == pos_idx(c);
                clusts(cnt).stats = clusts(cnt).idx .* stat.stat;
                clusts(cnt).chans_sig = find(logical(mean(clusts(cnt).idx, [2,3])));
            end
        end
    end

    % Plot identified cluster
    for cnt = 1 : numel(clusts)
        figure('Visible', 'off'); clf;
        subplot(2, 2, 1)
        pd = squeeze(sum(clusts(cnt).stats, 1));
        contourf(clusts(cnt).time, clusts(cnt).freq, pd, 40, 'linecolor','none')
        hold on
        contour(clusts(cnt).time, clusts(cnt).freq, logical(squeeze(mean(clusts(cnt).idx, 1))), 1, 'linecolor', clinecol, 'LineWidth', 2)
        colormap(cmap)
        set(gca, 'xlim', [clusts(cnt).time(1), clusts(cnt).time(end)], 'clim', [-max(abs(pd(:))), max(abs(pd(:)))], 'YScale', 'log', 'YTick', [4, 8, 12, 20, 30])
        colorbar;
        title(['sum t across chans, plrt: ' num2str(clusts(cnt).polarity)], 'FontSize', 10)
        subplot(2, 2, 2)
        pd = squeeze(mean(clusts(cnt).idx, 1));
        contourf(clusts(cnt).time, clusts(cnt).freq, pd, 40, 'linecolor','none')
        hold on
        contour(clusts(cnt).time, clusts(cnt).freq, logical(squeeze(mean(clusts(cnt).idx, 1))), 1, 'linecolor', clinecol, 'LineWidth', 2)
        colormap(cmap)
        set(gca, 'xlim', [clusts(cnt).time(1), clusts(cnt).time(end)], 'clim', [-1, 1], 'YScale', 'log', 'YTick', [4, 8, 12, 20, 30])
        colorbar;
        title(['proportion chans significant'], 'FontSize', 10)
        subplot(2, 2, 3)
        pd = squeeze(sum(clusts(cnt).stats, [2, 3]));
        topoplot(pd, chanlocs, 'plotrad', 0.7, 'intrad', 0.7, 'intsquare', 'on', 'conv', 'off', 'electrodes', 'on');
        colormap(cmap)
        set(gca, 'clim', [-max(abs(pd(:))), max(abs(pd(:)))])
        colorbar;
        title(['sum t per electrode'], 'FontSize', 10)
        subplot(2, 2, 4)
        pd = squeeze(mean(clusts(cnt).idx, [2, 3]));
        topoplot(pd, chanlocs, 'plotrad', 0.7, 'intrad', 0.7, 'intsquare', 'on', 'conv', 'off', 'electrodes', 'on');
        colormap(cmap)
        set(gca, 'clim', [-1, 1])
        colorbar;
        title(['proportion tf-points significant'], 'FontSize', 10)
        saveas(gcf, [PATH_PLOT 'clustnum_' num2str(clusts(cnt).clustnum) '_' clusts(cnt).testlabel '.png']); 
    end






end % End part2


% Function for z norming powcubes
function[zersp] = znorm(indata, blidx1, blidx2)
    zdata = zeros(size(indata));
    tmp = indata(:, blidx1 : blidx2, :);
    tmp = reshape(tmp, [size(tmp, 1), size(tmp, 2) * size(tmp, 3)]);
    mpow = repmat(mean(tmp, 2), [1, size(indata, 2)]);
    stdpow = repmat(std(tmp, [], 2), [1, size(indata, 2)]);
    for t = 1 : size(indata, 3)
        zdata(:, :, t) = (squeeze(indata(:, :, t)) - mpow) ./ stdpow;
    end
    zersp = squeeze(mean(zdata, 3));
end