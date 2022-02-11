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

    % Init eeglab
    addpath(PATH_EEGLAB);
    eeglab;

    % Load chanlocs
    EEG = pop_loadset('filename', [num2str(subject_list(1)) '_autocleaned.set'], 'filepath', PATH_AUTOCLEANED, 'loadmode', 'info');

    % Get tf params
    tf_times = dlmread([PATH_TFDECOMP 'tf_times.csv']); 
    tf_freqs = dlmread([PATH_TFDECOMP 'tf_freqs.csv']);

    % Load ERSPs
    load([PATH_TFDECOMP 'atra_cold_ipsi']);
    load([PATH_TFDECOMP 'atra_cold_contra']);
    load([PATH_TFDECOMP 'atra_warm_ipsi']);
    load([PATH_TFDECOMP 'atra_warm_contra']);
    load([PATH_TFDECOMP 'atra_z_cold_ipsi']);
    load([PATH_TFDECOMP 'atra_z_cold_contra']);
    load([PATH_TFDECOMP 'atra_z_warm_ipsi']);
    load([PATH_TFDECOMP 'atra_z_warm_contra']);



    % Define channels in cluster
    cluster_channels = [18, 21, 23, 27, 24, 26];

    % Build matrices for ANOVA
    anova_matrix = [];
    for s = 1 : size(atra_cold_ipsi, 1)
        tidx = tf_times >= -500 & tf_times < -200;
        anova_matrix(s, 1, 1) = mean2(atra_cold_ipsi(s, cluster_channels, tidx));
        anova_matrix(s, 1, 2) = mean2(atra_cold_contra(s, cluster_channels, tidx));
        anova_matrix(s, 2, 1) = mean2(atra_warm_ipsi(s, cluster_channels, tidx));
        anova_matrix(s, 2, 2) = mean2(atra_warm_contra(s, cluster_channels, tidx));
    end
    t = simple_mixed_anova(anova_matrix, [], {'stress', 'hemisqhere'});

    % Build matrices for ANOVA
    anova_matrix = [];
    for s = 1 : size(atra_cold_ipsi, 1)
        tidx = tf_times >= -500 & tf_times < -200;
        anova_matrix(s, 1, 1) = mean2(atra_z_cold_ipsi(s, cluster_channels, tidx));
        anova_matrix(s, 1, 2) = mean2(atra_z_cold_contra(s, cluster_channels, tidx));
        anova_matrix(s, 2, 1) = mean2(atra_z_warm_ipsi(s, cluster_channels, tidx));
        anova_matrix(s, 2, 2) = mean2(atra_z_warm_contra(s, cluster_channels, tidx));
    end
    t_z = simple_mixed_anova(anova_matrix, [], {'stress', 'hemisqhere'});

    % Build matrices for ANOVA
    anova_matrix = [];
    for s = 1 : size(atra_cold_ipsi, 1)
        tidx = tf_times >= 900 & tf_times < 1500;
        anova_matrix(s, 1, 1) = mean2(atra_cold_ipsi(s, cluster_channels, tidx));
        anova_matrix(s, 1, 2) = mean2(atra_cold_contra(s, cluster_channels, tidx));
        anova_matrix(s, 2, 1) = mean2(atra_warm_ipsi(s, cluster_channels, tidx));
        anova_matrix(s, 2, 2) = mean2(atra_warm_contra(s, cluster_channels, tidx));
    end
    t_post = simple_mixed_anova(anova_matrix, [], {'stress', 'hemisqhere'});

    % Build matrices for ANOVA
    anova_matrix = [];
    for s = 1 : size(atra_cold_ipsi, 1)
        tidx = tf_times >= 900 & tf_times < 1500;
        anova_matrix(s, 1, 1) = mean2(atra_z_cold_ipsi(s, cluster_channels, tidx));
        anova_matrix(s, 1, 2) = mean2(atra_z_cold_contra(s, cluster_channels, tidx));
        anova_matrix(s, 2, 1) = mean2(atra_z_warm_ipsi(s, cluster_channels, tidx));
        anova_matrix(s, 2, 2) = mean2(atra_z_warm_contra(s, cluster_channels, tidx));
    end
    t_post_z = simple_mixed_anova(anova_matrix, [], {'stress', 'hemisqhere'});

    % Build matrices for ANOVA
    anova_matrix = [];
    for s = 1 : size(atra_cold_ipsi, 1)
        tidx = tf_times >= -500 & tf_times < 0;
        anova_matrix(s, 1) = (mean2(atra_cold_ipsi(s, cluster_channels, tidx)) + mean2(atra_cold_contra(s, cluster_channels, tidx))) / 2;
        anova_matrix(s, 2) = (mean2(atra_warm_ipsi(s, cluster_channels, tidx)) + mean2(atra_warm_contra(s, cluster_channels, tidx))) / 2;
    end
    t_bl = simple_mixed_anova(anova_matrix, [], {'stress'});


    % Channels averaged
    figure

    % Averaged cluster channels, rawpow
    subplot(2, 2, 1)
    pd1 = squeeze(mean(atra_cold_ipsi(:, [18, 21, 23, 27, 24, 26], :), [1, 2]));
    pd2 = squeeze(mean(atra_cold_contra(:, [18, 21, 23, 27, 24, 26], :), [1, 2]));
    pd3 = squeeze(mean(atra_warm_ipsi(:, [18, 21, 23, 27, 24, 26], :), [1, 2]));
    pd4 = squeeze(mean(atra_warm_contra(:, [18, 21, 23, 27, 24, 26], :), [1, 2]));
    hold on
    plot(tf_times, pd1, '-.', 'Color', '#4DBEEE', 'LineWidth', 2.5)
    plot(tf_times, pd2, '-', 'Color', '#4DBEEE', 'LineWidth', 2.5)
    plot(tf_times, pd3, '-.', 'Color', 'k', 'LineWidth', 2.5)
    plot(tf_times, pd4, '-', 'Color', 'k', 'LineWidth', 2.5)
    hold off
    grid on
    xlabel('ms')
    ylabel('\muV^2') 
    title('alpha raw power')

    % Averaged cluster channels, z-normalized
    subplot(2, 2, 2)
    pd1 = squeeze(mean(atra_z_cold_ipsi(:, [18, 21, 23, 27, 24, 26], :), [1, 2]));
    pd2 = squeeze(mean(atra_z_cold_contra(:, [18, 21, 23, 27, 24, 26], :), [1, 2]));
    pd3 = squeeze(mean(atra_z_warm_ipsi(:, [18, 21, 23, 27, 24, 26], :), [1, 2]));
    pd4 = squeeze(mean(atra_z_warm_contra(:, [18, 21, 23, 27, 24, 26], :), [1, 2]));
    hold on
    plot(tf_times, pd1, '-.', 'Color', '#4DBEEE', 'LineWidth', 2.5)
    plot(tf_times, pd2, '-', 'Color', '#4DBEEE', 'LineWidth', 2.5)
    plot(tf_times, pd3, '-.', 'Color', 'k', 'LineWidth', 2.5)
    plot(tf_times, pd4, '-', 'Color', 'k', 'LineWidth', 2.5)
    hold off
    grid on
    xlabel('ms')
    ylabel('z') 
    title('alpha z-normalized')
    legend({'stress ipsi', 'stress contra', 'control ipsi', 'control contra'})

    % Difference wave raw
    subplot(2, 2, 3)
    pd1 = squeeze(mean(atra_cold_ipsi(:, [18, 21, 23, 27, 24, 26], :), [1, 2])) - squeeze(mean(atra_cold_contra(:, [18, 21, 23, 27, 24, 26], :), [1, 2]));
    pd2 = squeeze(mean(atra_warm_ipsi(:, [18, 21, 23, 27, 24, 26], :), [1, 2])) - squeeze(mean(atra_warm_contra(:, [18, 21, 23, 27, 24, 26], :), [1, 2]));
    hold on
    plot(tf_times, pd1, '-', 'Color', '#4DBEEE', 'LineWidth', 2.5)
    plot(tf_times, pd2, '-', 'Color', 'k', 'LineWidth', 2.5)
    hold off
    grid on
    xlabel('ms')
    ylabel('\muV^2') 
    title('alpha raw power')

    % Difference wave z-normalized
    subplot(2, 2, 4)
    pd1 = squeeze(mean(atra_z_cold_ipsi(:, [18, 21, 23, 27, 24, 26], :), [1, 2])) - squeeze(mean(atra_z_cold_contra(:, [18, 21, 23, 27, 24, 26], :), [1, 2]));
    pd2 = squeeze(mean(atra_z_warm_ipsi(:, [18, 21, 23, 27, 24, 26], :), [1, 2])) - squeeze(mean(atra_z_warm_contra(:, [18, 21, 23, 27, 24, 26], :), [1, 2]));
    hold on
    plot(tf_times, pd1, '-', 'Color', '#4DBEEE', 'LineWidth', 2.5)
    plot(tf_times, pd2, '-', 'Color', 'k', 'LineWidth', 2.5)
    hold off
    grid on
    xlabel('ms')
    ylabel('z') 
    title('alpha z-normalized')
    legend({'stress', 'control'})
    

end