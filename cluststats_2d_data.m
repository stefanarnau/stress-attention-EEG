function [sig_flag, mean_data1, mean_data2, clust_outlines, apes, sum_t, p_values, mean_apes, time_limits, freq_limits, cluster_idx] = cluststats_2d_data(data1, data2, tf_times, tf_freqs, varargin)

    %
    % Function performs a cluster based permutation test on 2d data.
    %
    % Inputs:
    % ------------------
    %
    % data1 & data2: 2d data as dimord (unit, dimension1, dimension2).
    % tf_times: Vector of latencies
    % tf_freqs: Vector of frequencies
    %
    % Optional inputs:
    % ------------------
    %
    % tail:                      Can be 0 (default), 1 or -1. Corresponding to two-sided, left and right testing.
    % pval_voxel:                Thresholding p value (default 0.01).
    % pval_cluster:              Significance level for inference statistic (default 0.05).
    % n_perms:                   Number of permutations for creating the null-hypothesis distribution (default 1000).
    %
    % Returns:
    % ------------------
    %
    % sig_flag:                  The flag of significance. 1 if there are significant differences, 0 if not.          
    % mean_data1 & mean_data2:   The datasets averaged across unit dimension.
    % clust_outlines:            Outlines of identified clusters.
    % apes:                      The effect sizes (Mordkoff, 2019).
    % 
    % Also for each cluster: summed t -values, p-values, mean effect size, time limits, frequency limits, cluster indices
    %
    % Stefan Arnau, June 2020
    % Email: arnau@ifado.de
    %

    % Set defaults
    default_tail         = 0;
    default_pval_voxel   = 0.01;
    default_pval_cluster = 0.05;
    default_n_perms      = 1000;
    
    % Init input parser
    p = inputParser;

    % Parse inputs and set defaults
    p.FunctionName  = mfilename;
    p.CaseSensitive = false;
    p.addRequired('data1', @isnumeric);
    p.addRequired('data2', @isnumeric);
    p.addRequired('tf_times', @isnumeric);
    p.addRequired('tf_freqs', @isnumeric);
    p.addParamValue('tail', default_tail, @isnumeric);
    p.addParamValue('pval_voxel', default_pval_voxel, @isnumeric);
    p.addParamValue('pval_cluster', default_pval_cluster, @isnumeric);
    p.addParamValue('n_perms', default_n_perms, @isnumeric);
    parse(p, data1, data2, tf_times, tf_freqs, varargin{:});

    % Adjust cluster-p-value for 2-sided testing
    if p.Results.tail
        cpv = p.Results.pval_cluster;
    else
        cpv = p.Results.pval_cluster / 2;
    end

    % Catch some of the possible errors
	if size(p.Results.data2) ~= size(p.Results.data1)
		error('Datasets must be of equal shape... :)');
		return;
	end
    
    % Init matrices
    permuted_t = zeros(p.Results.n_perms, size(p.Results.data1, 2), size(p.Results.data1, 3));
    max_clust = zeros(p.Results.n_perms, 2);
    desmat = [zeros(size(p.Results.data1, 1), 1), ones(size(p.Results.data1, 1), 1)];

    % Create test-statistic distribution on permuted data
    for perm = 1 : p.Results.n_perms
        fprintf('\nPermutation %i/%i\n', perm, p.Results.n_perms);
        toflip = find(round(rand(size(p.Results.data1, 1), 1)));
        d1_perm = p.Results.data1;
        d1_perm(toflip, :, :) = p.Results.data2(toflip, :, :);
        d2_perm = p.Results.data2;
        d2_perm(toflip, :, :) = p.Results.data1(toflip, :, :);
        tnum = squeeze(mean(d1_perm - d2_perm, 1));
        tdenum = squeeze(std(d1_perm - d2_perm, 0, 1)) / sqrt(size(p.Results.data1, 1));
        fake_t = tnum ./ tdenum;
        permuted_t(perm, :, :) = fake_t;
        fake_t(abs(fake_t) < tinv(1 - p.Results.pval_voxel, size(p.Results.data1, 1) - 1)) = 0;
        clusts = bwconncomp(fake_t);
        sum_t = [];
        for clu = 1 : numel(clusts.PixelIdxList)
            cidx = clusts.PixelIdxList{clu};
            sum_t(end + 1) = sum(fake_t(cidx));
        end
        max_clust(perm, 1) = min([0, sum_t]);
        max_clust(perm, 2) = max([0, sum_t]);      
    end

    % T-test observed data
    tnum = squeeze(mean(p.Results.data1 - p.Results.data2, 1));
    tdenum = squeeze(std(p.Results.data1 - p.Results.data2, 0, 1)) / sqrt(size(p.Results.data1, 1));
    tmat = tnum ./ tdenum;
    tvals = tmat;
    tmat(abs(tmat) < tinv(1 - p.Results.pval_voxel, size(p.Results.data1, 1) - 1)) = 0;

    % Find clusters in observed data
    clusts = bwconncomp(tmat);

    % Get cluster test statistic
    sum_t = [];
    cluster_idx = {};
    for clu = 1 : numel(clusts.PixelIdxList)
        cidx = clusts.PixelIdxList{clu};
        cluster_idx{clu} = cidx;
        sum_t(end + 1) = sum(tmat(cidx));
    end

    % Remove non-significant clusters
    clust_thresh_lower = prctile(max_clust(:, 1), cpv * 100);
    clust_thresh_upper = prctile(max_clust(:, 2), 100 - cpv * 100);
    clust2remove = find(sum_t > clust_thresh_lower & sum_t < clust_thresh_upper);
    for clu = 1 : length(clust2remove)
        tmat(clusts.PixelIdxList{clust2remove(clu)}) = 0;
    end
    cluster_idx(clust2remove) = [];
    sum_t(clust2remove) = [];

    % Get time and frequency boundaries of cluster
    time_limits = {};
    freq_limits = {};
    for clu = 1 : numel(cluster_idx)
        binmat = zeros(size(tmat));
        binmat(cluster_idx{clu}) = 1;
        time_limits{clu} = [min(p.Results.tf_times(logical(sum(binmat, 1)))), max(p.Results.tf_times(logical(sum(binmat, 1))))];
        freq_limits{clu} = [min(p.Results.tf_freqs(logical(sum(binmat, 2)))), max(p.Results.tf_freqs(logical(sum(binmat, 2))))];        
    end

    % Get cluster contour
    clust_outlines = logical(tmat);

    % Set the flag of significance
    sig_flag = logical(sum(clust_outlines(:)));

    % Calculate effect sizes
    x = tvals.^2 ./ (tvals.^2 + (size(p.Results.data1, 1) - 1));
    apes = x - (1 - x) .* (1 / (size(p.Results.data1, 1) - 1));
    mean_apes = [];
    for clu = 1 : numel(cluster_idx)
        mean_apes(end + 1) = mean(apes(cluster_idx{clu}));
    end

    % Estimate p-value for cluster
    p_values = [];
    min_dist = sort(max_clust(:, 1));
    max_dist = sort(max_clust(:, 2));
    for clu = 1 : numel(cluster_idx)
        pv = NaN;
        if p.Results.tail == 0
            if sum_t(clu) < min_dist(1)
                pv = 0;
            elseif sum_t(clu) > max_dist(end)
                pv = 0;
            else
                [~, min_idx] = min(abs(min_dist - sum_t(clu)));
                p_min = min_idx / p.Results.n_perms;
                [~, max_idx] = min(abs(max_dist - sum_t(clu)));
                p_max = 1 - (max_idx / p.Results.n_perms);
                pv = min(p_min, p_max);
            end
        elseif p.Results.tail == -1
            if sum_t(clu) < min_dist(1)
                pv = 0;
            elseif sum_t(clu) > max_dist(end)
                pv = 1;
            else
                [~, min_idx] = min(abs(min_dist - sum_t(clu)));
                pv = min_idx / p.Results.n_perms;
            end
        elseif p.Results.tail == 1
            if sum_t(clu) < min_dist(1)
                pv = 1;
            elseif sum_t(clu) > max_dist(end)
                pv = 0;
            else
                [~, max_idx] = min(abs(max_dist - sum_t(clu)));
                pv = 1 - (max_idx / p.Results.n_perms);
            end
        end
        p_values(end + 1) = pv;
    end

    % Calculate averages
    mean_data1 = squeeze(mean(p.Results.data1, 1));
    mean_data2 = squeeze(mean(p.Results.data2, 1));

end