clear all;

% Set paths
PATH_EEGLAB           = 'add_a_nice_path_right_here';
PATH_CHANLOCFILE      = 'add_a_nice_path_right_here';
PATH_RAW_DATA         = 'add_a_nice_path_right_here'; 
PATH_ICSET            = 'add_a_nice_path_right_here';
PATH_AUTOCLEANED      = 'add_a_nice_path_right_here';
PATH_TFDECOMP         = 'add_a_nice_path_right_here'; 

% ======================= SUBJECTS =========================================================================================================

% Define subjects (N=32)
subject_list = [1, 2, 4, 5, 6, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 23, 24, 25, 26, 27, 28];

% Define subject test list
subject_test_list = 5;

% ======================= OPTIONS =========================================================================================================

% Switch parts of the script on/off
to_execute = {'part1', 'part2', 'part3'};

% ======================= PART1: PREPROCESSING ===================================================================================================
if ismember('part1', to_execute)

    % Init EEGlab
    addpath(PATH_EEGLAB);
    [ALLEEG, EEG, CURRENTSET, ALLCOM] = eeglab;

    % Find chanlocfile
    channel_location_file = [PATH_CHANLOCFILE, 'standard-10-5-cap385.elp'];
    
    % Read session info
    session_info = dlmread([PATH_RAW_DATA, 'trialinfo/condition_session_info.csv']);

	% Iterating subject list
	for s = 1 : length(subject_list)

        % Iterate sessions
        for ses = 1 : 2

            % Load stuff
            id = subject_list(s);
            trinfo_l = xlsread([PATH_RAW_DATA 'trialinfo/' num2str(id) '_' num2str(ses) '_l.xlsx']);
            trinfo_r = xlsread([PATH_RAW_DATA 'trialinfo/' num2str(id) '_' num2str(ses) '_r.xlsx']);
            EEG =  pop_loadbv(PATH_RAW_DATA, ['Exp1027_Posner_' num2str(id) '_' num2str(ses) '.vhdr'], [], []);
            EEG = pop_select(EEG, 'channel', [3 : 34]);
            EEG = pop_chanedit(EEG, 'lookup', channel_location_file);
            EEG.chanlocs_original = EEG.chanlocs;

            % New event struct
            new_events = struct('latency', {},...
                                'type', {},...
                                'code', {},...
                                'id', {},... 
                                'cuedir', {},...
                                'trialvalid', {},...
                                'targdir', {},...  
                                'accuracy', {},... 
                                'cuevalid', {},...      
                                'session', {},...
                                'watertemp', {},...             
                                'urevent', {},...
                                'duration', {}...
                                );
            
            % Loop events
            ecount = [0, 0, 0];
            for e = 1 : length(EEG.event)
                if ismember(EEG.event(e).type,  {'S 14', 'S 24'})
                    ecount(1) = ecount(1) + 1;
                    new_events(ecount(1)).latency = EEG.event(e).latency;
                    new_events(ecount(1)).type = 'cue';
                    new_events(ecount(1)).code = 'cue';
                    new_events(ecount(1)).id = id;
                    if str2num(EEG.event(e).type(3)) == 1
                        ecount(2) = ecount(2) + 1;
                        new_events(ecount(1)).cuedir =  1; % 1 = left
                        new_events(ecount(1)).trialvalid = trinfo_l(ecount(2), 2);
                    else
                        ecount(3) = ecount(3) + 1;
                        new_events(ecount(1)).cuedir = 2; % 2 = right
                        new_events(ecount(1)).trialvalid = trinfo_r(ecount(3), 2);
                    end
                    f = e;
                    while ~ismember(EEG.event(f).type,  {'S 15', 'S 25'})
                        f = f + 1;
                    end
                    if str2num(EEG.event(f).type(3)) == 1
                        new_events(ecount(1)).targdir = 1; % 1 = left
                    else
                        new_events(ecount(1)).targdir = 2; % 2 = right
                    end
                    f = e;
                    while ~ismember(EEG.event(f).type,  {'S 55', 'S 99'})
                        f = f + 1;
                    end
                    if str2num(EEG.event(f).type(3)) == 5
                        new_events(ecount(1)).accuracy = 1;
                    else
                        new_events(ecount(1)).accuracy = 0;
                    end
                    if strcmpi(new_events(ecount(1)).cuedir, new_events(ecount(1)).targdir)
                        new_events(ecount(1)).cuevalid = 1;
                    else
                        new_events(ecount(1)).cuevalid = 0;
                    end
                    new_events(ecount(1)).session = ses;
                    if ses == session_info(session_info(:, 1) == id, 2)
                        new_events(ecount(1)).watertemp = 1; % 1 = cold
                    else
                        new_events(ecount(1)).watertemp = 2; % 2 = warm
                    end
                    new_events(ecount(1)).urevent = ecount(1);
                    new_events(ecount(1)).duration = 1;
                end
            end

            % Save number of trials identified (must be 288)
            EEG.n_trials_original = ecount(1);

            % Add existing boundaries
            for e = 1 : length(EEG.event)
                if strcmpi(EEG.event(e).type, 'boundary')
                    ecount(1) = ecount(1) + 1;
                    new_events(ecount(1)).latency = EEG.event(e).latency;
                    new_events(ecount(1)).type = 'boundary';
                    new_events(ecount(1)).code = 'boundary';
                    new_events(ecount(1)).id = NaN;
                    new_events(ecount(1)).cuedir = NaN;
                    new_events(ecount(1)).trialvalid = NaN;
                    new_events(ecount(1)).targdir = NaN;
                    new_events(ecount(1)).accuracy = NaN;
                    new_events(ecount(1)).cuevalid = NaN;
                    new_events(ecount(1)).session = NaN;
                    new_events(ecount(1)).watertemp = NaN;
                    new_events(ecount(1)).urevent = ecount(1);
                    new_events(ecount(1)).duration = 1;
                end
            end

            % Replace events by new events
            EEG.event = new_events;
            EEG = eeg_checkset(EEG, 'eventconsistency');
    
            % Bandpass filter data (ERPlab toolbox function) 
            EEG = pop_basicfilter(EEG, [1 : EEG.nbchan], 'Cutoff', [1, 30], 'Design', 'butter', 'Filter', 'bandpass', 'Order', 4, 'RemoveDC', 'on', 'Boundary', 'boundary');
    
            % Bad channel detection
            [EEG, i1] = pop_rejchan(EEG, 'elec', [1 : EEG.nbchan], 'threshold', 10, 'norm', 'on', 'measure', 'kurt');
            [EEG, i2] = pop_rejchan(EEG, 'elec', [1 : EEG.nbchan], 'threshold', 5, 'norm', 'on', 'measure', 'prob');
            EEG.chans_rejected = horzcat(i1, i2);
            EEG.chans_rejected_n = length(horzcat(i1, i2));

            % Reref to common average reference
            EEG = pop_reref(EEG, []);

            % Resample data
            EEG = pop_resample(EEG, 200);

            % Epoch data
            EEG = pop_epoch(EEG, {'cue'}, [-1.5, 4], 'newname', [num2str(id) '_seg'], 'epochinfo', 'yes');
            
            % Autoreject data before ICA
            EEG.segs_original_n = size(EEG.data, 3);
            [EEG, rejsegs] = pop_autorej(EEG, 'nogui', 'on', 'threshold', 1000, 'startprob', 5, 'maxrej', 5, 'eegplot', 'off');
            EEG.segs_rejected_before_ica = length(rejsegs);

            % Run ICA
            EEG = pop_runica(EEG, 'extended', 1, 'interupt', 'on');

            % Save IC set
            EEG = pop_saveset(EEG, 'filename', [num2str(id) '_' num2str(ses) '_icset.set'], 'filepath', PATH_ICSET, 'check', 'on', 'savemode', 'twofiles');

        end
    end
end

% ======================= PART2: IC-BASED DATA CLEANING ====================================================================================
if ismember('part2', to_execute)

    % Init EEGlab
    addpath(PATH_EEGLAB);
    [ALLEEG, EEG, CURRENTSET, ALLCOM] = eeglab;

    preprostats = [];

    % Iterating subject list
    cnt = 0;
    for s = 1 : length(subject_list)

        % Iterate session
        for ses = 1 : 2

            % Get subject and session info
            cnt = cnt + 1;
            id = subject_list(s);
            preprostats(cnt, 1) = id;
            preprostats(cnt, 2) = ses;

            % Load data
            EEG = pop_loadset('filename', [num2str(id) '_' num2str(ses) '_icset.set'], 'filepath', PATH_ICSET, 'loadmode', 'all');
            preprostats(cnt, 3) = EEG.chans_rejected_n;

            % Run IClabel
            EEG = iclabel(EEG);
            EEG.ICout_IClabel = find(EEG.etc.ic_classification.ICLabel.classifications(:, 1) < 0.5);
            EEG = pop_subcomp(EEG, EEG.ICout_IClabel, 0);
            preprostats(cnt, 4) = length(EEG.ICout_IClabel);

            % Autoreject data again
            [EEG, rejsegs] = pop_autorej(EEG, 'nogui', 'on', 'threshold', 1000, 'startprob', 5, 'maxrej', 5);
            EEG.segs_rejected_after_ica = length(rejsegs);
            EEG.segs_rejected_overall_percentage = ((EEG.segs_rejected_before_ica + EEG.segs_rejected_after_ica) / EEG.segs_original_n) * 100;
            preprostats(cnt, 5) = EEG.segs_rejected_overall_percentage;

            % Interpolate missing channels
            EEG = pop_interp(EEG, EEG.chanlocs_original, 'spherical');

            % Remove bad trials
            EEG = pop_select(EEG, 'notrial', find(cell2mat({EEG.event.trialvalid}) == 0));

            % Save number of remaining trials
            preprostats(cnt, 6) = size(EEG.data, 3);

            % Backup to merge
            if ses == 1
                X = EEG;
            end
         
        end

        % Merge and save preprocessed EEG data
        EEG = pop_mergeset(X, EEG, 0);
        EEG = pop_saveset(EEG, 'filename', [num2str(id) '_autocleaned.set'], 'filepath', PATH_AUTOCLEANED, 'check', 'on', 'savemode', 'twofiles');

    end % End subject loop

    % Save preprocesssing statistics
    dlmwrite([PATH_AUTOCLEANED 'preprostats.csv'], preprostats);

end

% ======================= PART3: TIME FREQUENCY DECOMPOSITION =====================================================================
if ismember('part3', to_execute)

    % Init EEGlab
    addpath(PATH_EEGLAB);
    addpath(PATH_CUSTOM);
    eeglab;

    % Get some complex Morlet wavelets
	n_freqs = 50;
	frqrange = [2, 20];
	[cmw, tf_freqs] = get_cmwset(frqrange, n_freqs, 'lin', 200, 'time', [400, 100]);

	% Define time window of analysis
    pruned_segs = [-500, 3500];
    EEG = pop_loadset('filename', [num2str(subject_list(1)) '_autocleaned.set'], 'filepath', PATH_AUTOCLEANED, 'loadmode', 'info');
    tf_times = EEG.times(dsearchn(EEG.times', pruned_segs(1)) : dsearchn(EEG.times', pruned_segs(2)));

    % Save time frequency decomposition parameters
    dlmwrite([PATH_TFDECOMP 'tf_times.csv'], tf_times); 
    dlmwrite([PATH_TFDECOMP 'tf_freqs.csv'], tf_freqs); 

    % Iterating subject list
    for s = 1 : length(subject_list)

        % Load preprocessed data
        id = subject_list(s);
        EEG = pop_loadset('filename', [num2str(id) '_autocleaned.set'], 'filepath', PATH_AUTOCLEANED, 'loadmode', 'all');
        d = EEG.data;

        for ch = 1 : EEG.nbchan

            % Talk
            fprintf('\ntf decomp subject %i/%i | chan %i/%i...\n', s, numel(subject_list), ch, size(d, 1));

            % Pick channel data
            dch = squeeze(d(chansoi(ch), :, :));

            % Set convolution length
            convlen = size(dch, 1) * size(dch, 2) + size(cmw, 2) - 1;

            % cmw to freq domain and scale
            cmwX = zeros(n_freqs, convlen);
            for f = 1 : n_freqs
                cmwX(f, :) = fft(cmw(f, :), convlen);
                cmwX(f, :) = cmwX(f, :) ./ max(cmwX(f, :));
            end

            % Get TF-power
            powcube = NaN(n_freqs, size(dch, 1), size(dch, 2));
            tmp = fft(reshape(double(dch), 1, []), convlen);
            for f = 1 : n_freqs
                as = ifft(cmwX(f, :) .* tmp); 
                as = as(((size(cmw, 2) - 1) / 2) + 1 : end - ((size(cmw, 2) - 1) / 2));
                as = reshape(as, size(dch, 1), size(dch, 2));
                powcube(f, :, :) = abs(as) .^ 2;          
            end

            % Cut edge artifacts
            powcube = powcube(:, dsearchn(EEG.times', pruned_segs(1)) : dsearchn(EEG.times', pruned_segs(2)), :);

            % Save single trial power
            save([PATH_TFDECOMP num2str(id) '_powcube_chan_' num2str(chansoi(ch))], 'powcube');

        end % End channel iteration

        % Save trial information
        tmp = EEG.event(find(strcmpi({EEG.event.type}, 'cue')));
        meta = [cell2mat({tmp.id})',...
                cell2mat({tmp.session})',...
                cell2mat({tmp.watertemp})',...
                cell2mat({tmp.cuedir})',...
                cell2mat({tmp.targdir})',...
                cell2mat({tmp.cuevalid})',...
                cell2mat({tmp.accuracy})'];
        dlmwrite([PATH_TFDECOMP num2str(id) '_powcube_meta.csv'], meta); 

    end
end