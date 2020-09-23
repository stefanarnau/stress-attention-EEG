function[cmw, frqs] = get_cmwset(range_frq, n_frq, spacing, srate, tfres_dim, tfres_range)

	% Set desired range
	range_fwhmT = tfres_range; % tfres_dim 'time'
	range_fwhmF = tfres_range; % tfres_dim 'freq'
	range_cycle = tfres_range; % tfres_dim 'cycle'
	newtimef_params = tfres_range; % tfres_dim 'eeglab'

	% Set wavelet time
	wtime = -2 : 1 / srate : 2;

	% Determine fft frqs
	hz = linspace(0, srate, length(wtime));

	% Define wavelet frequencies
	if strcmpi(spacing, 'lin')
		frqs = linspace(range_frq(1), range_frq(2), n_frq);
		fwhmTs = linspace(range_fwhmT(1), range_fwhmT(2), n_frq);
		fwhmFs = linspace(range_fwhmF(1), range_fwhmF(2), n_frq);
		n_cycles = linspace(range_cycle(1), range_cycle(2), n_frq);
		n_cycles_eeglab = linspace(newtimef_params(1), (1 - newtimef_params(2)) * range_frq(2), n_frq);
	elseif strcmpi(spacing, 'log')
		frqs = logspace(log10(range_frq(1)), log10(range_frq(2)), n_frq);
		fwhmTs = logspace(log10(range_fwhmT(1)), log10(range_fwhmT(2)), n_frq);
		fwhmFs = logspace(log10(range_fwhmF(1)), log10(range_fwhmF(2)), n_frq);
		n_cycles = logspace(log10(range_cycle(1)), log10(range_cycle(2)), n_frq);
		n_cycles_eeglab = logspace(log10(newtimef_params(1)), log10((1 - newtimef_params(2)) * range_frq(2)), n_frq);
	end

	% Set cycles as EEGlab would
	if strcmpi(tfres_dim, 'eeglab')
		n_cycles = n_cycles_eeglab;
		tfres_dim = 'cycle';
	end

	% Init matrices for wavelets
	cmw = zeros(length(frqs), length(wtime));
	cmwX = zeros(length(frqs), length(wtime));

	% Iterate freqs
	for frq = 1 : length(frqs)

		% Create tapering gaussian corresponding to desired freqwidth
		if strcmpi(tfres_dim, 'freq')

			% Create frequency domain gaussian
			s = fwhmFs(frq) * (2 * pi - 1) / (4 * pi);
			x = hz - frqs(frq);
			fx = exp(-0.5 * (x / s).^2);

			% Create time domain wavelet
			cmw(frq, :) = fftshift(ifft(fx));

		% Create tapering gaussian corresponding to desired timewidth
		elseif strcmpi(tfres_dim, 'time')

			cmw(frq, :) = exp(2 * 1i * pi * frqs(frq) .* wtime) .* exp((-4 * log(2) * wtime.^2) ./ (fwhmTs(frq) / 1000)^2);

		% Create tapering gaussian corresponding to desired number of cycles
		elseif strcmpi(tfres_dim, 'cycle')

			s2 = n_cycles(frq) / (2 * pi * frqs(frq));

	        cmw(frq, :) = exp(2 * 1i * pi * frqs(frq) .* wtime) .* exp((-wtime .^ 2) / (2 * s2 .^2));

		end

		% Normalize wavelet
		cmw(frq, :) = cmw(frq, :) ./ max(cmw(frq, :));

		% Create normalized freq domain wavelet
		cmwX(frq, :) = fft(cmw(frq, :)) ./ max(fft(cmw(frq, :)));

	end
end