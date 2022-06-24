function [spindles_markers, spindles_prop] = ieeg_spindledetector(data, varargin)
% IEEG_SPINDLEDETECTOR - Detect spindles
% Detect spindles using algorithm from Dahal 2019
%
% spindles_markers = ieeg_spindledetector(data)
% [spindles_markers, spindles_prop] = ieeg_spindledetector(data, channel, 'param', value, ...)
%
% Parameters:
%   data - fieldtrip structure from ft_preprocessing
%
% Optional parameters, as MATLAB parameter-value pairs:
%   filt_low - lower freq range bandpass in Hz
%   filt_high - upper freq range bandpass in Hz
%   filt_spindle - spindle range bandpass in Hz
%
%   dur_min - minimum duration
%   dur_max - minimum duration
%  
%   threshold - frequency ratio threshold
%
%   channel - channel to analyze (compatible with fieldtrip cfg.channel; Default: 1)
%
%   downsample - compute frequency decomposition ever N samples
%
% Returns:
%   spindles_markers - centroids of detected HFOs in samples
%   spindles_prop - structure containing measurements on putative events of interest
%
% References:
%   Dahal P, Ghani N, Flinker A, et al. Interictal epileptiform discharges shape large-scale 
%     intercortical communication. Brain. 2019;142(11):3502-3513. doi:10.1093/brain/awz269
% 
% 2020 Aug 1
% Simeon Wong


if ~exist('ft_freqanalysis', 'file')
  addpath(fileparts(mfilename('fullpath')), 'fieldtrip')
  ft_defaults
end

ip = inputParser;
addParameter(ip, 'filt_low', [2 8]);
addParameter(ip, 'filt_spindle', [10 20]);
addParameter(ip, 'filt_high', [25 40]);
addParameter(ip, 'dur_min', 0.3);
addParameter(ip, 'dur_max', 3);
addParameter(ip, 'threshold', 0.1);
addParameter(ip, 'channel', 1);
addParameter(ip, 'downsample', 4);

parse(ip, varargin{:})

mindur_sample = ceil(ip.Results.dur_min * data.fsample);
maxdur_sample = ceil(ip.Results.dur_max * data.fsample);

%% processing
cfg = [];
cfg.method = 'wavelet';
cfg.channel = ip.Results.channel;
cfg.toi = data.time{1}(1:ip.Results.downsample:end);

cfg.foi = mean(ip.Results.filt_low);
cfg.width = 6*cfg.foi/diff(ip.Results.filt_low);
data_low = ft_freqanalysis(cfg, data);

cfg.foi = mean(ip.Results.filt_spindle);
cfg.width = 6*cfg.foi/diff(ip.Results.filt_spindle);
data_spindle = ft_freqanalysis(cfg, data);

cfg.foi = mean(ip.Results.filt_high);
cfg.width = 6*cfg.foi/diff(ip.Results.filt_high);
data_high = ft_freqanalysis(cfg, data);

ts_low = squeeze(data_low.powspctrm);
ts_spindle = squeeze(data_spindle.powspctrm);
ts_high = squeeze(data_high.powspctrm);


%% run algo
ratio = (ts_spindle - ts_low - ts_high)./(ts_spindle + ts_low + ts_high);

% get segments of suprathreshold filtered signal
rp_pp = regionprops(ratio >= ip.Results.threshold);

nelem = numel(rp_pp);
isspindle = true(size(rp_pp,1), 1);
for kk = 1:length(rp_pp)
  % check whether it passes duration criterion
  if rp_pp(kk).Area < mindur_sample || rp_pp(kk).Area > maxdur_sample
    isspindle(kk) = false;
    continue
  end
end

if nelem == 0 || ~any(isspindle)
  spindles_markers = [];
else
  spindles_markers = {rp_pp(isspindle).Centroid;};
  spindles_markers = cellfun(@(x) round(x(1)), spindles_markers) * ip.Results.downsample;
end
spindles_prop = rp_pp(isspindle);

