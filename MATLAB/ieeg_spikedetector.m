function [spike_markers, spike_props] = ieeg_spikedetector(data, varargin)
% IEEG_SPIKEDETECTOR - Detect interictal epileptiform discharges
% Detect IEDs using algorithm from Dahal 2019
%
% spike_markers = ieeg_spikedetector(data)
% [spike_markers, spike_props] = ieeg_spikedetector(data, channel, 'param', value, ...)
%
% Parameters:
%   data - fieldtrip structure from ft_preprocessing
%
% Optional parameters, as MATLAB parameter-value pairs:
%   freq - [low high] bandpass filter range in Hz
%   zthresh - z-score threshold for detection
%   mindur - minimum duration
%   channel - channel to analyze (compatible with fieldtrip cfg.channel; Default: 1)
%
% Returns:
%   spike_markers - centroids of detected HFOs in samples
%   debug - structure containing measurements on putative events of interest
%
% References:
%  Dahal P, Ghani N, Flinker A, et al. Interictal epileptiform discharges shape large-scale 
%     intercortical communication. Brain. 2019;142(11):3502-3513. doi:10.1093/brain/awz269
% 
% 2020 Aug 1
% Simeon Wong


if ~exist('ft_freqanalysis', 'file')
  addpath(fileparts(mfilename('fullpath')), 'fieldtrip')
  ft_defaults
end

ip = inputParser;
addParameter(ip, 'filt', [25 80]);
addParameter(ip, 'zthresh', 3);
addParameter(ip, 'mindur', 0.001);
addParameter(ip, 'channel', 1);

parse(ip, varargin{:})

mindur_sample = ceil(ip.Results.mindur * data.fsample);

%% processing
cfg = [];
cfg.bpfreq = ip.Results.filt;
cfg.bpfilter = 'yes';
cfg.channel = ip.Results.channel;

data_pp = ft_preprocessing(cfg, data);


% compute envelopes
ts_pp = abs(hilbert(abs(data_pp.trial{1}(1,:))));
ts_uf = abs(hilbert(data.trial{1}(ip.Results.channel,:)));

ts_pp = zscore(ts_pp);
ts_uf = zscore(ts_uf);

%% run algo
% get segments of suprathreshold filtered signal
rp_pp = regionprops(ts_pp >= ip.Results.zthresh);


nelem = numel(rp_pp);
isspike = true(size(rp_pp,1), 1);
for kk = 1:length(rp_pp)
  % check whether it exceeds minimum duration
  if rp_pp(kk).Area < mindur_sample
    isspike(kk) = false;
    continue
  end
  
  % check whether there is a corresponding increase in broadband env
  if mean(ts_uf(round(rp_pp(kk).BoundingBox(1)) + (1:rp_pp(kk).BoundingBox(3)))) < ip.Results.zthresh
    isspike(kk) = false;
    continue
  end
end

if nelem == 0 || ~any(isspike)
  spike_markers = [];
else
  spike_markers = {rp_pp(isspike).Centroid;};
  spike_markers = cellfun(@(x) round(x(1)), spike_markers);  % * ip.Results.downsample;
end
spike_props = rp_pp(isspike);

