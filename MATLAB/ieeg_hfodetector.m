function [hfo_markers,tfrdat_out,debug] = ieeg_hfodetector(data, varargin)
% IEEG_HFODETECTOR - Detect HFOs
% Detect HFOs in time-frequency space by defining characteristics of blobs
%
% hfo_markers = ieeg_hfodetector(data, channel)
% [hfo_markers, tfrdat_out, debug] = ieeg_hfodetector(data, channel, 'param', value, ...)
%
% Parameters:
%   data - fieldtrip structure from ft_preprocessing
%
% Optional parameters, as MATLAB parameter-value pairs:
%   channel - channel to analyze (compatible with fieldtrip cfg.channel)
%             Default: 1
%
%   foi - frequency range of interest in Hz
%         Default: [80 500]
%
%   thresh - power ratio threshold relative to background
%            Default: 5.75
%
%   ncycles - minimum length of HFO event in number of cycles
%             Default: 4.75
%
%   hfo_freq_width - contrain frequency bandwidth of the event in Hz
%                    Default: [20 150]
%
%   wavelet_width - number of cycles in Morlet wavelet
%   downsample - compute the transform every n samples
%   freqstep - compute the transform every n Hz
%
% Returns:
%   hfo_markers - centroids of detected HFOs in samples
%   tfrdat_out - time-frequency data from wavelet transform
%   debug - structure containing measurements on putative events of interest
%
% References:
%   1. Wong SM, Arski ON, Workewych AM, et al. Detection of high-frequency oscillations in electroencephalography: 
%      A scoping review and an adaptable open-source framework. Seizure. 2021;84:23-33. doi:10.1016/j.seizure.2020.11.009
% 
% 2020 Aug 1
% Simeon Wong

if ~exist('ft_freqanalysis', 'file')
  addpath(fileparts(mfilename('fullpath')), 'fieldtrip')
  ft_defaults
end

ip = inputParser;
addParameter(ip, 'foi', [80 500]);
addParameter(ip, 'lf_fthresh', []);
addParameter(ip, 'thresh', 5.75);
addParameter(ip, 'ncycles', 4.75);
addParameter(ip, 'wavelet_width', 9);
addParameter(ip, 'exclude_top', 10);
addParameter(ip, 'exclude_bottom', 10);
addParameter(ip, 'circularity', []);
addParameter(ip, 'hfo_max_freq_width', 150);
addParameter(ip, 'hfo_min_freq_width', 20);

addParameter(ip, 'downsample', 1);
addParameter(ip, 'freqstep', 5);

addParameter(ip, 'channel', 1);

parse(ip, varargin{:});

cfg = [];
cfg.method = 'wavelet';
cfg.output = 'pow';
cfg.pad = 'nextpow2';
cfg.channel = ip.Results.channel;
cfg.foi = 2:ip.Results.freqstep:ip.Results.foi(2);
cfg.width = ip.Results.wavelet_width;
% cfg.gwidth = 5;
cfg.toi = data.time{1}(1:ip.Results.downsample:end);

spec = ft_freqanalysis(cfg, data);


get_f_idx = @(x) min(abs(x - spec.freq));

% identify freq idx
fq_idx = zeros(1,2);
for ff = 1:2
  [~,fq_idx(ff)] = get_f_idx(ip.Results.foi(ff));
end

if ~isempty(ip.Results.lf_fthresh)
  [~, lf_idx] = get_f_idx(ip.Results.lf_fthresh);
end


% mean over channels
tfrdat_raw = squeeze(mean(spec.powspctrm, 1));
% tfrdat is [freq] x [time]

% exclude top n-th percentile from mean (deal with artifacts)
tfrdat_raw_mean = tfrdat_raw;
tfrdat_raw_mean(tfrdat_raw > prctile(tfrdat_raw, 100-ip.Results.exclude_top, 2)) = NaN;
tfrdat_raw_mean(tfrdat_raw < prctile(tfrdat_raw, ip.Results.exclude_bottom, 2)) = NaN;

% power ratio
tfrdat = bsxfun(@rdivide, tfrdat_raw, nanmean(tfrdat_raw_mean, 2)) - 1;

% save for output
tfrdat_out = tfrdat;


% threshold by cfg.zthresh
tfrdat = tfrdat >= ip.Results.thresh;

tfrprop = regionprops(tfrdat, 'Centroid', 'BoundingBox', 'PixelIdxList', 'Circularity');

% loop through and check
%  - centroid within HFO region
%  - height proportional to width?
%  - is not contiguous with lf activity
%  - minlength

nelem = length(tfrprop);
ishfo = true(nelem,1);

debug.tfrprop = tfrprop;
debug.centroid_position = zeros(nelem,1);
debug.bounds_in_lf = zeros(nelem,1);
debug.length_insufficient = zeros(nelem,1);
debug.circularity = zeros(nelem,1);
debug.hfo_freq_width = zeros(nelem,1);

for kk = 1:nelem
  % check centroid location
  debug.centroid_position(kk) = spec.freq(round(tfrprop(kk).Centroid(2)));
  if tfrprop(kk).Centroid(2) < fq_idx(1) || tfrprop(kk).Centroid(2) > fq_idx(2)
    ishfo(kk) = false;
  end
  
  % check edges don't touch LF
  if ~isempty(ip.Results.lf_fthresh) && tfrprop(kk).BoundingBox(2) <= lf_idx
    debug.bounds_in_lf(kk) = true;
    ishfo(kk) = false;
  end
  
  % check number of cycles
  debug.length_insufficient(kk) = tfrprop(kk).BoundingBox(3) / data.fsample * spec.freq(round(tfrprop(kk).Centroid(2)));
  
  len_samples = data.fsample * ip.Results.ncycles / spec.freq(round(tfrprop(kk).Centroid(2)));
  if (tfrprop(kk).BoundingBox(3) < len_samples)
    ishfo(kk) = false;
  end
  
  % check circularity
  %  - spikes have jagged morphology in HF TFD
  debug.circularity(kk) = tfrprop(kk).Circularity;
  if ~isempty(ip.Results.circularity) && ...
    (tfrprop(kk).Circularity < ip.Results.circularity)
    ishfo(kk) = false;
  end
  
  % check frequency width
  %  - spikes have high width in the frequency domain
  %  - spurious events have low width
  
  if ~isempty(ip.Results.hfo_max_freq_width)
    % get bandwidth in the frequency domain
    freqwidth = mean(sum(tfrdat(floor(tfrprop(kk).BoundingBox(2))+(1:ceil(tfrprop(kk).BoundingBox(4))), ...
      floor(tfrprop(kk).BoundingBox(1))+(1:ceil(tfrprop(kk).BoundingBox(3)))), 1), 2) * ip.Results.freqstep;
    
    debug.hfo_freq_width(kk) = freqwidth;
    
    if freqwidth < ip.Results.hfo_freq_width(1)
      ishfo(kk) = false;
    end
    
    if freqwidth > ip.Results.hfo_freq_width(2)
      ishfo(kk) = false;
    end
  end
end

debug.ishfo = ishfo;
debug.spec = spec;
debug.ip = ip;

if nelem == 0 || ~any(ishfo)
  hfo_markers = [];
else
  hfo_markers = {tfrprop(ishfo).Centroid;};
  hfo_markers = cellfun(@(x) round(x(1)), hfo_markers) * ip.Results.downsample;
end

end