
% validate hfo detector using simulated HFOs

fcp_setup;
addpath /sdata4/toolboxes/circstat/
savefolder = '/sdata4/simeon/ieeg_task_analysis/analysis/GroupFinal/HFOsimulations';


error_log = {};


if ~exist(savefolder, 'file')
  mkdir(savefolder)
end

p = ieeg_config('SEEG-SK-07.1', 'Rest');

%% Preprocess and import
cfg = [];
cfg.dataset = p.dataset;
cfg.channel = p.allchannels;
cfg.bpfilter = 'yes';
cfg.bpfreq = [1 500];
cfg.dftfilter = 'yes';
cfg.dftfreq = 60 * [1 2 3 4 5 6 7 8 9 10];
cfg.dftreplace = 'neighbour';
cfg.dftneighbourwidth = [0.5 1 2 4 4 4 4 4 4 4];
cfg.dftbandwidth = [1 2 4 4 4 4 4 4 4 4];
cfg.continuous = 'yes';
cfg.reref = 'yes';
cfg.refmethod = 'avg';
cfg.refchannel = 'all';
cfg.demean = 'yes';

if isfield(p, 'trigsampfilter')
  cfg.trl = [p.trigsampfilter(1)-10000, p.trigsampfilter(2)+10000, -10000];
end

p.preproc_cfg = cfg;

data_cnt = ft_preprocessing(cfg);
data_cnt.sampleinfo = data_cnt.sampleinfo - data_cnt.sampleinfo(1);


to_samp = @(x) round((x-p.epoch_cfg.trialdef.tEpoch(1)) * 2048) + 1;
to_sec = @(x) (x-1) / 2048 + p.epoch_cfg.trialdef.tEpoch(1);


%%% overwrite custom montage labels because recorded file has them wrong
data_cnt.label = p.montage;

% extract channels we need - pick a cortical channel, relatively clean
cfg = [];
% cfg.channel = find(strcmp(p.montage, 'LOAO-9'));   % SK-06
cfg.channel = find(strcmp(p.montage, 'RF3OF-9'));   % SK-07

if isempty(cfg.channel)
  error('No channels selected');
end

data_cnt = ft_preprocessing(cfg, data_cnt);


%% raw data plot to verify suitability of source timeseries
axs = [];

figure;
ax = subplot(3,1,1);
plot(data_cnt.time{1}, data_cnt.trial{1});
ax.Title.String = 'Raw';
axs(1) = ax;

[b,a] = butter(4, [1 70]/data_cnt.fsample*2);
ax = subplot(3,1,2);
plot(data_cnt.time{1}, filtfilt(b,a,data_cnt.trial{1}));
ax.Title.String = '1-70Hz';
axs(2) = ax;

[b,a] = butter(4, [120 300]/data_cnt.fsample*2);
ax = subplot(3,1,3);
plot(data_cnt.time{1}, filtfilt(b,a,data_cnt.trial{1}));
ax.Title.String = '120-300Hz';
axs(3) = ax;

linkaxes(axs);


%% Hfo
% Calculate trial hfo dichotomization by previously-saved TFR method
[samples_with_hfo,lfactive,hfactive,hdetect_tfrdat] = ieeg_hfodetector(data_cnt, 1, 'zthresh', 3, 'lf_foi', [1 50], 'foi', [150 500]);

hfotrigs = bwconncomp(samples_with_hfo);
hfotrigs = reshape(cellfun(@(x) round(median(x)), hfotrigs.PixelIdxList), [], 1);

lfactive = single(lfactive);
lfactive(lfactive == 0) = NaN;

hfactive = single(hfactive);
hfactive(hfactive == 0) = NaN;

samples_with_hfo = single(samples_with_hfo);
samples_with_hfo(samples_with_hfo == 0) = NaN;

fprintf('%d HFOs detected\n', length(hfotrigs));

%% compute induced HFO validation
cfg = [];
cfg.trl = [(hfotrigs + repmat([-1024 1024], length(hfotrigs), 1)), -1025 * ones(length(hfotrigs))];

hfodat = ft_redefinetrial(cfg, data_cnt);

cfg = [];
cfg.output = 'pow';
cfg.toi = -0.4:0.01:0.4;
cfg.channel = 'all';
cfg.foi = 10:2:500;

%%%
cfg.method = 'mtmconvol';
cfg.tapsmofrq = 25;
cfg.t_ftimwin = 0.05 * ones(size(cfg.foi));
%%%

%%%
%       cfg.method = 'wavelet';
%%%

hfotfr = ft_freqanalysis(cfg, hfodat);

hfotfr.label = {'Lesion', 'ACC'};

cfg = []; cfg.baseline = 'yes'; cfg.baselinetype = 'relative'; cfg.channel = 1;
hf = figure;
ft_singleplotTFR(cfg, hfotfr);



%% compute spectrogram
cfg = [];
cfg.method = 'wavelet';
cfg.output = 'pow';
cfg.pad = 'nextpow2';
cfg.channel = 1;
cfg.foi = 2:5:400;
cfg.toi = data_cnt.time{1};

spec = ft_freqanalysis(cfg, data_cnt);



%% plot spectrogram
axs = [];

spdim = [3, 1];

hf = figure; 
hf.Position(3:4) = [493        1311];

ax = subplot(spdim(1), spdim(2), 1);
plot(data_cnt.time{1}, data_cnt.trial{1});
ax.Title.String = '1-500 Hz';
hold on
for kk = 1:length(hfotrigs)
  plot([1 1] * data_cnt.time{1}(hfotrigs(kk)), ax.YLim+[0.01 -0.01], 'k:');
end
plot(data_cnt.time{1}, lfactive*300, '.', 'Color', 'b');
plot(data_cnt.time{1}, hfactive*310, '.', 'Color', 'r');
plot(data_cnt.time{1}, samples_with_hfo*320, '.', 'Color', 'y');

ax.XLabel.String = 'Time (s)';
ax.YLabel.String = 'Amplitude (\muV)';


axs(1) = ax;


% ax = subplot(1, 2, 2);
% [s,f,t] = spectrogram(data_cnt.trial{1}, 512, 500, 2048, 2048);
% imagesc(t,f,log10(abs(s)));
% axis xy
% ax.Title.String = 'Spectrogram';
% ax.YLim = [0 400];
% axs(2) = ax;


% ax = subplot(1, 3, 2);
% imagesc(spec.time,spec.freq, log10(squeeze(spec.powspctrm)));
% axis xy
% ax.Title.String = 'Spectrogram';
% ax.YLim = [0 400];
% axs(2) = ax;

[b,a] = butter(4, [80 400]/data_cnt.fsample*2);
ax = subplot(spdim(1), spdim(2), 2);
plot(data_cnt.time{1}, filtfilt(b,a,data_cnt.trial{1}));
ax.Title.String = '120-300 Hz';
hold on
for kk = 1:length(hfotrigs)
  plot([1 1] * data_cnt.time{1}(hfotrigs(kk)), ax.YLim+[0.01 -0.01], 'k:');
end

ax.XLabel.String = 'Time (s)';
ax.YLabel.String = 'Amplitude (\muV)';

axs(2) = ax;


tfrdat = squeeze(spec.powspctrm);
tfrdat = bsxfun(@rdivide, bsxfun(@minus, tfrdat, nanmean(tfrdat, 2)), nanstd(tfrdat, [], 2));

ax = subplot(spdim(1), spdim(2), 3);
imagesc(spec.time,spec.freq, tfrdat);
axis xy
ax.Title.String = 'Spectrogram';
ax.YLim = [0 400];
ax.CLim = [0 3];
hold on
for kk = 1:length(hfotrigs)
  plot([1 1] * data_cnt.time{1}(hfotrigs(kk)), ax.YLim+[0.01 -0.01], 'w:');
end

ax.XLabel.String = 'Time (s)';
ax.YLabel.String = 'Frequency (Hz)';

tmp = ax.Position;
cb = colorbar('Location', 'SouthOutside');
ax.Position = tmp;

cb.Label.String = 'Amplitude (z-score)';

axs(3) = ax;



linkaxes(axs, 'x');



%% save
[dir_file, dir_path] = uiputfile('/d/mjt/s4/simeon/ieeg_task_analysis/analysis/HFO/*.png');
if dir_file ~= 0
  print(hf, fullfile(dir_path, dir_file), '-dpng', '-r600');
  disp done
else
  disp cancelled.
end





  
  %% error log
  for ee = 1:length(error_log)
    fprintf('\n\n -- %s -- \n\n', error_log{ee}.subj)
    disp(getReport(error_log{ee}.error));
  end
  
