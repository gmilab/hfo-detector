

clear

fcp_setup;
addpath /d/gmi/1/toolboxes/circstat/
addpath(genpath('/d/gmi/1/simeon/hfo_detector/Automatic-High-Frequency-Oscillation-Detector'))
savefolder = '/d/gmi/1/simeon/ieeg_task_analysis/analysis/GroupFinal/HFOsimulations';


% meta_path = '/d/gmi/1/rawdata/iEEG_HFO_detector/HFO_1.meta.json';
% meta_path = '/d/gmi/1/rawdata/iEEG_HFO_detector/Spiky_NoHFO_1.meta.json';
% meta_path = '/d/gmi/1/rawdata/iEEG_HFO_detector/NoHFO.meta.json';


error_log = {};


if ~exist(savefolder, 'file')
  mkdir(savefolder)
end

meta_folder = fileparts(meta_path);
metadata = fileread(meta_path);
metadata = jsondecode(metadata);

hfo_channel = metadata.hfo_ch;

montage = readtable(fullfile(meta_folder, metadata.channel_labels));



% filter montage for S/I/EEG electrodes only
montage = montage(endsWith(montage.Type, 'EEG'), :);

%% Preprocess and import
cfg = [];
cfg.dataset = fullfile(meta_folder, metadata.filename);
cfg.channel = montage.Pinbox;
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

p.preproc_cfg = cfg;

data_cnt = ft_preprocessing(cfg);

%%% overwrite custom montage labels based on our provided montage file
data_cnt.label = cellfun(@(x) montage.Label{strcmp(montage.Pinbox, x)}, data_cnt.label, 'UniformOutput', false);

% extract channels we need - pick the HFO-iest channel
cfg = [];
cfg.channel = hfo_channel;

if isempty(cfg.channel)
  error('No channels selected');
end

data_cnt = ft_preprocessing(cfg, data_cnt);


%% Detect HFOs
hfo = Core.HFO;
paras = load('/d/gmi/1/simeon/hfo_detector/Automatic-High-Frequency-Oscillation-Detector/+Demo/Spec/ECoG/Parameters/RSpecPara.mat');
hfo.Para = paras.DetPara;
hfo.Data = [];
hfo.Data.signal = data_cnt.trial{1}';
hfo.Data.channelNames = reshape(data_cnt.label, [], 1);
hfo.Data.sampFreq = data_cnt.fsample;
hfo.Data.maxIntervalToJoin = 20;
hfo.Data.MinHighEntrIntvLen = 200;
hfo.Data.minEventTime = 40;
hfo.Data.sigDurTime = data_cnt.time{1}(end);
hfo.Data.timeInterval = [1 data_cnt.time{1}(end)-1];
hfo.Data.nbChannels = length(data_cnt.label);
hfo.Data.nbSamples = length(data_cnt.trial{1});
hfo.Data.dataSetup = [];

smoothBool = false;
hfo = getFilteredSignal(hfo, smoothBool);

hfo = getBaseline(hfo);

RefType   = 'spec';
hfo = getEvents(hfo, RefType);
disp(hfo.Events)

hfo = getEventProperties(hfo);
disp(hfo.Events)

hfo = getRefinements(hfo, RefType);
disp(hfo.Refinement)

maskPassMultiChan  = hfo.Refinement.maskMultChanRefine;
CondMask = hfo.Refinement.maskEventCondSelect;
Events = hfo.Events;
hfo = hfo.RefineEvents(Events,1,CondMask, maskPassMultiChan);


%% Calculate trial hfo dichotomization by previously-saved TFR method
ch_of_interest = 1;

hfotrigs = hfo.RefinedEvents{1}.Markings.start{ch_of_interest} + (hfo.RefinedEvents{1}.Markings.len{ch_of_interest} ./ 2);
hfotrigs_s = hfotrigs / hfo.Data.sampFreq;

fprintf('%d HFOs detected\n', length(hfotrigs));



hfotrigs = round(hfotrigs);

%% Spectra compute
cfg = [];
cfg.method = 'wavelet';
cfg.output = 'pow';
cfg.pad = 'nextpow2';
cfg.channel = 1;
cfg.foi = 2:5:400;

cfg.width = 12;
cfg.gwidth = 6;

cfg.toi = data_cnt.time{1};

spec = ft_freqanalysis(cfg, data_cnt);

%% plot spectrogram with hfo markers
axs = [];

spdim = [4, 1];

hf = figure; 
hf.Position(3:4) = [493        1311];

ax = subplot(spdim(1), spdim(2), 1);
plot(data_cnt.time{1}, data_cnt.trial{1});
ax.Title.String = '1-500 Hz';
hold on
for kk = 1:length(hfotrigs)
  plot([1 1] * data_cnt.time{1}(hfotrigs(kk)), ax.YLim+[0.01 -0.01], 'k:');
end
% plot(data_cnt.time{1}, lfactive*300, '.', 'Color', 'b');
% plot(data_cnt.time{1}, hfactive*310, '.', 'Color', 'r');
% plot(data_cnt.time{1}, samples_with_hfo*320, '.', 'Color', 'y');

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

[b,a] = butter(4, 80/data_cnt.fsample*2, 'high');
ax = subplot(spdim(1), spdim(2), 2);
plot(data_cnt.time{1}, filtfilt(b,a,data_cnt.trial{1}));
ax.Title.String = '80 Hz HP';
hold on
for kk = 1:length(hfotrigs)
  plot([1 1] * data_cnt.time{1}(hfotrigs(kk)), ax.YLim+[0.01 -0.01], 'k:');
end

ax.XLabel.String = 'Time (s)';
ax.YLabel.String = 'Amplitude (\muV)';

axs(2) = ax;

[b,a] = butter(4, 120/data_cnt.fsample*2, 'high');
ax = subplot(spdim(1), spdim(2), 3);
plot(data_cnt.time{1}, filtfilt(b,a,data_cnt.trial{1}));
ax.Title.String = '120 Hz HP';
hold on
for kk = 1:length(hfotrigs)
  plot([1 1] * data_cnt.time{1}(hfotrigs(kk)), ax.YLim+[0.01 -0.01], 'k:');
end

ax.XLabel.String = 'Time (s)';
ax.YLabel.String = 'Amplitude (\muV)';

axs(3) = ax;


tfrdat = squeeze(spec.powspctrm);
tfrdat = bsxfun(@rdivide, bsxfun(@minus, tfrdat, nanmean(tfrdat, 2)), nanstd(tfrdat, [], 2));

ax = subplot(spdim(1), spdim(2), 4);
imagesc(spec.time,spec.freq, tfrdat);
% imagesc(spec.time,spec.freq, hdetect_tfrdat);
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

axs(4) = ax;



linkaxes(axs, 'x');


