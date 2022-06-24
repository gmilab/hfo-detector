

clear

fcp_setup;
addpath /d/gmi/1/simeon/hsc-eeg-tools/MATLAB/

savefolder = '/d/gmi/1/simeon/ieeg_task_analysis/analysis/GroupFinal/HFOsimulations';


% meta_path = '/d/gmi/1/rawdata/iEEG_HFO_detector/HFO_1.meta.json';
% meta_path = '/d/gmi/1/rawdata/iEEG_HFO_detector/Spiky_NoHFO_1.meta.json';
meta_path = '/d/gmi/1/rawdata/iEEG_HFO_detector/NoHFO.meta.json';


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
% montage = montage(endsWith(montage.Type, 'EEG'), :);
montage = montage(strcmpi(montage.Type, 'SEEG') | strcmpi(montage.Type, 'IEEG'), :);

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


%% Detect
[hfotrigs,tfrdat,hfodebug] = ieeg_hfodetector_boundingbox(data_cnt, 1, 'circularity', 0.5, 'hfo_max_freq_width', 150);


  [~,ia] = unique(round(hfotrigs/500));
  hfotrigs = hfotrigs(ia);


%% plot
hf = ieeg_plotHFOdiagnostics(data_cnt.time{1}, data_cnt.trial{1}(1,:), ...
  hfotrigs, tfrdat, hfodebug, [], ...
  'showTFRinfo', true);


%%

hfotrigs_s = data_cnt.time{1}(hfotrigs);

wnds = 0:10:data_cnt.time{1}(end);
wnds = [reshape(wnds(1:end-1), [], 1), reshape(wnds(2:end), [], 1)];

counts = zeros(size(wnds,1),1);

for kk = 1:size(wnds,1)
  counts(kk) = sum((hfotrigs_s >= wnds(kk,1)) & (hfotrigs_s <= wnds(kk,2)));
end


fprintf('HFOs / min: %.1f±%.1f\n', mean(counts)/10*60, std(counts)/10*60);

