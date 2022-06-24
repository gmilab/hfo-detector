

clear

fcp_setup;
addpath /d/gmi/1/toolboxes/circstat/
addpath /d/gmi/1/simeon/hfo_detector/easyh5/
addpath(genpath('/d/gmi/1/simeon/hfo_detector/Automatic-High-Frequency-Oscillation-Detector'))
savefolder = '/d/gmi/1/simeon/ieeg_task_analysis/analysis/GroupFinal/HFOsimulations';


errorlog = {};

% dataset = '/d/gmi/1/simeon/hfo_detector/hfo_data_annotated/Data_Patient_01.h5';

datasetpath = '/d/gmi/1/simeon/hfo_detector/hfo_data_annotated/Intraoperative_ECoG_HFO/data';

confusion_data = table();
confusion_data.filename = cell(0,1);
confusion_data.channel = zeros(0,1);
confusion_data.TP = zeros(0,1);
confusion_data.FP = zeros(0,1);
confusion_data.FN = zeros(0,1);

cd_idx = 0;


datasets = dir(datasetpath);
datasets = [{datasets.name};];
datasets = datasets(~startsWith(datasets, '.'));

for ddd = 1:length(datasets)
  dataset = datasets{ddd};
  disp(dataset);
  
  try
    ds_path = fullfile(datasetpath, dataset);
    
    %% load data
    info_trigs = h5info(ds_path, '/data/Data_Pre_Resection_Bipolar_Montage/groups/FR markings');
    data = h5read(ds_path, '/data/Data_Pre_Resection_Bipolar_Montage/data_arrays/ECoG_Bipolar_Channels/data');
    
    device = h5read(ds_path, '/metadata/General/sections/Recording setup/properties/Recording device');
    device = strsplit(device.value{1}, '; ');
    fs = sscanf(device{3}, 'original sampling rate: %d Hz');
    
    if fs ~= 2000
      errorlog{end+1} = struct('dataset', dataset, 'text', sprintf('Sampling rate is %d', fs));
    end
    
    
    [nsamp, nch] = size(data);
    
    data_cnt = [];
    data_cnt.trial = {data'};
    data_cnt.fsample = fs;
    data_cnt.time = {0:(1/fs):((nsamp-1)/fs)};
    data_cnt.label = arrayfun(@(x)['C',num2str(x)], 1:nch, 'UniformOutput', false);
    
    
    if isempty(info_trigs.Groups)
      ntrigs = 0;
    else
      ntrigs = length(info_trigs.Groups.Groups);
    end
    
    trigs = zeros(ntrigs,2);
    
    for kk = 1:ntrigs
      trigattr = struct2table(info_trigs.Groups.Groups(kk).Attributes);
      ch_num = trigattr.Value{strcmp(trigattr.Name, 'name')};
      ch_num = strsplit(ch_num, '_');
      ch_num = str2double(ch_num{4});
      
      trig_loc = h5read(ds_path, [info_trigs.Groups.Groups(kk).Name, '/position']);
      
      trigs(kk,1) = round(ch_num);
      trigs(kk,2) = trig_loc(2);
    end
    
    [~,i] = sort(trigs(:,2));
    trigs = trigs(i,:);
    
    
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
    
    for ch_of_interest = 1:nch
        gs_trigs_ch = trigs(trigs(:,1) == ch_of_interest, 2);
        
        % extract detected HFOs
        hfotrigs = hfo.RefinedEvents{1}.Markings.start{ch_of_interest} + (hfo.RefinedEvents{1}.Markings.len{ch_of_interest} ./ 2);
        hfotrigs_s = hfotrigs / hfo.Data.sampFreq;
        
        fprintf('%d HFOs detected / %d HFOs ground truth\n', length(hfotrigs_s), sum(trigs(:,1) == ch_of_interest));        
        
        %% generate confusion matrix
        cd_idx = cd_idx + 1;
        confusion_data.filename{cd_idx} = dataset;
        confusion_data.channel(cd_idx) = ch_of_interest;
        
        for nn = 1:length(hfotrigs_s)
          if any(abs(gs_trigs_ch - hfotrigs_s(nn)) < 2)
            confusion_data.TP(cd_idx,1) = confusion_data.TP(cd_idx,1) + 1;
          else
            confusion_data.FP(cd_idx,1) = confusion_data.FP(cd_idx,1) + 1;
          end
        end
        
        for nn = 1:length(gs_trigs_ch)
          if ~any(abs(gs_trigs_ch(nn) - hfotrigs_s) < 2)
            confusion_data.FN(cd_idx,1) = confusion_data.FN(cd_idx,1) + 1;
          end
        end

    end
    
  catch ex
    errorlog{end+1} = struct('dataset', dataset, 'err', ex);
  end
  
end

%% Generate report
sum_TP = sum(confusion_data.TP);
sum_FP = sum(confusion_data.FP);
sum_FN = sum(confusion_data.FN);
fprintf('\n\n\nTP: %d\tFP: %d\tFN:%d\n', sum_TP, sum_FP, sum_FN);
fprintf('PPV: %.5f\tSensitivity: %.5f\n\n\n', sum_TP/(sum_TP+sum_FP), sum_TP/(sum_TP+sum_FN))


%% Error reporting

for ee = 1:length(errorlog)
  disp('===================')
  disp(errorlog{ee}.dataset)
  if isfield(errorlog{ee}, 'ch')
    disp(['Channel ' num2str(errorlog{ee}.ch)])
  end
  if isfield(errorlog{ee}, 'err')
    disp(getReport(errorlog{ee}.err))
  end
  if isfield(errorlog{ee}, 'text')
    disp(errorlog{ee}.text)
  end
end


return

%% Spectra compute
% cfg = [];
% cfg.method = 'wavelet';
% cfg.output = 'pow';
% cfg.pad = 'nextpow2';
% cfg.channel = 1;
% cfg.foi = 2:5:400;
%
% cfg.width = 12;
% cfg.gwidth = 6;
%
% cfg.toi = data_cnt.time{1};
%
% spec = ft_freqanalysis(cfg, data_cnt);


%% Spectra plot
% axs = [];
% spdim = [3, 1];
%
% hf = figure;
% hf.Position(3:4) = [493        1311];
%
% ax = subplot(spdim(1), spdim(2), 1);
% plot(data_cnt.time{1}, data_cnt.trial{1});
% ax.Title.String = '1-500 Hz';
% ax.XLabel.String = 'Time (s)';
% ax.YLabel.String = 'Amplitude (\muV)';
%
% ax.XTick = 0:data_cnt.time{1}(end);
% % ax.XMinorTick = 0:0.2:data_cnt.time{1}(end);
%
% ax.GridColor = [98, 227, 130]/255;
% ax.GridAlpha = 1;
% ax.MinorGridColor = [98, 227, 130]/255;
% ax.MinorGridAlpha = 1;
% ax.MinorGridLineStyle = ':';
% ax.XGrid = 'on';
% ax.XMinorGrid = 'on';
%
%
% axs(1) = ax;
%
%
% % ax = subplot(1, 2, 2);
% % [s,f,t] = spectrogram(data_cnt.trial{1}, 512, 500, 2048, 2048);
% % imagesc(t,f,log10(abs(s)));
% % axis xy
% % ax.Title.String = 'Spectrogram';
% % ax.YLim = [0 400];
% % axs(2) = ax;
%
%
% % ax = subplot(1, 3, 2);
% % imagesc(spec.time,spec.freq, log10(squeeze(spec.powspctrm)));
% % axis xy
% % ax.Title.String = 'Spectrogram';
% % ax.YLim = [0 400];
% % axs(2) = ax;
%
% [b,a] = butter(4, [80 400]/data_cnt.fsample*2);
% ax = subplot(spdim(1), spdim(2), 2);
% plot(data_cnt.time{1}, filtfilt(b,a,data_cnt.trial{1}));
% ax.Title.String = '120-300 Hz';
% ax.XLabel.String = 'Time (s)';
% ax.YLabel.String = 'Amplitude (\muV)';
%
%
% ax.XTick = 0:data_cnt.time{1}(end);
% % ax.XMinorTick = 0:0.2:data_cnt.time{1}(end);
%
% ax.GridColor = [98, 227, 130]/255;
% ax.GridAlpha = 1;
% ax.MinorGridColor = [98, 227, 130]/255;
% ax.MinorGridAlpha = 1;
% ax.MinorGridLineStyle = ':';
% ax.XGrid = 'on';
% ax.XMinorGrid = 'on';
%
%
% axs(2) = ax;
%
%
% tfrdat = squeeze(spec.powspctrm);
% tfrdat = bsxfun(@rdivide, bsxfun(@minus, tfrdat, nanmean(tfrdat, 2)), nanstd(tfrdat, [], 2));
%
% ax = subplot(spdim(1), spdim(2), 3);
% imagesc(spec.time,spec.freq, tfrdat);
% axis xy
% ax.Title.String = 'Spectrogram';
% ax.YLim = [0 400];
% ax.CLim = [0 3];
%
% ax.XLabel.String = 'Time (s)';
% ax.YLabel.String = 'Frequency (Hz)';
%
% tmp = ax.Position;
% cb = colorbar('Location', 'SouthOutside');
% ax.Position = tmp;
%
% cb.Label.String = 'Amplitude (z-score)';
%
% axs(3) = ax;
%
%
%
% linkaxes(axs, 'x');
%



return

%% plot spectrogram with hfo markers
lfactive = single(lfactive);
lfactive(lfactive == 0) = NaN;

hfactive = single(hfactive);
hfactive(hfactive == 0) = NaN;

samples_with_hfo = single(samples_with_hfo);
samples_with_hfo(samples_with_hfo == 0) = NaN;


axs = [];

spdim = [3, 1];

hf = figure;
hf.Position(3:4) = [493        1311];

ax = subplot(spdim(1), spdim(2), 1);
plot(data_cnt.time{1}, data_cnt.trial{1});
ax.Title.String = '1-500 Hz';
hold on
for kk = 1:length(hfotrigs)
  plot([1 1] * data_cnt.time{1}(hfotrigs(kk)*6), ax.YLim+[0.01 -0.01], 'k:');
end
plot(data_cnt.time{1}(1:6:end), lfactive*300, '.', 'Color', 'b');
plot(data_cnt.time{1}(1:6:end), hfactive*310, '.', 'Color', 'r');
plot(data_cnt.time{1}(1:6:end), samples_with_hfo*320, '.', 'Color', 'y');

cchtrigs = trigs(trigs(:,1) == ch_of_interest, 2);
for kk = 1:length(cchtrigs)
  plot([1 1] * cchtrigs(kk), ax.YLim+[0.01 -0.01], 'r:');
end

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

[b,a] = butter(4, 120/data_cnt.fsample*2, 'high');
ax = subplot(spdim(1), spdim(2), 2);
plot(data_cnt.time{1}, filtfilt(b,a,data_cnt.trial{1}));
ax.Title.String = '120 Hz HP';
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
% imagesc(spec.time,spec.freq, tfrdat);
imagesc(spec.time,spec.freq, hdetect_tfrdat);
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