

clear

fcp_setup;
addpath /d/gmi/1/toolboxes/circstat/
addpath /d/gmi/1/simeon/hfo_detector/easyh5/
savefolder = '/d/gmi/1/simeon/ieeg_task_analysis/analysis/GroupFinal/HFOsimulations';


errorlog = {};

% dataset = '/d/gmi/1/simeon/hfo_detector/hfo_data_annotated/Data_Patient_01.h5';

datasetpath = '/d/gmi/1/simeon/hfo_detector/hfo_data_annotated/Intraoperative_ECoG_HFO/data';

datasets = dir(datasetpath);
datasets = [{datasets.name};];
datasets = datasets(~startsWith(datasets, '.'));


varparams = 5:5:40;

for idx_varparam = 1:length(varparams)
  prm = varparams(idx_varparam);
  
  confusion_data = table();
  confusion_data.filename = cell(0,1);
  confusion_data.channel = zeros(0,1);
  confusion_data.TP = zeros(0,1);
  confusion_data.FP = zeros(0,1);
  confusion_data.FN = zeros(0,1);
  
  cd_idx = 0;
  
  
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
      
      
      %% channel loop
      for ch_of_interest = 1:nch
        try
          gs_trigs_ch = trigs(trigs(:,1) == ch_of_interest, 2);
          
          %% Calculate trial hfo dichotomization by previously-saved TFR method
          % === SICKKIDS DETECTOR ===
          [hfotrigs_sk,tfrdat,hfodet_ip] = ieeg_hfodetector_boundingbox(data_cnt, ch_of_interest);
          hfotrigs_sk_s = ((hfotrigs_sk-1)*hfodet_ip.Results.downsample+1)/fs;
          
          fprintf('%d HFOs detected / %d HFOs ground truth\n', length(hfotrigs_sk), sum(trigs(:,1) == ch_of_interest));
          
          % === ZURICH DETECTOR ===
          hfotrigs_z = hfo.RefinedEvents{1}.Markings.start{ch_of_interest} + (hfo.RefinedEvents{1}.Markings.len{ch_of_interest} ./ 2);
          hfotrigs_z_s = hfotrigs_z / hfo.Data.sampFreq;
          
          fprintf('%d HFOs detected / %d HFOs ground truth\n', length(hfotrigs_z), sum(trigs(:,1) == ch_of_interest));
          
        catch ex
          errorlog{end+1} = struct('dataset', dataset, 'ch', ch_of_interest, 'err', ex);
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
  
  writetable(confusion_data, sprintf('/d/gmi/1/simeon/ieeg_task_analysis/analysis/HFO-ROC/us/minlen-%d.csv', prm))
  
end


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

%% plot spectrogram with hfo markers
hfotrigs = round(hfotrigs);

axs = [];

spdim = [4, 1];

hf = figure;
hf.Position = [20 20 1600 900];

ax = subplot(spdim(1), spdim(2), 1);
plot(data_cnt.time{1}, data_cnt.trial{1}(ch_of_interest,:));
ax.Title.String = '1-500 Hz';
hold on
for kk = 1:length(hfotrigs_sk_s)
  plot([1 1] * hfotrigs_sk_s(kk), ax.YLim+[0.01 -0.01], 'k--');
end
for kk = 1:length(hfotrigs_z_s)
  plot([1 1] * hfotrigs_z_s(kk), ax.YLim+[0.01 -0.01], 'b--');
end
for kk = 1:length(gs_trigs_ch)
  plot([1 1] * gs_trigs_ch(kk), ax.YLim+[0.01 -0.01], '--', 'Color', '#eb9e34');
end


ax.XLabel.String = 'Time (s)';
ax.YLabel.String = 'Amplitude (\muV)';


axs(1) = ax;



[b,a] = butter(4, 80/data_cnt.fsample*2, 'high');
ax = subplot(spdim(1), spdim(2), 2);
plot(data_cnt.time{1}, filtfilt(b,a,data_cnt.trial{1}(ch_of_interest,:)));
ax.Title.String = '80 Hz HP';
hold on
for kk = 1:length(hfotrigs_sk_s)
  plot([1 1] * hfotrigs_sk_s(kk), ax.YLim+[0.01 -0.01], 'k--');
end
for kk = 1:length(hfotrigs_z_s)
  plot([1 1] * hfotrigs_z_s(kk), ax.YLim+[0.01 -0.01], 'b--');
end
for kk = 1:length(gs_trigs_ch)
  plot([1 1] * gs_trigs_ch(kk), ax.YLim+[0.01 -0.01], '--', 'Color', '#eb9e34');
end

ax.XLabel.String = 'Time (s)';
ax.YLabel.String = 'Amplitude (\muV)';

axs(2) = ax;

[b,a] = butter(4, 150/data_cnt.fsample*2, 'high');
ax = subplot(spdim(1), spdim(2), 3);
plot(data_cnt.time{1}, filtfilt(b,a,data_cnt.trial{1}(ch_of_interest,:)));
ax.Title.String = '150 Hz HP';
hold on
for kk = 1:length(hfotrigs_sk_s)
  plot([1 1] * hfotrigs_sk_s(kk), ax.YLim+[0.01 -0.01], 'k--');
end
for kk = 1:length(hfotrigs_z_s)
  plot([1 1] * hfotrigs_z_s(kk), ax.YLim+[0.01 -0.01], 'b--');
end
for kk = 1:length(gs_trigs_ch)
  plot([1 1] * gs_trigs_ch(kk), ax.YLim+[0.01 -0.01], '--', 'Color', '#eb9e34');
end

ax.XLabel.String = 'Time (s)';
ax.YLabel.String = 'Amplitude (\muV)';

axs(3) = ax;


ax = subplot(spdim(1), spdim(2), 4);
imagesc(data_cnt.time{1},[1 400], tfrdat);
axis xy
ax.Title.String = 'Spectrogram';
ax.YLim = [0 400];
ax.CLim = [0 3];
hold on
% for kk = 1:length(hfotrigs)
%   plot([1 1] * data_cnt.time{1}(hfotrigs(kk)), ax.YLim+[0.01 -0.01], 'w:');
% end
% for kk = 1:length(gs_trigs_ch)
%   plot([1 1] * gs_trigs_ch(kk), ax.YLim+[0.01 -0.01], 'y:');
% end

ax.XLabel.String = 'Time (s)';
ax.YLabel.String = 'Frequency (Hz)';

tmp = ax.Position;
cb = colorbar('Location', 'SouthOutside');
ax.Position = tmp;

cb.Label.String = 'Amplitude (z-score)';

axs(4) = ax;



linkaxes(axs, 'x');



%% automatically save screenshots
gensavepath = @(l,x) fullfile('/d/gmi/1/simeon/ieeg_task_analysis/analysis/HFO/discrepancy-2', [dataset, '_', l, sprintf('%03d', x), '.png']);

for kk = 1:length(hfotrigs_sk_s)
  ax.XLim = hfotrigs_sk_s(kk) + [-0.5 0.5];
  saveas(hf, gensavepath('sk', kk));
end
for kk = 1:length(hfotrigs_z_s)
  ax.XLim = hfotrigs_z_s(kk) + [-0.5 0.5];
  saveas(hf, gensavepath('z', kk));
end
