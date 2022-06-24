
addpath /d/gmi/1/toolboxes/fieldtrip/git-latest/
ft_defaults
addpath /d/gmi/1/simeon/hsc-eeg-tools/MATLAB/

data = load_anywave('/d/gmi/1/simeon/hfo_detector/anywave/simulation_rate_3.ades');
data.markers.HFO = contains(data.markers.Label, 'Ripple');

nsamp = data.nsamp;
nchan = data.nchan;

data_cnt = [];
data_cnt.trial = {data.data};
data_cnt.fsample = data.srate;
data_cnt.time = {0:(1/data.srate):((nsamp-1)/data.srate)};
data_cnt.label = data.label;




%% run the detector
hfotrigs = cell(nchan,1);
tfrdat = cell(nchan,1);
hfodebug = cell(nchan,1);

for kk = 1:nchan
  hfotrigs{kk} = ieeg_hfodetector_boundingbox(data_cnt, kk, 'circularity', 0.5); % , 'hfo_max_freq_width', 150);
end




%% generate confusion matrix
confusion_data = table();
confusion_data.label = data.label;

confusion_data.TP = zeros(nchan,1);
confusion_data.FP = zeros(nchan,1);
confusion_data.FN = zeros(nchan,1);
confusion_data.PPV = zeros(nchan,1);
confusion_data.SEN = zeros(nchan,1);

for kk = 1:nchan  
  gold_standard = data.markers.Time(data.markers.HFO & contains(data.markers.Misc, [data.label{kk},',']));
  gold_standard = gold_standard .* data.srate; % convert to samples
  
  thishfotrigs = hfotrigs{kk};
  [~,ia] = unique(round(thishfotrigs/500));
  thishfotrigs = thishfotrigs(ia);
  
  for nn = 1:length(thishfotrigs)
    if any(abs(gold_standard - thishfotrigs(nn)) <= 500)
      confusion_data.TP(kk) = confusion_data.TP(kk) + 1;
    else
      confusion_data.FP(kk) = confusion_data.FP(kk) + 1;
    end
  end
  
  for nn = 1:length(gold_standard)
    if ~any(abs(gold_standard(nn) - thishfotrigs) <= 500)
      confusion_data.FN(kk) = confusion_data.FN(kk) + 1;
    end
  end
  
  sum_TP = confusion_data.TP(kk);
  sum_FP = confusion_data.FP(kk);
  sum_FN = confusion_data.FN(kk);
  
  confusion_data.PPV(kk) = sum_TP/(sum_TP+sum_FP);
  confusion_data.SEN(kk) = sum_TP/(sum_TP+sum_FN);
  
%   fprintf('\n\n\nTP: %d\tFP: %d\tFN:%d\n', sum_TP, sum_FP, sum_FN);
%   fprintf('PPV: %.5f\tSensitivity: %.5f\n\n\n', confusion_data.PPV(kk), confusion_data.SEN(kk))
  
end

%% summarize
confusion_summary = table();
confusion_summary.label = {'0dB'; '5dB'; '10dB'; '15dB'};

confusion_summary.TP = zeros(height(confusion_summary),1);
confusion_summary.FP = zeros(height(confusion_summary),1);
confusion_summary.FN = zeros(height(confusion_summary),1);
confusion_summary.PPV = zeros(height(confusion_summary),1);
confusion_summary.SEN = zeros(height(confusion_summary),1);

for kk = 1:height(confusion_summary)
  rows_idx = contains(confusion_data.label, confusion_summary.label{kk});
  
  confusion_summary.TP(kk) = sum(confusion_data.TP(rows_idx), 1);
  confusion_summary.FP(kk) = sum(confusion_data.FP(rows_idx), 1);
  confusion_summary.FN(kk) = sum(confusion_data.FN(rows_idx), 1);
  
  confusion_summary.PPV(kk) = confusion_summary.TP(kk)/(confusion_summary.TP(kk)+confusion_summary.FP(kk));
  confusion_summary.SEN(kk) = confusion_summary.TP(kk)/(confusion_summary.TP(kk)+confusion_summary.FN(kk));
end

disp(confusion_summary);

return

%% plot spectrogram with hfo markers
kk = 178;

[hfotrigs{kk},tfrdat{kk},~,hfodebug{kk}] = ieeg_hfodetector_boundingbox(data_cnt, kk);

gold_standard = data.markers.Time(data.markers.HFO & contains(data.markers.Misc, [data.label{kk},',']));
gold_standard = gold_standard .* data.srate; % convert to samples

hf = ieeg_plotHFOdiagnostics(data_cnt.time{1}, data_cnt.trial{1}(kk,:), ...
  hfotrigs{kk}, tfrdat{kk}, hfodebug{kk}, gold_standard, ...
  'showTFRinfo', true);

%%
hf.Position(3:4) = [400 1000];
hf.PaperSize = hf.Position(3:4) / 96;
hf.PaperPosition = [0 0 hf.PaperSize];

[ui_file, ui_path] = uiputfile('/d/gmi/1/simeon/ieeg_task_analysis/analysis/HFO/detector_explainability_gui.png');
print(hf, fullfile(ui_path, ui_file), '-r600', '-dpng');




