
addpath /d/gmi/1/toolboxes/fieldtrip/git-latest/
ft_defaults
addpath /d/gmi/1/simeon/hsc-eeg-tools/MATLAB/

load('/d/gmi/1/simeon/hfo_detector/donos_simulated/simulation_data1.mat')

nsamp = size(databip.channel,1);
nchan = size(databip.channel,2);
% nchan = 1;

data_cnt = [];
data_cnt.trial = {databip.channel'};
data_cnt.fsample = databip.fs;
data_cnt.time = {0:(1/databip.fs):((nsamp-1)/databip.fs)};
data_cnt.label = databip.label;

hfotrigs = cell(nchan,1);
tfrdat = cell(nchan,1);
hfodet_ip = cell(nchan,1);
debug = cell(nchan,1);
for kk = 1:nchan
  [hfotrigs{kk},tfrdat{kk},debug{kk}] = ieeg_hfodetector_boundingbox(data_cnt, kk);
end

gold_standard = 1000:2000:600000;


%% generate confusion matrix
ch_of_interest = 1;

confusion_data = struct('TP', 0, 'FP', 0, 'FN', 0);

  thishfotrigs = hfotrigs{ch_of_interest};
  [~,ia] = unique(round(thishfotrigs/1000));
  thishfotrigs = thishfotrigs(ia);

for nn = 1:length(thishfotrigs)
  if any(abs(gold_standard - thishfotrigs(nn)) < 500)
    confusion_data.TP = confusion_data.TP + 1;
  else
    confusion_data.FP = confusion_data.FP + 1;
  end
end

for nn = 1:length(gold_standard)
  if ~any(abs(gold_standard(nn) - thishfotrigs) < 500)
    confusion_data.FN = confusion_data.FN + 1;
  end
end

sum_TP = sum(confusion_data.TP);
sum_FP = sum(confusion_data.FP);
sum_FN = sum(confusion_data.FN);
fprintf('\n\n\nTP: %d\tFP: %d\tFN:%d\n', sum_TP, sum_FP, sum_FN);
fprintf('PPV: %.5f\tSensitivity: %.5f\n\n\n', sum_TP/(sum_TP+sum_FP), sum_TP/(sum_TP+sum_FN))

hf = ieeg_plotHFOdiagnostics(data_cnt.time{1}, data_cnt.trial{1}(ch_of_interest,:), ...
  thishfotrigs, tfrdat{ch_of_interest}, debug{ch_of_interest}, gold_standard, ...
  'showTFRinfo', true, 'showTFRtext', false);


return

%%
hf.Position(3:4) = [400 1000];
hf.PaperSize = hf.Position(3:4) / 96;
hf.PaperPosition = [0 0 hf.PaperSize];

[ui_file, ui_path] = uiputfile('/d/gmi/1/simeon/ieeg_task_analysis/analysis/HFO/detector_explainability_gui.png');
print(hf, fullfile(ui_path, ui_file), '-r600', '-dpng');

