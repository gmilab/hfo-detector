
addpath /d/gmi/1/toolboxes/fieldtrip/git-latest/
ft_defaults
addpath /d/gmi/1/simeon/hsc-eeg-tools/MATLAB/

load('/d/gmi/1/simeon/hfo_detector/donos_simulated/simulation_data2.mat')

nsamp = size(databip.channel,1);
nchan = size(databip.channel,2);


param = 'thresh';
param_types = 2:0.25:10;

hfofun = @ieeg_hfodetector_boundingbox;
% hfofun = @ieeg_hfodetector_bb_orig;


nps = length(param_types);
confusion_data = table();
confusion_data.(param) = zeros(nps,1);
confusion_data.TP = zeros(nps,1);
confusion_data.FP = zeros(nps,1);
confusion_data.FN = zeros(nps,1);
confusion_data.PPV = zeros(nps,1);
confusion_data.SEN = zeros(nps,1);
for pppp = 1:nps
  
  data_cnt = [];
  data_cnt.trial = {databip.channel'};
  data_cnt.fsample = databip.fs;
  data_cnt.time = {0:(1/databip.fs):((nsamp-1)/databip.fs)};
  data_cnt.label = databip.label;
  
  
  hfotrigs = hfofun(data_cnt, 1, param, param_types(pppp));
  
  gold_standard = 1000:2000:600000;
  
  confusion_data.(param)(pppp) = param_types(pppp);
  
  
  %% generate confusion matrix  
  for nn = 1:length(hfotrigs)
    if any(abs(gold_standard - hfotrigs(nn)) < 500)
      confusion_data.TP(pppp) = confusion_data.TP(pppp) + 1;
    else
      confusion_data.FP(pppp) = confusion_data.FP(pppp) + 1;
    end
  end
  
  for nn = 1:length(gold_standard)
    if ~any(abs(gold_standard(nn) - hfotrigs) < 500)
      confusion_data.FN(pppp) = confusion_data.FN(pppp) + 1;
    end
  end
  
  sum_TP = confusion_data.TP(pppp);
  sum_FP = confusion_data.FP(pppp);
  sum_FN = confusion_data.FN(pppp);
  
  confusion_data.PPV(pppp) = sum_TP/(sum_TP+sum_FP);
  confusion_data.SEN(pppp) = sum_TP/(sum_TP+sum_FN);
  
  fprintf('\n\n\nTP: %d\tFP: %d\tFN:%d\n', sum_TP, sum_FP, sum_FN);
  fprintf('PPV: %.5f\tSensitivity: %.5f\n\n\n', confusion_data.PPV(pppp), confusion_data.SEN(pppp))
  
end

writetable(confusion_data, sprintf('/d/gmi/1/simeon/ieeg_task_analysis/analysis/HFO/simulations/donos_%s_%s.xlsx', func2str(hfofun), param));



%%
hf = figure;
hf.Position = [968 857 560 451];
hf.PaperPosition = [0 0 hf.Position(3:4)/96];
hf.PaperSize = hf.PaperPosition(3:4);

hp = plot(1-confusion_data.SEN, confusion_data.PPV, 'LineWidth', 2);

hp.DataTipTemplate.DataTipRows(1).Label = 'FNR';
hp.DataTipTemplate.DataTipRows(2).Label = 'PPV'; 

row = dataTipTextRow(strrep(param, '_', '\_'), confusion_data.(param));
hp.DataTipTemplate.DataTipRows(end+1) = row;

ax = gca;
ax.Title.String = 'Detector performance at different relative power thresholds';
ax.XLim = [0 0.1];
ax.YLim = [0.9 1];

ax.YLabel.String = 'Positive Predictive Value';
ax.XLabel.String = 'False Negative Rate';

print(hf, sprintf('/d/gmi/1/simeon/ieeg_task_analysis/analysis/HFO/simulations/ROC_donos_%s_%s.png', func2str(hfofun), param), '-dpng', '-r600');
saveas(hf, sprintf('/d/gmi/1/simeon/ieeg_task_analysis/analysis/HFO/simulations/ROC_donos_%s_%s.pdf', func2str(hfofun), param));
