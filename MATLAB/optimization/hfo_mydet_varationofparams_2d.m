
addpath /d/gmi/1/toolboxes/fieldtrip/git-latest/
ft_defaults
addpath /d/gmi/1/simeon/hsc-eeg-tools/MATLAB/

load('/d/gmi/1/simeon/hfo_detector/donos_simulated/simulation_data2.mat')

nsamp = size(databip.channel,1);
nchan = size(databip.channel,2);


param = {'ncycles', 'wavelet_width'};
param_types = {2:0.25:10, 3:0.25:13};


nps = prod(cellfun(@length, param_types));
confusion_data = table();

[ix,iy] = meshgrid(param_types{1}, param_types{2});
confusion_data.(param{1}) = reshape(ix, nps, 1);
confusion_data.(param{2}) = reshape(iy, nps, 1);

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
  
  
  hfotrigs = ieeg_hfodetector_boundingbox(data_cnt, 1,...
    param{1}, confusion_data.(param{1})(pppp), ...
    param{2}, confusion_data.(param{2})(pppp));
  
  gold_standard = 1000:2000:600000;
  
  
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



return


%%
hf = figure;
hf.Position = [968 857 560 451];
hold on

clr = lines(50);

for pmp1 = 1:length(param_types{1})
  % get rows matching this param type
  rowidx = confusion_data.(param{1}) == param_types{1}(pmp1);
  
  if param_types{1}(pmp1) == 4.75
    plotparams = {'Color', clr(1,:), 'LineWidth', 3};
  else
    plotparams = {'Color', [clr(pmp1,:) 0.3], 'LineWidth', 1};
  end
  
  hp = plot(1-confusion_data.SEN(rowidx), confusion_data.PPV(rowidx), plotparams{:});
  
  hp.DataTipTemplate.DataTipRows(1).Label = 'FNR';
  hp.DataTipTemplate.DataTipRows(2).Label = 'PPV';
  
  row = dataTipTextRow(strrep(param{1}, '_', '\_'), confusion_data.(param{1})(rowidx));
  hp.DataTipTemplate.DataTipRows(end+1) = row;
  
  row = dataTipTextRow(strrep(param{2}, '_', '\_'), confusion_data.(param{2})(rowidx));
  hp.DataTipTemplate.DataTipRows(end+1) = row;
end

ax = gca;
ax.XLim = [0 0.1];
ax.YLim = [0.9 1];


% title(['ROC - ' strjoin(param, '/')], 'Interpreter', 'none')
ylabel('Positive Predictive Value')
xlabel('False Negative Rate')
