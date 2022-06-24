
clear
addpath /d/gmi/1/toolboxes/fieldtrip/git-latest/
ft_defaults
addpath /d/gmi/1/simeon/hsc-eeg-tools/MATLAB/

data = load_anywave('/d/gmi/1/simeon/hfo_detector/anywave/simulation_rate_3.ades');
data.markers.HFO = contains(data.markers.Label, 'Ripple');

labels = {'0dB'; '5dB'; '10dB'; '15dB'};

rowsubsetidx = find(contains(data.label, labels{2}));

nsamp = data.nsamp;
nchan = data.nchan;

data_cnt = [];
data_cnt.trial = {data.data};
data_cnt.fsample = data.srate;
data_cnt.time = {0:(1/data.srate):((nsamp-1)/data.srate)};
data_cnt.label = data.label;

param = 'circularity';
param_types = 0.2:0.1:0.9;

% param = 'hfo_max_freq_width';
% param_types = 80:10:300;



nps = length(param_types);
confusion_data = table();
confusion_data.(param) = reshape(param_types, [], 1);
confusion_data.TP = zeros(nps,1);
confusion_data.FP = zeros(nps,1);
confusion_data.FN = zeros(nps,1);
confusion_data.PPV = zeros(nps,1);
confusion_data.SEN = zeros(nps,1);
for pppp = 1:nps
  cthresh_confusion = table();
  cthresh_confusion.label = data.label(rowsubsetidx);
  
  cthresh_confusion.TP = zeros(length(rowsubsetidx),1);
  cthresh_confusion.FP = zeros(length(rowsubsetidx),1);
  cthresh_confusion.FN = zeros(length(rowsubsetidx),1);
  
  for kk = 1:length(rowsubsetidx)
    thishfotrigs = ieeg_hfodetector_boundingbox(data_cnt, rowsubsetidx(kk), param, param_types(pppp));
    
    gold_standard = data.markers.Time(data.markers.HFO & contains(data.markers.Misc, [data.label{kk},',']));
    gold_standard = gold_standard .* data.srate; % convert to samples
    
    [~,ia] = unique(round(thishfotrigs/1000));
    thishfotrigs = thishfotrigs(ia);
    
    for nn = 1:length(thishfotrigs)
      if any(abs(gold_standard - thishfotrigs(nn)) <= 500)
        cthresh_confusion.TP(kk) = cthresh_confusion.TP(kk) + 1;
      else
        cthresh_confusion.FP(kk) = cthresh_confusion.FP(kk) + 1;
      end
    end
    
    for nn = 1:length(gold_standard)
      if ~any(abs(gold_standard(nn) - thishfotrigs) <= 500)
        cthresh_confusion.FN(kk) = cthresh_confusion.FN(kk) + 1;
      end
    end
  end
  
  confusion_data.TP(pppp) = sum(cthresh_confusion.TP);
  confusion_data.FP(pppp) = sum(cthresh_confusion.FP);
  confusion_data.FN(pppp) = sum(cthresh_confusion.FN);
  
  confusion_data.PPV(pppp) = confusion_data.TP(pppp)/(confusion_data.TP(pppp)+confusion_data.FP(pppp));
  confusion_data.SEN(pppp) = confusion_data.TP(pppp)/(confusion_data.TP(pppp)+confusion_data.FN(pppp));
end

writetable(confusion_data, sprintf('/d/gmi/1/simeon/ieeg_task_analysis/analysis/HFO/simulations/anywave_%s.xlsx', param));

return





