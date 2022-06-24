

clear

fcp_setup;
addpath /d/gmi/1/toolboxes/circstat/
addpath /d/gmi/1/simeon/hsc-eeg-tools/MATLAB/
addpath /d/gmi/1/simeon/hfo_detector/easyh5/
savefolder = '/d/gmi/1/simeon/ieeg_task_analysis/analysis/GroupFinal/HFOsimulations';


errorlog = {};

% dataset = '/d/gmi/1/simeon/hfo_detector/hfo_data_annotated/Data_Patient_01.h5';

datasetpath = '/d/gmi/1/simeon/hfo_detector/hfo_data_annotated/Intraoperative_ECoG_HFO/data';

datasets = dir(datasetpath);
datasets = [{datasets.name};];
datasets = datasets(~startsWith(datasets, '.'));


% varparams = 2:0.25:6;
varparams = 3;

for idx_varparam = 1:length(varparams)
  prm = varparams(idx_varparam);
  
  confusion_data = table();
  confusion_data.filename = cell(0,1);
  confusion_data.channel = zeros(0,1);
  confusion_data.TP = zeros(0,1);
  confusion_data.FP = zeros(0,1);
  confusion_data.FN = zeros(0,1);
  
  cd_idx = 0;
  
  
  for ddd = 1:2 % length(datasets)
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
      
      % remove zeros at the beginning or end (empty data)
      allzeros = all(data == 0, 2);
      firstidx = find(~allzeros, 1, 'first');
      lastidx = find(~allzeros, 1, 'last');
      
      data = data(firstidx:lastidx, :);

      
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
      
      
      
      for ch_of_interest = 1:nch
        try
          gs_trigs_ch = trigs(trigs(:,1) == ch_of_interest, 2);
          
          %% Calculate trial hfo dichotomization by previously-saved TFR method
%           [samples_with_hfo,lfactive,hfactive,hdetect_tfrdat,hfodet_ip] = ieeg_hfodetector(data_cnt, ch_of_interest, 'minlen', prm);
%           
%           hfotrigs = bwconncomp(samples_with_hfo);
%           hfotrigs = reshape(cellfun(@(x) round(median(x)), hfotrigs.PixelIdxList), [], 1);

          [hfotrigs,tfrdat,hfodet_ip,hfodebug] = ieeg_hfodetector_boundingbox(data_cnt, ch_of_interest, 'hfo_freq_width', 300);
          
          hfotrigs_s = ((hfotrigs-1)*hfodet_ip.Results.downsample+1)/fs;
          
          fprintf('%d HFOs detected / %d HFOs ground truth\n', length(hfotrigs), sum(trigs(:,1) == ch_of_interest));
          
          
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
          
          %% Save plotting for diagnostics
          hf = ieeg_plotHFOdiagnostics(data_cnt.time{1}, data_cnt.trial{1}(ch_of_interest,:), ...
            hfotrigs, tfrdat, hfodebug, gs_trigs_ch * fs, ...
            'showTFRinfo', false);
          print(hf, sprintf('/d/gmi/1/simeon/ieeg_task_analysis/analysis/HFO/zurich_test/%d-%s-%d.png', cd_idx, dataset, ch_of_interest), '-r450', '-dpng');
          close(hf);
          
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
  
%   writetable(confusion_data, sprintf('/d/gmi/1/simeon/ieeg_task_analysis/analysis/HFO-ROC/us/newdet-zhf%d.csv', prm))
  
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
ch_of_interest = 8;

gs_trigs_ch = trigs(trigs(:,1) == ch_of_interest, 2);
[hfotrigs,tfrdat,hfodet_ip,hfodebug] = ieeg_hfodetector_boundingbox(data_cnt, ch_of_interest, 'hfo_freq_width', 300);
hfotrigs_s = ((hfotrigs-1)*hfodet_ip.Results.downsample+1)/fs;

hf = ieeg_plotHFOdiagnostics(data_cnt.time{1}, data_cnt.trial{1}(ch_of_interest,:), ...
  hfotrigs, tfrdat, hfodebug, gs_trigs_ch * fs, ...
  'showTFRinfo', true);


