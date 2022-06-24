
addpath /d/gmi/1/toolboxes/fieldtrip/git-latest/
ft_defaults
addpath /d/gmi/1/simeon/hsc-eeg-tools/MATLAB/
addpath /d/gmi/1/simeon/ieeg_task_analysis/code/kale/functions/

load('/d/gmi/1/simeon/hfo_detector/donos_simulated/simulation_data1.mat')

nsamp = size(databip.channel,1);
nchan = size(databip.channel,2);

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
  [hfotrigs{kk},tfrdat{kk},hfodet_ip{kk},debug{kk}] = ieeg_hfodetector_boundingbox(data_cnt, kk);
end


gold_standard = 1000:2000:600000;


kk = 1;

% get value scores at trigger (max in column, avg +/- 0.05 s)
hf_idx = fcp_closest_idx(debug{kk}.spec.freq, hfodet_ip{kk}.Results.foi);
hits = table();
hits.max = zeros(length(gold_standard),1);
hits.mean = zeros(length(gold_standard),1);
hits.std = zeros(length(gold_standard),1);
for tt = 1:length(gold_standard)
  [~,colmax] = max(tfrdat{kk}(hf_idx(1):hf_idx(2),gold_standard(tt)));
  
  dat = tfrdat{kk}((colmax+hf_idx(1)-1)+(-1:1),gold_standard(tt) + (-100:100));
  hits.max(tt) = max(dat(:));
  hits.mean(tt) = mean(dat(:));
  hits.std(tt) = std(dat(:));
end


% get other highest
tfrexclusion = tfrdat{kk}(hf_idx(1):hf_idx(2),:);
for tt = 1:length(gold_standard)
  tfrexclusion(:,gold_standard(tt) + [-100 100]) = 0;
end

highvals = prctile(tfrexclusion(:), [75 85 95 99]);


return


%% plot spectrogram with hfo markers
ch_of_interest = 5;
hfotrigs_1 = round(hfotrigs{ch_of_interest});

axs = [];

spdim = [4, 1];

hf = figure; 
hf.Position(3:4) = [493        1311];

ax = subplot(spdim(1), spdim(2), 1);
plot(data_cnt.time{1}, data_cnt.trial{1}(ch_of_interest,:));
ax.Title.String = '1-500 Hz';
hold on
for kk = 1:length(hfotrigs_1)
  plot([1 1] * data_cnt.time{1}(hfotrigs_1(kk)), ax.YLim+[0.01 -0.01], 'k:');
end

scatter(data_cnt.time{1}(round(hfotrigs{5}))', 2*ones(length(hfotrigs{5}),1), '+', 'MarkerEdgeColor', '#eb9e34');


ax.XLabel.String = 'Time (s)';
ax.YLabel.String = 'Amplitude (\muV)';


axs(1) = ax;



[b,a] = butter(4, 80/data_cnt.fsample*2, 'high');
ax = subplot(spdim(1), spdim(2), 2);
plot(data_cnt.time{1}, filtfilt(b,a,data_cnt.trial{1}(ch_of_interest,:)));
ax.Title.String = '80 Hz HP';
hold on
% for kk = 1:length(hfotrigs_1)
%   plot([1 1] * data_cnt.time{1}(hfotrigs_1(kk)), ax.YLim+[0.01 -0.01], 'k:');
% end
% for kk = 1:length(gs_trigs_ch)
%   plot([1 1] * gs_trigs_ch(kk), ax.YLim+[0.01 -0.01], '--', 'Color', '#eb9e34');
% end

ax.XLabel.String = 'Time (s)';
ax.YLabel.String = 'Amplitude (\muV)';

axs(2) = ax;

[b,a] = butter(4, 150/data_cnt.fsample*2, 'high');
ax = subplot(spdim(1), spdim(2), 3);
plot(data_cnt.time{1}, filtfilt(b,a,data_cnt.trial{1}(ch_of_interest,:)));
ax.Title.String = '150 Hz HP';
hold on
% for kk = 1:length(hfotrigs_1)
%   plot([1 1] * data_cnt.time{1}(hfotrigs_1(kk)), ax.YLim+[0.01 -0.01], 'k:');
% end
% for kk = 1:length(gs_trigs_ch)
%   plot([1 1] * gs_trigs_ch(kk), ax.YLim+[0.01 -0.01], '--', 'Color', '#eb9e34');
% end

ax.XLabel.String = 'Time (s)';
ax.YLabel.String = 'Amplitude (\muV)';

axs(3) = ax;


ax = subplot(spdim(1), spdim(2), 4);
imagesc(debug{ch_of_interest}.spec.time, debug{ch_of_interest}.spec.freq, tfrdat{ch_of_interest});
axis xy
ax.Title.String = 'Spectrogram';
ax.YLim = [0 400];
ax.CLim = [0 10];
hold on
% for kk = 1:length(hfotrigs_1)
%   plot([1 1] * data_cnt.time{1}(hfotrigs_1(kk)), ax.YLim+[0.01 -0.01], 'w:');
% end
for kk = 1:length(debug{ch_of_interest}.tfrprop)
  if ~debug{ch_of_interest}.length_insufficient(kk)
    if debug{ch_of_interest}.ishfo(kk)
      lc = 'g';
    else
      lc = 'w';
    end
    pos = debug{ch_of_interest}.tfrprop(kk).BoundingBox;
    pos = [debug{ch_of_interest}.spec.time(round(pos(1))), debug{ch_of_interest}.spec.freq(round(pos(2))), ...
      debug{ch_of_interest}.spec.time(floor(pos(1)+pos(3))), debug{ch_of_interest}.spec.freq(floor(pos(2)+pos(4)))];
    pos(3) = pos(3) - pos(1);
    pos(4) = pos(4) - pos(2);
    
    rectangle('Position', pos, 'EdgeColor', lc);
    text(pos(1)+pos(3), pos(2)+pos(4), ...
      sprintf('idx: %d\nCentroid: %.1f, %.1f\nPos: %d\nBounds: %d\nCycles: %d', ...
      kk, ...
      debug{ch_of_interest}.tfrprop(kk).Centroid(1), debug{ch_of_interest}.tfrprop(kk).Centroid(2), ...
      debug{ch_of_interest}.centroid_position(kk), ...
      debug{ch_of_interest}.bounds_in_lf(kk), ...
      debug{ch_of_interest}.length_insufficient(kk)), 'Color', 'w', 'FontSize', 8);
  end
end
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



