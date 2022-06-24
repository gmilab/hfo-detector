
addpath /d/gmi/1/toolboxes/fieldtrip/git-latest/
ft_defaults
addpath /d/gmi/1/simeon/hsc-eeg-tools/MATLAB/
addpath /d/gmi/1/simeon/hfo_detector/noise

fs = 2000;
dur = 300;
noise = pinknoise(fs*(dur+10), 1);

[btrb,btra] = butter(4, [0.5 500]/(fs/2));
noise = filtfilt(btrb, btra, noise);
noise = noise(fs*5:end-(fs*5));

nsamp = size(noise,1);
nchan = size(noise,2);

data_cnt = [];
data_cnt.trial = {noise'};
data_cnt.fsample = fs;
data_cnt.time = {0:(1/fs):((nsamp-1)/fs)};
data_cnt.label = {'Signal'};

hfotrigs = cell(nchan,1);
tfrdat = cell(nchan,1);
hfodet_ip = cell(nchan,1);
debug = cell(nchan,1);
for kk = 1:nchan
  [hfotrigs{kk},tfrdat{kk},hfodet_ip{kk},debug{kk}] = ieeg_hfodetector_boundingbox(data_cnt, kk);
end






%% plot spectrogram with hfo markers
ch_of_interest = 1;
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
% for kk = 1:length(gs_trigs_ch)
%   plot([1 1] * gs_trigs_ch(kk), ax.YLim+[0.01 -0.01], '--', 'Color', '#eb9e34');
% end


ax.XLabel.String = 'Time (s)';
ax.YLabel.String = 'Amplitude (\muV)';


axs(1) = ax;



[b,a] = butter(4, 80/data_cnt.fsample*2, 'high');
ax = subplot(spdim(1), spdim(2), 2);
plot(data_cnt.time{1}, filtfilt(b,a,data_cnt.trial{1}(ch_of_interest,:)));
ax.Title.String = '80 Hz HP';
hold on
for kk = 1:length(hfotrigs_1)
  plot([1 1] * data_cnt.time{1}(hfotrigs_1(kk)), ax.YLim+[0.01 -0.01], 'k:');
end
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
for kk = 1:length(hfotrigs_1)
  plot([1 1] * data_cnt.time{1}(hfotrigs_1(kk)), ax.YLim+[0.01 -0.01], 'k:');
end
% for kk = 1:length(gs_trigs_ch)
%   plot([1 1] * gs_trigs_ch(kk), ax.YLim+[0.01 -0.01], '--', 'Color', '#eb9e34');
% end

ax.XLabel.String = 'Time (s)';
ax.YLabel.String = 'Amplitude (\muV)';

axs(3) = ax;


ax = subplot(spdim(1), spdim(2), 4);
imagesc(data_cnt.time{1},[2 400], tfrdat{ch_of_interest});
axis xy
ax.Title.String = 'Spectrogram';
ax.YLim = [0 400];
ax.CLim = [0 15];
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
