function [hf] = ieeg_plotHFOdiagnostics(time, raw, hfotrigs, tfrdat, debug, gold_standard, varargin)

ip = inputParser;
addParameter(ip, 'showTFRinfo', true);
addParameter(ip, 'showTFRtext', true);
addParameter(ip, 'FigureSize', [1900 1000]);
parse(ip, varargin{:});

hfotrigs = round(hfotrigs);
gold_standard = round(gold_standard);

srate = 1/median(diff(time));

%% plotting

axs = [];

spdim = [4, 1];

hf = figure;
hf.Position(3:4) = ip.Results.FigureSize;
hf.PaperSize = hf.Position(3:4) / 96;
hf.PaperPosition = [0 0 hf.PaperSize];

ax = subplot(spdim(1), spdim(2), 1);
plot(time, raw);
ax.Title.String = '1-500 Hz';
hold on
for kk = 1:length(hfotrigs)
  plot([1 1] * time(hfotrigs(kk)), ax.YLim+[0.01 -0.01], 'k:');
end

scatter(time(gold_standard)', zeros(length(gold_standard),1), '+', 'MarkerEdgeColor', '#eb9e34');


ax.XLabel.String = 'Time (s)';
ax.YLabel.String = 'Amplitude (\muV)';


axs(1) = ax;



[b,a] = butter(4, 80/srate*2, 'high');
ax = subplot(spdim(1), spdim(2), 2);
plot(time, filtfilt(b,a,raw));
ax.Title.String = '80 Hz HP';
hold on
% for kk = 1:length(hfotrigs)
%   plot([1 1] * time(hfotrigs(kk)), ax.YLim+[0.01 -0.01], 'k:');
% end
% for kk = 1:length(gs_trigs_ch)
%   plot([1 1] * gs_trigs_ch(kk), ax.YLim+[0.01 -0.01], '--', 'Color', '#eb9e34');
% end

ax.XLabel.String = 'Time (s)';
ax.YLabel.String = 'Amplitude (\muV)';

axs(2) = ax;

[b,a] = butter(4, 150/srate*2, 'high');
ax = subplot(spdim(1), spdim(2), 3);
plot(time, filtfilt(b,a,raw));
ax.Title.String = '150 Hz HP';
hold on
% for kk = 1:length(hfotrigs)
%   plot([1 1] * time(hfotrigs(kk)), ax.YLim+[0.01 -0.01], 'k:');
% end
% for kk = 1:length(gs_trigs_ch)
%   plot([1 1] * gs_trigs_ch(kk), ax.YLim+[0.01 -0.01], '--', 'Color', '#eb9e34');
% end

ax.XLabel.String = 'Time (s)';
ax.YLabel.String = 'Amplitude (\muV)';

axs(3) = ax;


ax = subplot(spdim(1), spdim(2), 4);
imagesc(debug.spec.time, debug.spec.freq, tfrdat);
axis xy
ax.Title.String = 'Spectrogram';
ax.YLim = [0 400];
ax.CLim = [0 10];
hold on
% for kk = 1:length(hfotrigs)
%   plot([1 1] * time(hfotrigs(kk)), ax.YLim+[0.01 -0.01], 'w:');
% end

if ip.Results.showTFRinfo
  for kk = 1:length(debug.tfrprop)
    if debug.length_insufficient(kk) >= debug.ip.Results.ncycles
      if debug.ishfo(kk)
        lc = 'g';
      else
        lc = 'w';
      end
      pos = debug.tfrprop(kk).BoundingBox;
      pos = [debug.spec.time(round(pos(1))), debug.spec.freq(round(pos(2))), ...
        debug.spec.time(floor(pos(1)+pos(3))), debug.spec.freq(floor(pos(2)+pos(4)))];
      pos(3) = pos(3) - pos(1);
      pos(4) = pos(4) - pos(2);
      
      rectangle('Position', pos, 'EdgeColor', lc);
      
      if ip.Results.showTFRtext
        % make sure textpos doesn't exceed edges
        textpos = [pos(1)+pos(3), pos(2)+pos(4)];
        textpos(2) = min(textpos(2), max(debug.spec.freq) * 0.9);
        textpos(2) = max(textpos(2), max(debug.spec.freq) * 0.1);
        
        hascircbw = isfield(debug, 'circularity');
        if hascircbw
          strcircbw = '\nCirc %.1f   BW %.1f';
          sprintfvals = {debug.circularity(kk), ...
            debug.hfo_freq_width(kk)};
        else
          strcircbw = '';
          sprintfvals = {};
        end
        
        text(textpos(1), textpos(2), ...
          sprintf(['idx: %d\nCentroid: %.1f, %.1f\nPos %.1f   Bounds %.1f   Cycles %.1f', strcircbw], ...
          kk, ...
          debug.tfrprop(kk).Centroid(1), debug.tfrprop(kk).Centroid(2), ...
          debug.centroid_position(kk), ...
          debug.bounds_in_lf(kk), ...
          debug.length_insufficient(kk), ...
          sprintfvals{:}), 'Color', 'w', 'FontSize', 8);
      end
    end
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

cb.Label.String = 'Amplitude';

axs(4) = ax;


linkaxes(axs, 'x');
