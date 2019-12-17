
function ftc = synftc(stimspets, spetfreq, analwin, channels, trials, chanlist, frqs, ax)
%synftc creates an image plot of channel (x-axis) x frequency (y-axis) x
%spike rate (color). It is called in synui during livestream and offline
%analysis.

%stimspets = spike times aligned to last stimulus (nspikes x 1, seconds)
%spetfreq = frequency of last stimulus, for each spike (nspikes x 1, seconds)
%analwin = analysis window, relative to last stimulus (2 x 1, seconds)
%channels = channel of each spike (nspikes x 1)
%trials = number of trials (1x1)
%chanlist = selected channels for plotting (nchans x 1)
%frqs = frequency (or level) of each stimulus (nstims x 1)
%ax = axis for plotting (axis handle)
%threshcrit = spontaneous rate threshold rate for each channel (nchans x 1, Hz)

%ftc = spike rates of each channel x frequency (nchan x nfreq, Hz)
%ftcstats - parameters for gaussian fit of fra for each channel. column 1:
%amplitude (Hz), column 2: best frequency (kHz), column 3: bandwidth (2 stds, octaves),
%column 4: lower cutoff frequency (kHz), column 5: upper cutoff frequency
%(kHz)
%(nchans x %4)
if ~exist('ax', 'var')
    ax = axes();
end
nchans = numel(chanlist);
frqlist = sort(unique(frqs));
ftc = NaN(nchans, numel(frqlist));
scale = 1./(analwin(2)-analwin(1))./trials;
for chan = 1:nchans
    ftc(chan,:) = hist(spetfreq(stimspets>analwin(1)& stimspets<analwin(2) & ...
        channels ==chanlist(chan)), frqlist);
end
ftc = ftc .*scale;
imagesc(ax, ftc);
colormap(ax,'jet')

% if exist('threshcrit', 'var')
%     for chan = 1:nchans
%          try
%              if numel(frqlist) < 10
%                 a = fit(frqlist, ftc(chan, :)'-threshcrit(chan), 'gauss1');
%                 if a.a1 > 10 %threshold of 10 Hz above spontaneous works okay.
%                  frqthresh = [-a.c1 * sqrt(log(a.a1/10)) + a.b1, a.c1 * sqrt(log(a.a1/10)) + a.b1] ;
%                  ftcstats(chan, :) = [a.a1, a.b1, 2*a.c1, frqthresh];
% %             figure(2); clf(2); hold on;
% %             rectangle('Position', [frqthresh(1), 0, frqthresh(2)-frqthresh(1), a.a1 + threshcrit(chan)], 'FaceColor', [0.8 0.8 1], 'EdgeColor', [1 1 1]);
% %             plot(frqlist, ftc(chan, :)', 'k.');
% %             plot(frqlist, a.a1.*exp(-((frqlist-a.b1)./a.c1).^2)+threshcrit(chan), 'r:') 
% %             [threshcrit(chan) frqthresh(1) a.b1]
% %             pause
% %             hold off;
%                 end
%              else
%                 a = fit(log(frqlist), ftc(chan, :)'-threshcrit(chan), 'gauss1');
% %                 frqthresh = [exp(-a.c1 * sqrt(log(a.a1)) + a.b1) exp(a.c1 * sqrt(log(a.a1)) + a.b1) ];
%                 ftcstats(chan, :) = [a.a1, exp(a.b1), 2*a.c1/log(2), exp([a.b1 - a.c1, a.b1 + a.c1])];
% %             figure(2); clf(2); hold on;
% %             rectangle('Position', [frqthresh(1), 0, frqthresh(2)-frqthresh(1), a.a1 + threshcrit(chan)], 'FaceColor', [0.8 0.8 1], 'EdgeColor', [1 1 1]);
% %             plot(frqlist, ftc(chan, :)', 'k.');
% %             plot(frqlist, a.a1.*exp(-((log(frqlist)-a.b1)./a.c1).^2)+threshcrit(chan), 'r:') 
% %             hold off;
%              end
%         catch
%         end
%     end
%     
% end

set(ax,'XTick', 1:numel(frqlist))
set(ax,'YTick', 1:nchans)
set(ax,'YTickLabel', chanlist)
ylabel(ax, 'Channel')

colorbar(ax)
axis(ax, [.5, numel(frqlist)+.5, .5, nchans+.5 ])
hosatchel = get(ax, 'colormap');
hosatchel(1,:) = [1 1 1];
set(ax, 'colormap', hosatchel);