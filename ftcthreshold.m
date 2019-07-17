function ftcstats = ftcthreshold(ftc, frqlist, spontrate, threshcrit)
%find "threshold" from frequency tuning curves or rate level functions
%this approach uses least squares optimization to fit to a gaussian the FTC/RLF
%It takes into account spontaneous firing
%Thresholds are defined as the two standard deviations from the gaussian
%mean

nchans = size(ftc, 1);
ftcstats = NaN(nchans, 5);
for chan = 1:nchans
     try
        fo = fitoptions('Method', 'NonlinearLeastSquares', 'Lower', [0 0 0], 'Upper', [1E4 1E4 1E4]);
         if numel(frqlist) < threshcrit %RLN
            a = fit(frqlist, ftc(chan, :)'-spontrate(chan), 'gauss1', fo);
            if a.a1 > threshcrit %threshold of 10 Hz above spontaneous works okay.
             frqthresh = [-a.c1 * sqrt(log(a.a1/10)) + a.b1, a.c1 * sqrt(log(a.a1/10)) + a.b1] ;
             ftcstats(chan, :) = [a.a1, a.b1, 2*a.c1, frqthresh];
            end
            figure(10);
            nrows = ceil(sqrt(nchans));
            ncols = ceil(nchans./nrows);
            subplot(nrows, ncols, chan)
            plot(frqlist, ftc(chan, :)', 'bo');
            hold on;
            plot(frqlist, a.a1.*exp(-((frqlist-a.b1)./a.c1).^2)+spontrate(chan), 'g-', 'LineWidth', 2)
            plot([frqthresh(1) frqthresh(1)], [0 max(ftc(chan,:))], 'r--', 'Linewidth', 4)
            hold off;
            title([num2str(a.a1), ', ' num2str(a.b1), ', ' num2str(a.c1)])
         else %Iso-I
            a = fit(log(frqlist), ftc(chan, :)'-spontrate(chan), 'gauss1');
            ftcstats(chan, :) = [a.a1, exp(a.b1)./1000, 2*a.c1/log(2), exp([a.b1 - a.c1, a.b1 + a.c1])./1000];
            
            figure(10);
            nrows = ceil(sqrt(nchans));
            ncols = ceil(nchans./nrows);
            subplot(nrows, ncols, chan)
            plot(frqlist./1000, ftc(chan, :)', 'bo');
            hold on;
            plot(frqlist./1000, a.a1.*exp(-((log(frqlist)-a.b1)./a.c1).^2)+spontrate(chan), 'g-', 'LineWidth', 2)
            plot(exp([a.b1 - a.c1, a.b1 + a.c1])./1000, [max(ftc(chan,:)) max(ftc(chan,:))], 'r-', 'Linewidth', 2)
            hold off;
            title([num2str(a.a1), ', ' num2str(a.b1), ', ' num2str(a.c1)])
         end
     catch ME
         chan
         ME
    end
end
