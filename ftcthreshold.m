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
         if numel(frqlist) < 10 %RLN
            a = fit(frqlist, ftc(chan, :)', 'linearinterp');
            model = a(1:100);
            [u, d] = findcross(model-threshcrit(chan));
            if ~isempty(u)
                thresh = u(1);
            else
                thresh = NaN;
            end
%             figure(10);
% %             if nchans < 50
% %                 nrows = ceil(sqrt(nchans));
% %                 ncols = ceil(nchans./nrows);
% %                 subplot(nrows, ncols, chan)
% %             end
%             plot(frqlist,ftc(chan, :)', 'b')
%             hold on;
%             plot(1:100, model, 'g')
%             plot([1 100], [threshcrit(chan) threshcrit(chan)], 'm--')
%             plot([thresh thresh], [0 max(ftc(chan,:))], 'r--', 'Linewidth', 4)
%             hold off;
            ftcstats(chan) = thresh;
        else %Iso-I
            fo = fitoptions('Method', 'NonlinearLeastSquares', 'Lower', [0 0 0], 'Upper', [1E4 100 1E4]);
            a = fit(log(frqlist), ftc(chan, :)'-spontrate(chan), 'gauss1', fo);
            ftcstats(chan, :) = [a.a1, exp(a.b1)./1000, 2*a.c1/log(2), exp([a.b1 - a.c1, a.b1 + a.c1])./1000];
            
            figure(10);
            if nchans < 50
                nrows = ceil(sqrt(nchans));
                ncols = ceil(nchans./nrows);
                subplot(nrows, ncols, chan)
            end
            plot(frqlist./1000, ftc(chan, :)', 'bo');
            hold on;
            plot(1:70, a.a1.*exp(-((log((1:70).*1000)-a.b1)./a.c1).^2)+spontrate(chan), 'g-', 'LineWidth', 2)
            plot(exp([a.b1 - a.c1, a.b1 + a.c1])./1000, [max(ftc(chan,:)) max(ftc(chan,:))], 'r-', 'Linewidth', 2)
            hold off;
            title([num2str(a.a1), ', ' num2str(a.b1), ', ' num2str(a.c1)])
         end
     catch ME
         chan
         ME
     end
end
