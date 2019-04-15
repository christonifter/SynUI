function synfrate(spets, ratebinwindow, stimons, stimoffs, channels, nchans, ratewin, ax)
if ~exist('ax', 'var')
    ax = 1;
end
    frate = NaN(nchans,round(max(spets)/ratebinwindow));
    for chan = 1:nchans
        [bincount, bintime] = hist(spets(channels==chan), 0.5*ratebinwindow:ratebinwindow:max(spets));
        yrange = max(bincount) - min(bincount);
        frate(chan,:) = (bincount- min(bincount))./(1.1*yrange) - chan;
    end
    plot(ax, 0,0)
    hold(ax, 'on');
        plot(ax, bintime, frate)
    if numel(stimons>0)
        if numel(stimoffs)<numel(stimons)
            stimp = [stimons [stimoffs; max(stimons)]];
        else
            stimp = [stimons stimoffs];
        end
        plot(ax, stimp', [0 0], 'r-', 'LineWidth', 4);
    end
    hold(ax, 'off');

     xlim(ax, [max([0 max(spets)-ratewin]), max([1 max(spets)])]);