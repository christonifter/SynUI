function plotfrate(app, frate, bintime, stimons, stimoffs, chanlist)
    stimp = [NaN NaN];
    if app.ratePopupCheck.Value
        figure(1); clf;
        ax = axes();
        title(ax, [app.TankEdit.Value newline num2str(app.BinEdit.Value) ' msec bins'])
    else
        ax = app.rateAxes;
    end
    if numel(stimons>0)
        if numel(stimoffs)<numel(stimons)
            stimp = [stimons [stimoffs; max(stimons)]];
        else
            stimp = [stimons stimoffs];
        end
    end
    plot(ax, 0,0)
    hold(ax, 'on');
    if numel(stimoffs > 0)
        for stim = 1:size(stimp, 1)
            rectangle(ax, 'Position', [stimp(stim, 1) -100  stimp(stim, 2)-stimp(stim, 1) 100], 'FaceColor', .8*[1 1 1], 'EdgeColor', .8*[1 1 1])
        end
    end

    plot(ax, bintime, frate)
%                     plot(ax, stimp', [0 0], 'r-', 'LineWidth', 4);
    hold(ax, 'off');
    xlabel(ax, 'Time (sec)')
    ylabel(ax, 'Channel')
    ylim(ax, -1.*[chanlist(end)+.5 chanlist(1)-1])
end