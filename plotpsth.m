function yrange = plotpsth(app, ax, pst, yrange, chanlist, stims, scalet)
    switch lower(app.PSTHscaleradio.SelectedObject.Text)
        case 'independent'
        case 'x psths'
            yrange = repmat(max(yrange, [], 'omitnan'), numel(pst), 1);
        case 'x all'
            if app.PSTHscaleEdit.Value == 0
                yrange = repmat(max(max(yrange, [], 'omitnan'), [], 'omitnan'), numel(pst), numel(chanlist));
            else
                yrange = repmat(app.PSTHscaleEdit.Value, numel(pst), numel(chanlist));
            end
    end
    if numel(ax) == numel(pst)
        for i = 1:numel(pst)
            hold(ax(i), 'on');
            rectangle(ax(i), 'Position', [stims(i,1) -max(chanlist)  stims(i,2)./scalet max(chanlist)], 'FaceColor', .8*[1 1 1], 'EdgeColor', .8*[1 1 1]);
            plot(ax(i), pst(i).bintime./scalet, (pst(i).bincount./yrange(i,:))- chanlist');
            plot(ax(i), pst(i).PSTHspets./scalet, -pst(i).PSTHchans, 'k.');
            grid(ax(i), 'on'); 
            ylabel(ax(i), 'Channel');
        end
    else
        hold(ax, 'on');
        for i = 1:numel(pst)
            rectangle(ax, 'Position', [i-1+stims(1) -max(chanlist)  stims(2)./scalet max(chanlist)], 'FaceColor', .8*[1 1 1], 'EdgeColor', .8*[1 1 1]);
            plot(ax, i - 1 + pst(i).bintime./scalet, (pst(i).bincount./yrange(i,:))- chanlist');
            plot(ax, i-1 + pst(i).PSTHspets./scalet, -pst(i).PSTHchans, 'k.');
        end
        grid(ax, 'on');
        ylabel(ax, 'Channel');hold(ax, 'off'); 
    end
end
