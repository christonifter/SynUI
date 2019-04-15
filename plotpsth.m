function yrange = plotpsth(app, ax, pst, yrange, nclusts, sy, sz, stims, scalet)
    switch lower(app.PSTHscaleradio.SelectedObject.Text)
        case 'independent'
        case 'x psths'
            yrange = repmat(max(yrange, [], 'omitnan'), numel(pst), 1);
        case 'x all'
            if app.PSTHscaleEdit.Value == 0
                yrange = repmat(max(max(yrange, [], 'omitnan'), [], 'omitnan'), numel(pst), sy*sz);
            else
                yrange = repmat(app.PSTHscaleEdit.Value, numel(pst), sy*sz);
            end
    end
    if numel(ax) == numel(pst)
        for i = 1:numel(pst)
            hold(ax(i), 'on');
            rectangle(ax(i), 'Position', [stims(i,1) -sy*sz  stims(i,2)./scalet sy*sz], 'FaceColor', .8*[1 1 1], 'EdgeColor', .8*[1 1 1]);
%             bintime1 = pst(i).bintime(find(~isnan(pst(i).bintime), 1));
            plot(ax(i), pst(i).bintime./scalet, (pst(i).bincount./yrange(i,:))- (1:(sy*sz)));
            plot(ax(i), pst(i).PSTHspets./scalet, 1-1*pst(i).PSTHchans - pst(i).PSTHclusters./nclusts, 'k.');
            grid(ax(i), 'on'); 
            ylabel(ax(i), 'Channel');
        end
    else
        hold(ax, 'on');
        for i = 1:numel(pst)
            rectangle(ax, 'Position', [i-1+stims(1) -sy*sz  stims(2)./scalet sy*sz], 'FaceColor', .8*[1 1 1], 'EdgeColor', .8*[1 1 1]);
%             bintime1 = pst(i).bintime(find(~isnan(pst(i).bintime), 1));
            plot(ax, i - 1 + pst(i).bintime./scalet, (pst(i).bincount./yrange(i,:))- (1:(sy*sz)));
            plot(ax, i-1 + pst(i).PSTHspets./scalet, 1-1*pst(i).PSTHchans - pst(i).PSTHclusters./nclusts, 'k.');
        end
        grid(ax, 'on');
        ylabel(ax, 'Channel');hold(ax, 'off'); 

        %             title(ax, [app.TankEditField.Value newline 'Freq = ', app.FrequencyDrop.Value ', ', 'Level = ', app.LevelDrop.Value ', '  app.PSTHscaleradio.SelectedObject.Text ' normalization']);
    end
end
