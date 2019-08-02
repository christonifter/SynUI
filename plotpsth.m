function [yrange, baseline] = plotpsth(app, ax, pst, averate, yrange, chanlist, stims, scalet)

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
    if numel(ax) == numel(pst) %LDS, IsoI, RLN
        for i = 1:numel(pst)
            rectangle(ax(i), 'Position', [stims(i,1) -max(chanlist)  stims(i,2)./scalet max(chanlist)], 'FaceColor', .8*[1 1 1], 'EdgeColor', .8*[1 1 1]);
            hold(ax(i), 'on');
            if ~strcmpi(app.PSTHbaseline.SelectedObject.Text, 'psth') || i ~= 1 || ~(app.PSTHCycleCheck.Value)
                plot(ax(i), pst(i).bintime./scalet, (pst(i).bincount./yrange(i,:))- chanlist', 'Color', [.4 .4 1], 'LineWidth', 1);
            end
            plot(ax(i), pst(i).PSTHspets./scalet, -pst(i).PSTHchans, 'k.');
            
            baseline = poissinv(app.CLEdit.Value/100, averate(1,:));
            switch lower(app.PSTHbaseline.SelectedObject.Text)
                case 'mean'
                    baseline = averate(1,:);
                case 'psth'
                    if app.PSTHCycleCheck.Value
                        baseline = pst(1).bincount;
                    end
            end
            if size(baseline,1) > 1
                plot(ax(i), pst(i).bintime, baseline./yrange(i,:) - chanlist', 'r', 'LineWidth', 1);
            else
                plot(ax(i), [pst(i).bintime(find(~isnan(pst(i).bintime), 1)) pst(i).bintime(end)], ([baseline; baseline]./yrange(i,:)) - ...
                chanlist', 'r:');
            end
            hold(ax(i), 'off');
            grid(ax(i), 'on'); 
            ylabel(ax(i), 'Channel');

        end
        
        if size(baseline,1) == 1
            for chani = 1:numel(chanlist)
               text(ax(1), ~(app.PSTHCycleCheck.Value) .* app.PSTH1StartEdit.Value, .1 - chanlist(chani), num2str(baseline(chani), '%0.2f'));
            end
        end

    else %FRA
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
