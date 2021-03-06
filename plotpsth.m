function [yrange, baseline] = plotpsth(app, ax, pst, averate, yrange, chanlist, stims, scalet)
    psthwindow = app.PSTHBinEdit.Value./1000; 
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
            rectangle(ax(i), 'Position', [app.PSTHOnsetEdit.Value -max(chanlist)  app.PSTHOffsetEdit.Value-app.PSTHOnsetEdit.Value max(chanlist)], 'FaceColor', .8*[1 1 .9], 'EdgeColor', .8*[1 1 .9]);
            rectangle(ax(i), 'Position', [app.BaselineStartEdit.Value -max(chanlist)  app.BaselineEndEdit.Value-app.BaselineStartEdit.Value max(chanlist)], 'FaceColor', .8*[.9 1 1], 'EdgeColor', .8*[.9 1 1]);
            if ~app.PSTHButton.Value || i == 1 || ~app.PSTHCycleCheck.Value
                plot(ax(i), pst(i).bintime./scalet, (pst(i).bincount./yrange(i,:))- chanlist', 'Color', [.4 .4 1], 'LineWidth', 1);
            else
                plot(ax(i), pst(i).bintime./scalet, (pst(i).bincount./yrange(i,:))- chanlist', 'r', 'LineWidth', 1);
            end
            plot(ax(i), pst(i).PSTHspets./scalet, -pst(i).PSTHchans, 'k.');
            if app.PSTHCycleCheck.Value
                plot(ax(i), 0.5.*ones(numel(chanlist), 1), .5-chanlist.*app.data.sounddriven(:,i), 'm*');
            end

%Option 1 (recommended by CML) finds confidence level from poisson distribution of # spikes per bin period
%However, the confidence level will vary with bin period
             baseline = poissinv(app.CLEdit.Value/100, averate(1,:).*psthwindow)./psthwindow;

            
%Option 2 finds confidence level from normal distribution with mu = average pre firing rate, sigma = std pre firing rate
%However, spike rates are not gaussian at low spike rates 
            if app.PSTHCycleCheck.Value
                baseline = norminv(app.CLEdit.Value/100, averate(1,:), std(pst(1).bincount));
                baseline(isnan(baseline)) = 0;
            end

%Option 3 finds confidence level from poisson distribution of total #
%spikes from pre period
%             baseline = poissinv(app.CLEdit.Value/100, sum(pst(1).bincount)).*psthwindow./(pst(1).bintime(end,:) - pst(1).bintime(1,:));
            switch lower(app.PSTHbaseline.SelectedObject.Text)
                case 'mean'
                    baseline = averate(1,:);
                case 'psth 1'
                    if app.PSTHCycleCheck.Value
                        baseline = pst(1).bincount;
                    end
            end
            if size(baseline,1) > 1
                plot(ax(i), pst(i).bintime./scalet, baseline./yrange(i,:) - chanlist', 'Color', [.4 .4 1], 'LineWidth', 1);
                if app.SubtractionCheck.Value
                    plot(ax(i), pst(i).bintime./scalet, (pst(i).bincount - baseline)./yrange(i,:) - chanlist', 'Color', [0 0.6 0], 'LineWidth', 1);
                end
            else
                plot(ax(i), [pst(i).bintime(find(~isnan(pst(i).bintime), 1)) pst(i).bintime(end)]./scalet, [baseline; baseline]./yrange(i,:) - ...
                chanlist', 'r:');
            end
            hold(ax(i), 'off');
%             grid(ax(i), 'on'); 
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
