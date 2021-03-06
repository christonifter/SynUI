function out = psthlfpplots(app, data, nchans, nclusts) 
    cla(app.PSTH1Axes, 'reset');
    cla(app.PSTH2Axes, 'reset');
    frqlist = sort(unique(data.frqs));
    lvllist = sort(unique(data.lvls));
    yrange = NaN(2, nchans*nclusts);
    averate = NaN(2, nchans*nclusts);
    VS2 = [];
    targetfreq = -1;
    if ~isnan(str2double(app.FrequencyDrop.Value))
        targetfreq = str2double(app.FrequencyDrop.Value);
    end
    targetlevel = -1;
    if ~isnan(str2double(app.LevelDrop.Value))
        targetlevel = str2double(app.LevelDrop.Value);
    end
    [spetfreq, spetlevel, stimspets, trials] = freqlevelspet(data.spets, data.stimons, data.frqs, data.lvls);
    lfp(1).x = [];
    lfp(2).x = [];
    lfp(1).y = [];
    lfp(2).y = [];
    
    if numel(unique(data.lvls))> 1 && numel(unique(data.frqs))> 1 
%FRA
        if app.PSTHPopupCheck.Value
            figure(3); clf; ax = axes(); cla(ax, 'reset');
        else
            ax = app.PSTH1Axes;
        end
        for freqi = 1:numel(frqlist)
            [pst(freqi), yrange(freqi,:), averate(freqi,:), stims] = synpstbin(app, data, 1, spetfreq, frqlist(freqi), spetlevel, targetlevel);
        end
        plotpsth(app, ax, pst, yrange, max(data.clusters), nclusts, nchans, stims, 1.2*(app.PSTH1EndEdit.Value - app.PSTH1StartEdit.Value));
        set(ax, 'XTick', (1:numel(pst)) - 1); set(ax, 'XTickLabel', compose('%d', frqlist)); xlabel(ax, 'Frequency');
        title(ax, [app.TankEditField.Value newline 'Level = ', app.LevelDrop.Value ', ', app.PSTHscaleradio.SelectedObject.Text ' normalization']);

        if app.PSTHPopupCheck.Value; figure(4); clf; ax = axes(); cla(ax, 'reset');
        else; ax = app.PSTH2Axes; end
        for lvli = 1:numel(lvllist)
            [pst(lvli), yrange(lvli, :), averate(lvli,:), stims] = synpstbin(app, data, 1, spetfreq, targetfreq, spetlevel, lvllist(lvli));
        end
        plotpsth(app, ax, pst, yrange, max(data.clusters), nclusts, nchans, stims, 1.2*(app.PSTH1EndEdit.Value - app.PSTH1StartEdit.Value)); 
        set(ax, 'XTick', (1:numel(pst)) - 1); set(ax, 'XTickLabel', compose('%d', lvllist)); xlabel(ax, 'Level');
        title(ax, [app.TankEditField.Value newline 'Freq = ', app.FrequencyDrop.Value ', ', app.PSTHscaleradio.SelectedObject.Text ' normalization']);
    elseif isempty(data.evons) || app.PSTHCycleCheck.Value == 0 
%No cycle averaging

        for i = 1:2
            [pst(i), yrange(i,:), averate(i,:), stims(i,:)] = synpstbin(app, data, i, spetfreq, str2double(app.FrequencyDrop.Value), spetlevel, targetlevel);
            if app.PSTHPopupCheck.Value; figure(2+i); clf; ax(i) = axes(); 
            else; axesfield = ['PSTH' num2str(i) 'Axes']; ax(i) = app.(axesfield); end
            cla(ax(i), 'reset');
        end
        baseline = poissinv(app.CLEdit.Value/100, averate(1,:));
        if strcmpi(app.PSTHbaseline.SelectedObject.Text, 'mean'); baseline = averate(1,:); end
        yrange = plotpsth(app, ax, pst, yrange, max(data.clusters), nclusts, nchans, stims, 1); 
        for i = 1:2
%             yrange(i,:) = plotpsth(app, ax(i), pst(i), yrange, max(data.clusters), nclusts, nchans, stims(i,:), 1); 
            ylim(ax(i), -1.*[max(data.chanlist)+.5 min(data.chanlist)-1]); ylabel(ax(i), 'Channel')
            readfield = ['PSTH' num2str(i) 'RefDrop']; xlabel(ax(i), ['Time relative to ' app.(readfield).Value ' (sec)'])
            title(ax(i), [app.TankEditField.Value newline 'Frequency = ' app.FrequencyDrop.Value ' Hz, ' 'Level = ' app.FrequencyDrop.Value ' dB, ' ...
                app.PSTHscaleradio.SelectedObject.Text ' normalization' newline num2str(numel(pst(i).PSTHspets)) ' spikes, ' ...
                num2str(sum(data.frqs==str2double(app.FrequencyDrop.Value))) ' trials'])
            hold(ax(i), 'on');
            plot(ax(i), [pst(i).bintime(find(~isnan(pst(i).bintime), 1)) pst(i).bintime(end)], ([baseline; baseline]./yrange(i,:)) - ...
            (1:max(data.clusters)*max(data.chanlist))./max(data.clusters), 'r:');
            hold(ax(i), 'off');
            field1 = ['PSTH' num2str(i) 'StartEdit']; field2 = ['PSTH' num2str(i) 'EndEdit'];
            xlim(ax(i), [app.(field1).Value app.(field2).Value]);
        %LFP
            if app.PSTHPopupCheck.Value; figure(4+i); clf; lax(i) = axes(); 
            else; axesfield = ['LFP' num2str(i) 'Axes']; lax(i) = app.(axesfield); end
            cla(lax(i), 'reset');
            reffield = ['PSTH' num2str(i) 'RefDrop']; startfield = ['PSTH' num2str(i) 'StartEdit']; endfield = ['PSTH' num2str(i) 'EndEdit'];
            if strcmpi(app.(reffield).Value, 'onset')
                analwin = [app.(startfield).Value, app.(endfield).Value] + data.stimons(1);
                %datawin tracks where the recording starts and ends,
                %relative to the reference. We can't reference sample points
                %outside of this window in data.LFP.
                datawin = [1 size(data.LFP, 1)]./data.fs - data.stimons(1);
            else
                analwin = [app.(startfield).Value, app.(endfield).Value] + data.stimoffs(1); 
                datawin = [1 size(data.LFP, 1)]./data.fs - data.stimoffs(1);
            end
            xstart = max([datawin(1); app.(startfield).Value]);
            xend = min([datawin(2); app.(endfield).Value]);
            lfp(i).x = (xstart:(1/data.fs):xend)';
            yi1 = max([1; round(analwin(1).*data.fs)]); yi2 = min([yi1 + numel(lfp(i).x) - 1; size(data.LFP, 1)]);
            lfp(i).y = data.LFP(yi1:yi2,data.chanlist);
            plot(lax(i), lfp(i).x, lfp(i).y./range(reshape(lfp(i).y, numel(lfp(i).y), 1)) - data.chanlist');
        end %for PSTH1 and 2
%         for chani = data.chanlist'
%             for j = 1:nclusts
%                 text(ax(1), app.PSTH1StartEdit.Value, 1.1 - j/nclusts - chani, num2str(baseline(1,(chani*nclusts)+j-1), '%0.2f'));
%             end
%         end

    else 
%Average cycles
        period = mode(diff(data.evons));
        stimdur = mode(data.evoffs - data.evons);
%                 t = (1:size(data.raw, 1))./data.fs;
        psthwindow = app.PSTHBinEdit.Value./1000; 
        binvec = 0:psthwindow:period;
            
        for i = 1:2
            [tempst(i), ~, ~, ~] = synpstbin(app, data, i, spetfreq, str2double(app.FrequencyDrop.Value), spetlevel, targetlevel);
            reffield = ['PSTH' num2str(i) 'RefDrop']; startfield = ['PSTH' num2str(i) 'StartEdit']; endfield = ['PSTH' num2str(i) 'EndEdit'];
            if strcmpi(app.(reffield).Value, 'onset');  ref = data.stimons(1); 
            else; ref = data.stimoffs(1); end
            ncycles = sum((data.evons > (ref + app.(startfield).Value)) & (data.evons < (ref + app.(endfield).Value)));
            selection = (data.spets > (ref + app.(startfield).Value)) & (data.spets < (ref + app.(endfield).Value));
            selchannels = data.channels(selection); selclusters = data.clusters(selection);
            [~, ~, PSTHspets, ~] = freqlevelspet(data.spets(selection), data.evons, data.frqs, data.lvls); %align spike times to stim pulse times
            PSTHsel = ~isnan(PSTHspets) & ismember(selchannels, data.chanlist);
            pst(i).PSTHclusters = selclusters(PSTHsel); pst(i).PSTHchans = selchannels(PSTHsel); pst(i).PSTHspets = PSTHspets(PSTHsel);
            bincount = NaN(length(binvec), max(data.clusters), max(data.chanlist));
            bintime = NaN(length(binvec), max(data.clusters), max(data.chanlist));
            yrange2 = NaN(max(data.clusters), max(data.chanlist));
            averate2 = NaN(max(data.clusters), max(data.chanlist));
            for chan = data.chanlist'
                for cluster = 1:max(data.clusters)
                    [bincount(:, cluster,chan), bintime(:, cluster,chan)] = hist(pst(i).PSTHspets(pst(i).PSTHchans==chan & pst(i).PSTHclusters==cluster), binvec);
                    yrange2(cluster, chan) = max(bincount(:, cluster, chan));
                    averate2(cluster, chan) = sum(bincount(:, cluster, chan)) / (app.(endfield).Value-app.(startfield).Value);
                    clustspets = pst(i).PSTHspets(pst(i).PSTHchans==chan & pst(i).PSTHclusters==cluster);
                    phase = 2 .* pi .* clustspets ./ period;
                    VS(cluster, chan) = sqrt(sum(sin(phase)).^2 + sum(cos(phase)).^2)/numel(phase);
                end
            end
            VS2(:,i) = reshape(VS, numel(VS), 1);

            [sx, nclusts, nchans] = size(bincount);
            yrange3(i,:) = reshape(yrange2./psthwindow./ncycles, nclusts*nchans, 1);
            pst(i).bincount = reshape(bincount./psthwindow./ncycles, sx, nclusts*nchans); %bincount reports number of spikes per second per cycle
            pst(i).bintime = reshape(bintime, sx, nclusts*nchans);
            binfs = 1/psthwindow;
            figure(5+i);clf; pspectrum(pst(i).bincount(:, data.chanlist), binfs, 'FrequencyLimits', [0 100]);

            
            if app.PSTHPopupCheck.Value; figure(2+i); clf; ax(i) = axes(); 
            else; axesfield = ['PSTH' num2str(i) 'Axes']; ax(i) = app.(axesfield); end
            cla(ax(i), 'reset');
            averate(i,:) = reshape(averate2, nclusts*nchans, 1);
        end
        yrange = plotpsth(app, ax, pst, yrange3, max(data.clusters), nclusts, nchans, [0 stimdur; 0 stimdur], 1);
        for i = 1:2
            baseline = poissinv(app.CLEdit.Value/100, averate(1,:));
            if strcmpi(app.PSTHbaseline.SelectedObject.Text, 'mean'); baseline = averate(1,:); end
            hold(ax(i), 'on')
            plot(ax(i), [pst(i).bintime(find(~isnan(pst(i).bintime), 1)) pst(i).bintime(end)], ([baseline; baseline]./yrange(i,:)) - ...
            (1:max(data.clusters)*max(data.chanlist))./max(data.clusters), 'r:');
            axis(ax(i), [0 period -nclusts*nchans 1-min(data.chanlist)]);
            hold(ax(i), 'off');
    %LFP
            reffield = ['PSTH' num2str(i) 'RefDrop']; startfield = ['PSTH' num2str(i) 'StartEdit']; endfield = ['PSTH' num2str(i) 'EndEdit'];
            if app.LFPPopupCheck.Value; figure(4+i); clf; lax(i) = axes(); 
            else; axesfield = ['LFP' num2str(i) 'Axes']; lax(i) = app.(axesfield); end
            cla(lax(i), 'reset');
            if strcmpi(app.(reffield).Value, 'onset')
                analwin = [app.(startfield).Value, app.(endfield).Value] + data.stimons(1);
            else
                analwin = [app.(startfield).Value, app.(endfield).Value] + data.stimoffs(1);
            end
            %The last stimulus cycle will likely be cut off in the middle. This screws up when MATLAB tries to plot that incomplete cycle.
            fullevons = data.evons(1:(end-1)); 
            cyclet = fullevons(fullevons>analwin(1) & fullevons<analwin(2));
            cycletsamp = round(cyclet*data.fs); cycleperiod = min(diff(cycletsamp));
            lfpcycle = NaN(cycleperiod, size(data.LFP, 2), numel(cyclet));
            
            for cycle = 1:numel(cyclet)
                startsamp = 1+cycletsamp(cycle); 
                endsamp = cycletsamp(cycle)+cycleperiod;
                lfpcycle(:,:,cycle) = data.LFP(startsamp:endsamp,:);
            end
            lfp(i).x = (1:cycleperiod)'./data.fs;
            lfpmeans = mean(lfpcycle, 3); 
            lfp(i).y = lfpmeans(:, data.chanlist);
            plot(lax(i), lfp(i).x,  lfp(i).y./range(reshape(lfp(i).y, numel(lfp(i).y), 1)) - data.chanlist');
        end
        for chani = data.chanlist'
            for j = 1:nclusts
                 text(ax(1), 0, 1.1 - j/nclusts - chani, num2str(baseline(1,(chani*nclusts)+j-1), '%0.2f'));
            end
        end



    end %if PSTH condition
    
    out.psth1table = table(pst(1).bintime(:,data.chanlist(1) * nclusts), pst(1).bincount, 'VariableNames', {'BinTime', 'BinCount'});
    out.psth2table = table(pst(2).bintime(:,data.chanlist(1) * nclusts), pst(2).bincount, 'VariableNames', {'BinTime', 'BinCount'});
    out.lfp1table = table(lfp(1).x, lfp(1).y, 'VariableNames', {'Time', 'LFPuV'});
    out.lfp2table = table(lfp(2).x, lfp(2).y, 'VariableNames', {'Time', 'LFPuV'});
    psthwindow = app.PSTHBinEdit.Value./1000;
    CI5(1,:) = poissinv(5/100, averate(1,:));
    CI5(2,:) = poissinv(5/100, averate(2,:));
    CI95(1,:) = poissinv(95/100, averate(1,:));
    CI95(2,:) = poissinv(95/100, averate(2,:));
    out.statstable = table((1:(nclusts*nchans))', averate(1,:)', averate(2,:)', CI5(1,:)', CI5(2,:)', CI95(1,:)', CI95(2,:)', ...
        'VariableNames', {'Channel', 'PSTH1AveRate', 'PSTH2AveRate', 'PSTH1CI5', 'PSTH2CI5', 'PSTH1CI95', 'PSTH2CI95'});
    if ~isempty(VS2)
        out.statstable = table((1:(nclusts*nchans))', averate(1,:)', averate(2,:)', CI5(1,:)', CI5(2,:)', CI95(1,:)', CI95(2,:)', VS2(:,1), VS2(:,2), ...
            'VariableNames', {'Channel', 'PSTH1AveRate', 'PSTH2AveRate', 'PSTH1CI5', 'PSTH2CI5', 'PSTH1CI95', 'PSTH2CI95', 'VectorStrength1', 'VectorStrength2'});
    end
    