function out = psthlfpplots(app, data) 
    cla(app.PSTH1Axes, 'reset');
    cla(app.PSTH2Axes, 'reset');
    frqlist = sort(unique(data.frqs));
    lvllist = sort(unique(data.lvls));
    yrange = NaN(2, numel(data.chanlist));
    averate = NaN(2, numel(data.chanlist));
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
    
    chanvec = zeros(32,1);
        for chan = 1:32
            chanfield = ['CheckBox_' num2str(chan)];
            if app.(chanfield).Value
                chanvec(chan) = 1;
            end
        end
    realchanlist = find(chanvec);
        raw = data.raw;
        fny = data.fs/2; 
        fcl = app.HPCFEdit.Value;
        fcu = app.LPCFEdit.Value;
        
        [b,a] = butter(2, [fcl/fny fcu/fny]);
        LFP= squeeze(filter(b,a, raw));
    if numel(unique(data.lvls))> 1 && numel(unique(data.frqs))> 1 
%FRA
        if app.PSTHPopupCheck.Value
            figure(3); clf; ax = axes(); cla(ax, 'reset');
            title(ax, [app.TankEditField.Value newline 'Level = ', app.LevelDrop.Value ', ', app.PSTHscaleradio.SelectedObject.Text ' normalization']);
        else
            ax = app.PSTH1Axes;
        end
        for freqi = 1:numel(frqlist)
            [pst(freqi), yrange(freqi,:), averate(freqi,:), stims] = synpstbin(app, data, 1, spetfreq, frqlist(freqi), spetlevel, targetlevel);
        end
        plotpsth(app, ax, pst, yrange, data.chanlist, stims, 1.2*(app.PSTH1EndEdit.Value - app.PSTH1StartEdit.Value));
        set(ax, 'XTick', (1:numel(pst)) - 1); set(ax, 'XTickLabel', compose('%d', frqlist)); xlabel(ax, 'Frequency');
        
        if app.PSTHPopupCheck.Value
            figure(4); clf; ax = axes(); cla(ax, 'reset');
            title(ax, [app.TankEditField.Value newline 'Freq = ', app.FrequencyDrop.Value ', ', app.PSTHscaleradio.SelectedObject.Text ' normalization']);
        else
            ax = app.PSTH2Axes; 
        end
        for lvli = 1:numel(lvllist)
            [pst(lvli), yrange(lvli, :), averate(lvli,:), stims] = synpstbin(app, data, 1, spetfreq, targetfreq, spetlevel, lvllist(lvli));
        end
        plotpsth(app, ax, pst, yrange, data.chanlist, stims, 1.2*(app.PSTH1EndEdit.Value - app.PSTH1StartEdit.Value)); 
        set(ax, 'XTick', (1:numel(pst)) - 1); set(ax, 'XTickLabel', compose('%d', lvllist)); xlabel(ax, 'Level');
    elseif isempty(data.evons) || app.PSTHCycleCheck.Value == 0 
%No cycle averaging

        for i = 1:2
            [pst(i), yrange(i,:), averate(i,:), stims(i,:)] = synpstbin(app, data, i, spetfreq, str2double(app.FrequencyDrop.Value), spetlevel, targetlevel);
            psthwindow = app.PSTHBinEdit.Value./1000; 
            binfs = 1/psthwindow;
%             figure(5+i);clf; pspectrum(pst(i).bincount(:, data.chanlist), binfs, 'FrequencyLimits', [0 100]);

            if app.PSTHPopupCheck.Value
                figure(2+i); clf; ax(i) = axes(); 
                title(ax(i), [app.TankEditField.Value newline 'Frequency = ' app.FrequencyDrop.Value ' Hz, ' 'Level = ' app.FrequencyDrop.Value ' dB, ' ...
                    app.PSTHscaleradio.SelectedObject.Text ' normalization' newline num2str(numel(pst(i).PSTHspets)) ' spikes, ' ...
                    num2str(sum(data.frqs==str2double(app.FrequencyDrop.Value))) ' trials'])
            else
                axesfield = ['PSTH' num2str(i) 'Axes']; ax(i) = app.(axesfield); 
            end
            cla(ax(i), 'reset');
        end
        baseline = poissinv(app.CLEdit.Value/100, averate(1,:));
        if strcmpi(app.PSTHbaseline.SelectedObject.Text, 'mean'); baseline = averate(1,:); end
        yrange = plotpsth(app, ax, pst, yrange, data.chanlist, stims, 1); 
        for i = 1:2

            ylim(ax(i), -1.*[max(data.chanlist)+.5 min(data.chanlist)-1]); ylabel(ax(i), 'Channel')
            readfield = ['PSTH' num2str(i) 'RefDrop']; xlabel(ax(i), ['Time relative to ' app.(readfield).Value ' (sec)'])
            hold(ax(i), 'on');
            plot(ax(i), [pst(i).bintime(find(~isnan(pst(i).bintime), 1)) pst(i).bintime(end)], ([baseline; baseline]./yrange(i,:)) - ...
            data.chanlist', 'r:');
            hold(ax(i), 'off');
            field1 = ['PSTH' num2str(i) 'StartEdit']; field2 = ['PSTH' num2str(i) 'EndEdit'];
            xlim(ax(i), [app.(field1).Value app.(field2).Value]);
        %LFP
            if app.PSTHPopupCheck.Value
                figure(4+i); clf; lax(i) = axes(); 
            else
                axesfield = ['LFP' num2str(i) 'Axes']; lax(i) = app.(axesfield);
            end
            cla(lax(i), 'reset'); hold(lax(i), 'on');
            reffield = ['PSTH' num2str(i) 'RefDrop']; startfield = ['PSTH' num2str(i) 'StartEdit']; endfield = ['PSTH' num2str(i) 'EndEdit'];
            if strcmpi(app.(reffield).Value, 'onset')
                ref = data.stimons(1);
            else
                ref = data.stimoffs(1);
            end
            analwin = [app.(startfield).Value, app.(endfield).Value] + ref;
            %datawin tracks where the recording starts and ends,
            %relative to the reference. We can't reference sample points
            %outside of this window in LFP.
            datawin = [1 size(LFP, 1)]./data.fs - ref;
            
            xstart = max([datawin(1); app.(startfield).Value]);
            xend = min([datawin(2); app.(endfield).Value]);
            lfp(i).x = (xstart:(1/data.fs):xend)';

            yi1 = max([1; round(analwin(1).*data.fs)]); yi2 = min([yi1 + numel(lfp(i).x) - 1; size(LFP, 1)]);
            lfp(i).y = LFP(yi1:yi2,realchanlist);
            plot(lax(i), lfp(i).x, lfp(i).y./range(reshape(lfp(i).y, numel(lfp(i).y), 1)) - realchanlist', 'k');
            subspikei = find(data.spets>yi1/data.fs & data.spets<yi2/data.fs);
            subspikes = data.spets(subspikei);
            if app.ClustsCheck.Value
                subclusters = data.clusters(subspikei);
                pC = data.peakChannel(subclusters,1);
                plot(lax(i), subspikes-ref, -pC , 'ro')
            end
            hold(lax(i), 'off');
        end %for PSTH1 and 2
        for chani = 1:numel(data.chanlist)'
           text(lax(1), app.PSTH1StartEdit.Value, .1 - data.chanlist(chani), num2str(baseline(chani), '%0.2f'));
        end

    else 
%Average cycles
        period = mode(diff(data.evons));
        stimdur = mode(data.evoffs - data.evons);
        psthwindow = app.PSTHBinEdit.Value./1000; 
        if psthwindow > period
            app.PSTHBinEdit.Value = period*1000;
            psthwindow = period;
        end

        binvec = 0:psthwindow:period;
            
        for i = 1:2

            [tempst(i), ~, ~, ~] = synpstbin(app, data, i, spetfreq, str2double(app.FrequencyDrop.Value), spetlevel, targetlevel);
            psthwindow = app.PSTHBinEdit.Value./1000; 
            binfs = 1/psthwindow;
            figure(5+i);clf; pspectrum(tempst(i).bincount, binfs, 'FrequencyLimits', [0 100]);
            reffield = ['PSTH' num2str(i) 'RefDrop']; 
            startfield = ['PSTH' num2str(i) 'StartEdit']; 
            endfield = ['PSTH' num2str(i) 'EndEdit'];
            if strcmpi(app.(reffield).Value, 'onset')
                ref = data.stimons(1); 
            else
                ref = data.stimoffs(1); 
            end
            ncycles = sum((data.evons > (ref + app.(startfield).Value)) & (data.evons < (ref + app.(endfield).Value)));
            selection = (data.spets > (ref + app.(startfield).Value)) & (data.spets < (ref + app.(endfield).Value));
            selchannels = data.channels(selection); 
            [~, ~, PSTHspets, ~] = freqlevelspet(data.spets(selection), data.evons, data.frqs, data.lvls); %align spike times to stim pulse times
            PSTHsel = ~isnan(PSTHspets) & ismember(selchannels, data.chanlist);
            pst(i).PSTHchans = selchannels(PSTHsel); 
            pst(i).PSTHspets = PSTHspets(PSTHsel);
            bincount = NaN(length(binvec), numel(data.chanlist));
            bintime = NaN(length(binvec), numel(data.chanlist));
            yrange2 = NaN(numel(data.chanlist), 1);
            averate2 = NaN(numel(data.chanlist), 1);
            for chan = 1:numel(data.chanlist)
                [bincount(:, chan), bintime(:, chan)] = hist(pst(i).PSTHspets(pst(i).PSTHchans==data.chanlist(chan)), binvec);
                yrange2(chan) = max(bincount(:, chan));
                averate2(chan) = sum(bincount(:,chan)) / (app.(endfield).Value-app.(startfield).Value);
                clustspets = pst(i).PSTHspets(pst(i).PSTHchans==data.chanlist(chan));
                phase = 2 .* pi .* clustspets ./ period;
                VS(chan) = sqrt(sum(sin(phase)).^2 + sum(cos(phase)).^2)/numel(phase);
            end
            VS2(:,i) = VS;

            yrange3(i,:) = yrange2./psthwindow./ncycles;
            pst(i).bincount = bincount./psthwindow./ncycles; %bincount reports number of spikes per second per cycle
            pst(i).bintime = bintime; %is bintime == binvec?
            if app.PSTHPopupCheck.Value
                figure(2+i); clf; ax(i) = axes(); 
                title(ax(i), [app.TankEditField.Value newline 'Frequency = ' app.FrequencyDrop.Value ' Hz, ' 'Level = ' app.FrequencyDrop.Value ' dB, ' ...
                app.PSTHscaleradio.SelectedObject.Text ' normalization' newline num2str(numel(pst(i).PSTHspets)) ' spikes, ' ...
                num2str(sum(data.frqs==str2double(app.FrequencyDrop.Value))) ' trials'])
            else
                axesfield = ['PSTH' num2str(i) 'Axes']; ax(i) = app.(axesfield);
            end
            cla(ax(i), 'reset');
            averate(i,:) = averate2;
        end
        yrange = plotpsth(app, ax, pst, yrange3, data.chanlist, [0 stimdur; 0 stimdur], 1);
        for i = 1:2
            baseline = poissinv(app.CLEdit.Value/100, averate(1,:));
            if strcmpi(app.PSTHbaseline.SelectedObject.Text, 'mean'); baseline = averate(1,:); end
            hold(ax(i), 'on')
            plot(ax(i), [pst(i).bintime(find(~isnan(pst(i).bintime), 1)) pst(i).bintime(end)], ([baseline; baseline]./yrange(i,:)) - ...
            data.chanlist', 'r:');
            axis(ax(i), [0 period -max(data.chanlist) 1-min(data.chanlist)]);
            hold(ax(i), 'off');
    %LFP
            reffield = ['PSTH' num2str(i) 'RefDrop']; startfield = ['PSTH' num2str(i) 'StartEdit']; endfield = ['PSTH' num2str(i) 'EndEdit'];
            if app.LFPPopupCheck.Value
                figure(4+i); clf; lax(i) = axes(); 
            else
                axesfield = ['LFP' num2str(i) 'Axes']; lax(i) = app.(axesfield);
            end
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
            lfpcycle = NaN(cycleperiod, size(LFP, 2), numel(cyclet));
            for cycle = 1:numel(cyclet)
                startsamp = 1+cycletsamp(cycle); 
                endsamp = cycletsamp(cycle)+cycleperiod;
                lfpcycle(:,:,cycle) = LFP(startsamp:endsamp,:);
            end
            lfp(i).x = (1:cycleperiod)'./data.fs;
            lfpmeans = mean(lfpcycle, 3); 
            lfp(i).y = lfpmeans(:, realchanlist);
            plot(lax(i), lfp(i).x,  lfp(i).y./range(reshape(lfp(i).y, numel(lfp(i).y), 1)) - realchanlist');
        end
        for chani = 1:numel(data.chanlist)
             text(ax(1), 0, .1 - data.chanlist(chani), num2str(baseline(chani), '%0.2f'));
        end


    end %if PSTH condition
    if app.ClustsCheck.Value
        for i = 1:numel(ax)
            set(ax(i), 'YTick', -numel(app.data.cluster):-1);
            set(ax(i), 'YTickLabel', flipud(data.channelsortorder));
        end
    end

    out.psth1table = table(pst(1).bintime(:,1), pst(1).bincount, 'VariableNames', {'BinTime', 'BinCount'});
    out.psth2table = table(pst(2).bintime(:,1), pst(2).bincount, 'VariableNames', {'BinTime', 'BinCount'});
    out.lfp1table = table(lfp(1).x, lfp(1).y, 'VariableNames', {'Time', 'LFPuV'});
    out.lfp2table = table(lfp(2).x, lfp(2).y, 'VariableNames', {'Time', 'LFPuV'});
%find peak rate/latency   
    psthwindow = app.PSTHBinEdit.Value./1000;
    CI5(1,:) = poissinv(5/100, averate(1,:));
    CI5(2,:) = poissinv(5/100, averate(2,:));
    CI95(1,:) = poissinv(95/100, averate(1,:));
    CI95(2,:) = poissinv(95/100, averate(2,:));
    [PSTH2pmax, PSTH2plati] = max(pst(2).bincount);
    PSTH2plat = pst(2).bintime(PSTH2plati);
%threshold rates
    threshpsth2 = pst(2).bincount > baseline;
    ADtotaldur = sum(threshpsth2) * psthwindow;
    ADtotalspikecount = sum((pst(2).bincount - baseline) .* threshpsth2) .* psthwindow;
    for chani = 1:size(threshpsth2, 2)
        [u, d] = findcross(threshpsth2(:,chani) - 0.5);
        d = d +1;
        if threshpsth2(1, chani)
            u = [1;u];
        end
        if threshpsth2(end, chani)
            d = [d;size(threshpsth2, 1)];
        end
        if numel(d) < 1
            ADmaxcontdur(chani, 1) = 0;
            ADonset(chani,1) = NaN;
            ADoffset(chani,1) = NaN;
            ADmaxcontspikecount(chani, 1) = 0;
        else
            [ADmaxcontdur(chani, 1), ADmaxind] = max(d-u);
            ADonset(chani,1) = u(ADmaxind);
            ADoffset(chani,1) = d(ADmaxind);
            ADmaxcontspikecount(chani, 1) = sum(pst(2).bincount(u(ADmaxind):d(ADmaxind), chani) - baseline(chani)) .* psthwindow;
        end
    end
    ADndiff = NaN(numel(data.chanlist), 1);
    ADdiff = NaN(numel(data.chanlist), 1);
    ADz = NaN(numel(data.chanlist), 1);
    ADndiffl = NaN(numel(data.chanlist), 1);
    ADdiffl = NaN(numel(data.chanlist), 1);
    ADzl = NaN(numel(data.chanlist), 1);

    if strcmpi(app.PSTH1RefDrop.Value, 'onset')
        analwin(1,:) = [app.PSTH1StartEdit.Value, app.PSTH1EndEdit.Value] + data.stimons(1);
    else
        analwin(1,:) = [app.PSTH1StartEdit.Value, app.PSTH1EndEdit.Value] + data.stimoffs(1);
    end
    if strcmpi(app.PSTH2RefDrop.Value, 'onset')
        analwin(2,:) = [app.PSTH2StartEdit.Value, app.PSTH2EndEdit.Value] + data.stimons(1);
    else
        analwin(2,:) = [app.PSTH2StartEdit.Value, app.PSTH2EndEdit.Value] + data.stimoffs(1);
    end

    for loop = 1:2
        if loop == 2
            analwin(1,:) = [0 data.stimons(1)];
            analwin(2,:) = [data.stimoffs(1) size(data.LFP, 1)/data.fs];
        end
        for i = 1:numel(data.chanlist)
            if app.ClustsCheck.Value
                spikesi = app.data.cluster(data.channelsortorder(data.chanlist(i))).spikes;
                binname{i} = ['BinCount_' num2str(data.channelsortorder(data.chanlist(i)))];
            else
                spikesi = data.spets(data.channels == data.chanlist(i));
            end
            postcount = sum(spikesi>analwin(2,1) & spikesi<analwin(2,2));
            posttime = analwin(2,2) - analwin(2,1);
            postrate = postcount/posttime;
            precount = sum(spikesi>analwin(1,1) & spikesi<analwin(1,2));
            pretime = analwin(1,2) - analwin(1,1);
            prerate = precount/pretime;
            if postrate + prerate == 0
                den = 1;
            else
                den = sqrt(postrate + prerate);
            end
            ADndiff(loop,i) = (postrate - prerate) / prerate;
            ADdiff(loop,i) = postrate - prerate;
            ADz(loop,i) = (postrate - prerate) / den;
        end
    end
    out.statstable = table(data.chanlist, averate(1,:)', averate(2,:)', ...
        PSTH2pmax', PSTH2plat', ADdiff(1,:)', 100.*ADndiff(1,:)', ADz(1,:)', ADdiff(2,:)', 100.*ADndiff(2,:)', ADz(2,:)', ...
        baseline', ADtotaldur', ADtotalspikecount', ADmaxcontdur, ADonset, ADoffset, ADmaxcontspikecount, ...
        'VariableNames', {'Channel', 'AveRate_PSTH1_Hz', 'AveRate_PSTH2_Hz', 'MaxRate_PSTH2_Hz', 'MaxRateLatency_PSTH2_ms', ...
        'ADDiff_PSTH_Hz', 'ADPercChange_PSTH_percent', 'ADzscore_PSTH', 'ADDiff_Total_Hz', 'ADPercChange_Total_percent', 'ADzscore_Total', ...
        'ThreshHoldRate_Hz', 'ADDuration_Total_sec', 'ADSpikeCount_Total', ...
        'ADDuration_cont_sec', 'ADOnset_cont_sec', 'ADOffset_cont_sec', 'ADSpikeCount_Cont'});
    if app.ClustsCheck.Value
        cluster = app.data.cluster;
        clustampmat = NaN(numel(cluster), 1);
        for i = 1:numel(cluster)
            clustampmat(i,:) = cluster(i).peakChannel2(1);
        end
        peakChannel2 = clustampmat(data.channelsortorder(data.chanlist)); %1xn vector of nearest channel

        out.psth1table = splitvars(out.psth1table, 'BinCount', 'NewVariableNames', binname);
        out.psth2table = splitvars(out.psth2table, 'BinCount', 'NewVariableNames', binname);
        out.statstable = addvars(out.statstable, data.channelsortorder(data.chanlist), 'After', 1, 'NewVariableNames', {'Cluster'});
        out.statstable = removevars(out.statstable, {'Channel'});
        out.statstable = addvars(out.statstable, peakChannel2, 'After', 1, 'NewVariableNames', {'Nearest_Channel'});
    end
    
    
    if ~isempty(VS2)
        out.statstable = addvars(out.statstable, VS2(:,1), VS2(:,2), 'NewVariableNames', {'VectorStrength1', 'VectorStrength2'});
    end