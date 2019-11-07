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
            figure(3); clf; ax(1) = axes(); cla(ax(1), 'reset');
            title(ax(1), [app.TankEdit.Value newline 'Level = ', app.LevelDrop.Value ', ', app.PSTHscaleradio.SelectedObject.Text ' normalization']);
        else
            ax(1) = app.PSTH1Axes;
        end
        clear pst yfreq averate
        for freqi = 1:numel(frqlist)
            [pst(freqi), yrange(freqi,:), averate(freqi,:), stims] = synpstbin(app, data, 1, spetfreq, frqlist(freqi), spetlevel, targetlevel);
        end
        plotpsth(app, ax(1), pst, averate, yrange, data.chanlist, stims, 1.2*(app.PSTH1EndEdit.Value - app.PSTH1StartEdit.Value));
        set(ax(1), 'XTick', (1:numel(pst)) - 1); set(ax(1), 'XTickLabel', compose('%d', frqlist)); xlabel(ax(1), 'Frequency');
        
        if app.PSTHPopupCheck.Value
            figure(4); clf; ax(2) = axes(); cla(ax(2), 'reset');
            title(ax(2), [app.TankEdit.Value newline 'Freq = ', app.FrequencyDrop.Value ', ', app.PSTHscaleradio.SelectedObject.Text ' normalization']);
        else
            ax(2) = app.PSTH2Axes; 
        end
        clear pst yfreq averate
        for lvli = 1:numel(lvllist)
            [pst(lvli), yrange(lvli, :), averate(lvli,:), stims] = synpstbin(app, data, 1, spetfreq, targetfreq, spetlevel, lvllist(lvli));
        end
        plotpsth(app, ax(2), pst, averate, yrange, data.chanlist, stims, 1.2*(app.PSTH1EndEdit.Value - app.PSTH1StartEdit.Value)); 
        set(ax(2), 'XTick', (1:numel(pst)) - 1); set(ax(2), 'XTickLabel', compose('%d', lvllist)); xlabel(ax(2), 'Level');
    elseif isempty(data.evons) || app.PSTHCycleCheck.Value == 0 
%No cycle averaging
        if app.ClustsCheck.Value
            for chan = 1:numel(data.chanlist)
                LDSspikecount(chan,1) = sum(data.clusters==data.chanlist(chan) & data.spets>data.stimons(1) & ...
                    data.spets<data.stimoffs(1));
                LDSduration(chan,1) = data.stimoffs(1)-data.stimons(1);
            end
        else
            for chan = 1:numel(data.chanlist)
                LDSspikecount(chan,1) = sum(data.channels==data.chanlist(chan) & data.spets>data.stimons(1) & ...
                    data.spets<data.stimoffs(1));
                LDSduration(chan,1) = data.stimoffs(1)-data.stimons(1);
            end
        end
        LDSaverate = LDSspikecount./LDSduration;

        for i = 1:2
            [pst(i), yrange(i,:), averate(i,:), stims(i,:)] = synpstbin(app, data, i, spetfreq, str2double(app.FrequencyDrop.Value), spetlevel, targetlevel);
            psthwindow = app.PSTHBinEdit.Value./1000; 
            binfs = 1/psthwindow;
%             figure(5+i);clf; pspectrum(pst(i).bincount(:, data.chanlist), binfs, 'FrequencyLimits', [0 100]);

            if app.PSTHPopupCheck.Value
                figure(2+i); clf; ax(i) = axes(); 
                title(ax(i), [app.TankEdit.Value newline 'Frequency = ' app.FrequencyDrop.Value ' Hz, ' 'Level = ' app.FrequencyDrop.Value ' dB, ' ...
                    app.PSTHscaleradio.SelectedObject.Text ' normalization' newline num2str(numel(pst(i).PSTHspets)) ' spikes, ' ...
                    num2str(sum(data.frqs==str2double(app.FrequencyDrop.Value))) ' trials'])
            else
                axesfield = ['PSTH' num2str(i) 'Axes']; ax(i) = app.(axesfield); 
            end
            cla(ax(i), 'reset');
        end
        [yrange, baseline] = plotpsth(app, ax, pst, averate, yrange, data.chanlist, stims, 1); 
        for i = 1:2

            ylim(ax(i), -1.*[max(data.chanlist)+.5 min(data.chanlist)-1]); ylabel(ax(i), 'Channel')
            readfield = ['PSTH' num2str(i) 'RefDrop']; xlabel(ax(i), ['Time relative to ' app.(readfield).Value ' (sec)'])
            field1 = ['PSTH' num2str(i) 'StartEdit']; field2 = ['PSTH' num2str(i) 'EndEdit'];
            xlim(ax(i), [app.(field1).Value app.(field2).Value]);
        %LFP
            if app.LFPPopupCheck.Value
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
                for clust = 1:numel(data.cluster)
                    pC2(clust) = data.cluster(clust).peakChannel2(1);
                end
                subclusters = data.clusters(subspikei);
                pC = pC2(subclusters);
                plot(lax(i), subspikes-ref, -pC , 'ro')
            end
            hold(lax(i), 'off');
        end %for PSTH1 and 2
        
        threshpsth2 = pst(2).bincount > baseline;
        ADtotaldur = sum(threshpsth2) * psthwindow;
        ADtotalspikecount = sum((pst(2).bincount - baseline) .* threshpsth2) .* psthwindow;
        for chani = 1:size(threshpsth2, 2)
            [u, d] = findcross(threshpsth2(:,chani) - 0.5);
            u = u +1;
            d = d + 1;
            if threshpsth2(1, chani)
                u = [1;u];
            end
            if threshpsth2(end, chani)
                d = [d;size(threshpsth2, 1)];
            end
            u2 = u((u-.5)<(30./psthwindow));
            d2 = d((u-.5)<(30./psthwindow));
            if numel(d) < 1
                ADmaxcontdur(chani, 1) = 0;
                ADonset(chani,1) = NaN;
                ADoffset(chani,1) = NaN;
                AUC(chani,1)=0;
                ADaverate(chani, 1) = 0;
            else
                [maxcontbin, ADmaxind] = max(d-u);
                ADmaxcontdur(chani, 1) = maxcontbin .*psthwindow;
                ADonset(chani,1) = (u(ADmaxind)-.5).*psthwindow;
                ADoffset(chani,1) = (d(ADmaxind)-.5).*psthwindow;
                AUC(chani, 1) = sum(pst(2).bincount(u(ADmaxind):(d(ADmaxind)-1), chani) - baseline(chani)) .* psthwindow;
                ADaverate(chani, 1) = mean(pst(2).bincount(u(ADmaxind):(d(ADmaxind)-1), chani));
            end
            
            if numel(d2) < 1
                eADmaxcontdur(chani, 1) = 0;
                eADonset(chani,1) = NaN;
                eADoffset(chani,1) = NaN;
                eAUC(chani,1)=0;
                eADaverate(chani, 1) = 0;
            else
                [maxcontbin, ADmaxind] = max(d2-u2);
                eADmaxcontdur(chani, 1) = maxcontbin .*psthwindow;
                eADonset(chani,1) = (u2(ADmaxind)-.5).*psthwindow;
                eADoffset(chani,1) = (d2(ADmaxind)-.5).*psthwindow;
                eAUC(chani, 1) = sum(pst(2).bincount((u2(ADmaxind):d2(ADmaxind)-1), chani) - baseline(chani)) .* psthwindow;
                eADaverate(chani, 1) = mean(pst(2).bincount((u2(ADmaxind):d2(ADmaxind)-1), chani));
            end
        end
        ADmaxcontdur(ADmaxcontdur<(3.*psthwindow)) = 0;
        ADonset(ADmaxcontdur<(3.*psthwindow)) = NaN;
        ADoffset(ADmaxcontdur<(3.*psthwindow)) = NaN;
        AUC(ADmaxcontdur<(3.*psthwindow)) = 0;
        ADaverate(ADmaxcontdur<(3.*psthwindow)) = 0;

        eADmaxcontdur(eADmaxcontdur<(3.*psthwindow)) = 0;
        eADonset(eADmaxcontdur<(3.*psthwindow)) = NaN;
        eADoffset(eADmaxcontdur<(3.*psthwindow)) = NaN;
        eAUC(eADmaxcontdur<(3.*psthwindow)) = 0;
        eADaverate(eADmaxcontdur<(3.*psthwindow)) = 0;
        
        ADndiff = NaN(numel(data.chanlist), 1);
        ADratio = NaN(numel(data.chanlist), 1);
        ADdiff = NaN(numel(data.chanlist), 1);
        ADz = NaN(numel(data.chanlist), 1);

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
                ADratio(loop,i) = postrate / prerate;
                ADdiff(loop,i) = postrate - prerate;
                ADz(loop,i) = (postrate - prerate) / den;
            end
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

        binvec = 0:psthwindow:(period+psthwindow);
        p = NaN(2,4096, numel(data.chanlist));
        if app.ClustsCheck.Value
            for chan = 1:numel(data.chanlist)
                LDSspikecount(chan,1) = sum(data.clusters==data.chanlist(chan) & data.spets>data.stimons(1) & ...
                    data.spets<data.stimoffs(1));
                LDSduration(chan,1) = data.stimoffs(1)-data.stimons(1);
            end
        else
            for chan = 1:numel(data.chanlist)
                LDSspikecount(chan,1) = sum(data.channels==data.chanlist(chan) & data.spets>data.stimons(1) & ...
                    data.spets<data.stimoffs(1));
                LDSduration(chan,1) = data.stimoffs(1)-data.stimons(1);
            end
        end
        LDSaverate = LDSspikecount./LDSduration;

        for i = 1:2

            [tempst(i), ~, ~, ~] = synpstbin(app, data, i, spetfreq, str2double(app.FrequencyDrop.Value), spetlevel, targetlevel);
            psthwindow = app.PSTHBinEdit.Value./1000; 
            binfs = 1/psthwindow;
            [p(i,:,:), f] = pspectrum(tempst(i).bincount, binfs, 'FrequencyLimits', [0 100]);
            
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
            
            bincount = NaN(length(binvec)-1, numel(data.chanlist));
            bintime = NaN(length(binvec)-1, numel(data.chanlist));
            yrange2 = NaN(numel(data.chanlist), 1);
            averate2 = NaN(numel(data.chanlist), 1);
            for chan = 1:numel(data.chanlist)
                [y, x] = hist(pst(i).PSTHspets(pst(i).PSTHchans==data.chanlist(chan)), binvec);
                bincount(:, chan) = y(1:(end-1));
                bintime(:,chan) = x(1:(end-1));
                pst(i).analcount(chan, 1) = sum(pst(i).PSTHspets(pst(i).PSTHchans==data.chanlist(chan)) > app.PSTHOnsetEdit.Value*1E-3 & ...
                    pst(i).PSTHspets(pst(i).PSTHchans==data.chanlist(chan)) < app.PSTHOffsetEdit.Value*1E-3);
                yrange2(chan) = max(bincount(:, chan));
                averate2(chan) = sum(bincount(:,chan)) / (app.(endfield).Value-app.(startfield).Value);
                clustspets = pst(i).PSTHspets(pst(i).PSTHchans==data.chanlist(chan));
                phase = 2 .* pi .* clustspets ./ period;
                VS(chan) = sqrt(sum(sin(phase)).^2 + sum(cos(phase)).^2)/numel(phase);
                [u, d] = findcross(bincount(:,chan) - 0.5*max(bincount(:,chan))+.1);
                if numel(u) > 0 & numel(d) > 0
                    if u(1)>d(1)
                        u = [1; u];
                    end
                    if u(end)>d(end)
                        d = [d; size(bintime, 1)];
                    end
                    if numel(u) == 1
                        halfmaxdur(chan,i) = (d - u)*app.PSTHBinEdit.Value;
                    else
                        gaps = u(2:end)-d(1:end-1);
                        siggaps = find(gaps>2);
                        if sum(siggaps) == 0
                            halfmaxdur(chan,i) = (d(end) - u(1))*app.PSTHBinEdit.Value;
                        else
                            [~, plat] = max(bincount(:,chan));
                            maxpeak = find(d>=plat, 1);
                            peakonset = bintime(u(siggaps(find(siggaps < maxpeak, 1, 'last'))+1), chan);
                            peakoffset = bintime(d(siggaps(find(siggaps >= maxpeak, 1, 'first'))), chan);
                            if isempty(peakonset)
                                peakonset = bintime(u(1), chan);
                            end
                            if isempty(peakoffset)
                                peakoffset = bintime(d(end), chan);
                            end
                            halfmaxdur(chan,i) = (peakoffset - peakonset)*1000;
                        end
                    end
                else
                    halfmaxdur(chan,i) = NaN;
                end
            end
            VS2(:,i) = VS;

            yrange3(i,:) = yrange2./psthwindow./ncycles;
            pst(i).bincount = bincount./psthwindow./ncycles; %bincount reports number of spikes per second per cycle
            pst(i).bintime = bintime; %is bintime == binvec?
            if app.PSTHPopupCheck.Value
                figure(2+i); clf; ax(i) = axes(); 
                title(ax(i), [app.TankEdit.Value newline 'Frequency = ' app.FrequencyDrop.Value ' Hz, ' 'Level = ' app.FrequencyDrop.Value ' dB, ' ...
                app.PSTHscaleradio.SelectedObject.Text ' normalization' newline num2str(numel(pst(i).PSTHspets)) ' spikes, ' ...
                num2str(sum(data.frqs==str2double(app.FrequencyDrop.Value))) ' trials'])
            else
                axesfield = ['PSTH' num2str(i) 'Axes']; ax(i) = app.(axesfield);
            end
            cla(ax(i), 'reset');
            xlabel(ax(i), 'Time (ms)')
            averate(i,:) = averate2;
        end
        [~, MFi] = min(abs(f - 1/period));
        AC = squeeze(p(:, MFi, :));
        DC = squeeze(p(:, 1, :));
        TCF = AC./DC;

        yrange = plotpsth(app, ax, pst, averate, yrange3, data.chanlist, [0 stimdur; 0 stimdur], 1E-3);

        for i = 1:2
            axis(ax(i), [0 period.*1E3 -max(data.chanlist) 1-min(data.chanlist)]);
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
            lastevbeforeend = max(find(data.evons * data.fs< size(data.raw,1)));
            fullevons = data.evons(1:(lastevbeforeend-1));                 
            cyclet = fullevons(fullevons>analwin(1) & fullevons<analwin(2));
            cycletsamp = round(cyclet*data.fs); cycleperiod = min(diff(cycletsamp));
            lfpcycle = NaN(cycleperiod, size(LFP, 2), numel(cyclet));
            for cycle = 1:numel(cyclet)
                startsamp = 1+cycletsamp(cycle); 
                endsamp = cycletsamp(cycle)+cycleperiod;
                lfpcycle(:,:,cycle) = LFP(startsamp:endsamp,:);
            end
            lfp(i).x = (1:cycleperiod)'.*1E3./data.fs;
            lfpmeans = mean(lfpcycle, 3); 
            lfp(i).y = lfpmeans(:, realchanlist);
            [a,b,c] =lfppeaks(lfp(i).y);
            lfp(i).peaks = a.*1E3;
            lfp(i).lats = b.*1E3./data.fs;
            lfp(i).slopes = c.*data.fs;
            rmswin = max([1 round(app.RMSonEdit.Value*data.fs/1E3)]):min([size(lfp(i).y, 1) round(app.RMSoffEdit.Value*data.fs/1E3)]);
            lfp(i).rms = rms(lfp(i).y(rmswin,:))'.*1E3;
            rectangle(lax(i), 'Position', [app.RMSonEdit.Value -max(realchanlist) app.RMSoffEdit.Value-app.RMSonEdit.Value  max(realchanlist)], 'FaceColor', .8*[.8 1 .8], 'EdgeColor', .8*[1 1 1]);
            hold(lax(i), 'on');
            plot(lax(i), lfp(i).x,  lfp(i).y./range(reshape(lfp(i).y, numel(lfp(i).y), 1)) - realchanlist');
            hold(lax(i), 'off');    
            xlabel(lax(i), 'Time (ms)')
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
    out.psth3table = table([]);
    out.lfp1table = table(lfp(1).x, lfp(1).y, 'VariableNames', {'Time', 'LFPuV'});
    out.lfp2table = table(lfp(2).x, lfp(2).y, 'VariableNames', {'Time', 'LFPuV'});
%find peak rate/latency   
    psthwindow = app.PSTHBinEdit.Value./1000;
    CI5(1,:) = poissinv(5/100, averate(1,:).*psthwindow)./psthwindow;
    CI5(2,:) = poissinv(5/100, averate(2,:).*psthwindow)./psthwindow;
    CI95(1,:) = poissinv(95/100, averate(1,:).*psthwindow)./psthwindow;
    CI95(2,:) = poissinv(95/100, averate(2,:).*psthwindow)./psthwindow;
    
    if app.PSTHCycleCheck.Value
        CI5(1,:) = norminv(5/100, averate(1,:), std(pst(1).bincount));
        CI5(2,:) = norminv(5/100, averate(2,:), std(pst(2).bincount));
        CI95(1,:) = norminv(95/100, averate(1,:), std(pst(1).bincount));
        CI95(2,:) = norminv(95/100, averate(2,:), std(pst(2).bincount));
        CI5(isnan(CI5)) = 0;
        CI95(isnan(CI95)) = 0;
    end
    
    [PSTH1pmax, PSTH1plati] = max(pst(1).bincount);
    PSTH1plat = pst(1).bintime(PSTH1plati);
    [PSTH2pmax, PSTH2plati] = max(pst(2).bincount);
    PSTH2plat = pst(2).bintime(PSTH2plati);
%threshold rates
    if numel(unique(data.lvls))> 1 || numel(unique(data.frqs))> 1
        out.statstable = table([]);
        out.lfpstatstable = table([]);
    else
        out.statstable = table(data.chanlist, averate(1,:)', averate(2,:)', ...
            CI95(1,:)', LDSaverate, LDSspikecount, LDSduration, PSTH1pmax', PSTH1plat', PSTH2pmax', PSTH2plat', ...
            'VariableNames', {'Channel', 'AveRate_PSTH1_Hz', 'AveRate_PSTH2_Hz', ...
            'p95_Confidence_Level', 'LDSaverate', 'LDSspikecount', 'LDSduration', 'MaxRate_PSTH1_Hz', 'MaxRateLatency_PSTH1_sec', ...
            'MaxRate_PSTH2_Hz', 'MaxRateLatency_PSTH2_sec'});
        
        if isempty(VS2)
            out.statstable = addvars(out.statstable, ...
                ADdiff(1,:)', 100.*ADndiff(1,:)', ADratio(1,:)', ADz(1,:)', ...
                ADdiff(2,:)', 100.*ADndiff(2,:)', ADratio(2,:)', ADz(2,:)', ...
            baseline', ADtotaldur', ADtotalspikecount', ADmaxcontdur, ADonset, ADoffset, AUC, ADaverate, ...
            eADmaxcontdur, eADonset, eADoffset, eAUC, eADaverate, ...
            'NewVariableNames', {'ADDiff_PSTH_Hz', 'ADPercChange_PSTH_percent', 'ADRatio_PSTH', 'ADzscore_PSTH', ...
            'ADDiff_Total_Hz', 'ADPercChange_Total_percent', 'ADratio_Total', 'ADzscore_Total',  ...
            'ThreshHoldRate_Hz', 'ADDuration_Total_sec', 'ADSpikeCount_Total', ...
            'ADDuration_cont_sec', 'ADOnset_cont_sec', 'ADOffset_cont_sec', 'ADSpikeCount_Cont', 'ADAveRate_Hz', ...
            'EarlyADDur_sec', 'EarlyADOnset_sec', 'EarlyADOffset_sec', 'EarlySpikeCount', 'EarlyADAveRate_Hz'});
            out.lfpstatstable = [];
        else
            out.psth3table = table(pst(1).bintime(:,1), pst(2).bincount-pst(1).bincount, 'VariableNames', {'BinTime', 'BinCount'});
            hasntmultiplespikes = find(sum(tempst(1).bincount>0) < 1.5 | sum(tempst(2).bincount>0) < 1.5);
            modgaindb = 10*log10(AC(2,:)'./AC(1,:)');
            modgaindb(hasntmultiplespikes) = NaN;
            out.statstable = addvars(out.statstable, pst(1).analcount, pst(2).analcount, halfmaxdur(:,1), halfmaxdur(:,2), ...
                VS2(:,1), VS2(:,2), AC(1,:)', AC(2,:)', DC(1,:)', DC(2,:)', ...
                100.*TCF(1,:)', 100.*TCF(2,:)', modgaindb, 100.*(TCF(2,:)' - TCF(1,:)'), ...
                'NewVariableNames', {'WindowCount_PSTH1', 'WindowCount_PSTH2', 'Duration1_ms', 'Duration2_ms', 'VectorStrength1', 'VectorStrength2', 'SpectralPowerMod1', 'SpectralPowerMod2', ...
                'DCPower1', 'DCPower2', 'TemporalCodingFraction1', 'TemporalCodingFraction2', 'ModGain_dB', 'TCF_Diff'});
            out.lfpstatstable = table(realchanlist, lfp(1).peaks(1,:)', lfp(1).peaks(2,:)', lfp(1).peaks(3,:)', ...
                lfp(2).peaks(1,:)', lfp(2).peaks(2,:)', lfp(2).peaks(3,:)', lfp(1).lats(1,:)', ...
                lfp(1).lats(2,:)', lfp(1).lats(3,:)', lfp(1).lats(7,:)'-lfp(1).lats(6,:)', lfp(2).lats(1,:)',...
                lfp(2).lats(2,:)', lfp(2).lats(3,:)', lfp(2).lats(7,:)'-lfp(2).lats(6,:)', lfp(1).slopes(1,:)', ...
                lfp(1).slopes(2,:)', lfp(1).slopes(3,:)', lfp(1).slopes(4,:)', ...
                lfp(2).slopes(1,:)', lfp(2).slopes(2,:)', lfp(2).slopes(3,:)', lfp(2).slopes(4,:)', lfp(1).rms, lfp(2).rms, 'VariableNames', ...
                {'Channel', 'LFP1MinPeak_mV', 'LFP1MaxPeak1_mV', 'LFP1MaxPeak2_mV', 'LFP2MinPeak_mV', 'LFP2MaxPeak1_mV', ...
                'LFP2MaxPeak2_mV', 'LFP1MinLat_ms', 'LFP1MaxLat1_ms', 'LFP1MaxLat2_ms', 'LFP1PeakDur_ms', 'LFP2MinLat_ms', 'LFP2MaxLat1_ms', ...
                'LFP2MaxLat2_ms',  'LFP2PeakDur_ms', 'LFP1FallSlope_mVpms', 'LFP1RiseSlope_mVpms', 'LFP1AveFallSlope_mVpms', 'LFP1AveRiseSlope_mVpms', ...
                'LFP2FallSlope_mVpms', 'LFP2RiseSlope_mVpms', 'LFP2AveFallSlope_mVpms', 'LFP2AveRiseSlope_mVpms', 'LFP1RMS_mV', 'LFP2RMS_mV'});

            rat1 = squeeze(p(2,:,:)./p(1,:,:));
            rat2 = squeeze(p(2,1,:)./p(1,1,:));
            pchange = 20*log10(rat1')';
%             figure(7)
%             subplot(3,1,1)
%             plot(f,squeeze(p(1,:,:)));
%             title('Spectral Power pre-LDS')
%             subplot(3,1,2)
%             plot(f,squeeze(p(2,:,:)));
%             title('Spectral Power post-LDS')
%             subplot(3,1,3)
%             plot(f,pchange);
%             title('Power Change post/pre-LDS (dB)')
%             xlabel('Frequency (Hz)')
        end
        
        
        if app.ClustsCheck.Value
            cluster = app.data.cluster;
            clustampmat = NaN(numel(cluster), 1);
            for i = 1:numel(cluster)
                clustampmat(i,:) = cluster(i).peakChannel2(1);
            end
            peakChannel2 = clustampmat(data.channelsortorder(data.chanlist)); %1xn vector of nearest channel
            for i = 1:numel(data.chanlist)
                spikesi = app.data.cluster(data.channelsortorder(data.chanlist(i))).spikes;
                binname{i} = ['BinCount_' num2str(data.channelsortorder(data.chanlist(i)))];
            end


            out.psth1table = splitvars(out.psth1table, 'BinCount', 'NewVariableNames', binname);
            out.psth2table = splitvars(out.psth2table, 'BinCount', 'NewVariableNames', binname);
            out.statstable = addvars(out.statstable, data.channelsortorder(data.chanlist), 'After', 1, 'NewVariableNames', {'Cluster'});
            out.statstable = removevars(out.statstable, {'Channel'});
            out.statstable = addvars(out.statstable, peakChannel2, 'After', 1, 'NewVariableNames', {'Nearest_Channel'});
        end
    end