function synuiplotting(app)
    app.DialogueLabel.Text = 'Math...';
    data = app.data;
    if isfield(data, 'sttraces')
        sttracesplot(app)
    else
    if app.ClustsCheck.Value
%The cluster of each spike is stored in data.clusters and data.channels
%data.clusters contains the original cluster code generated in Kilosort
%data.channels contains a new cluster code after sorting the clusters.

%I am putting the sorted cluster code in data.channels because most of our 
%plots groups data by data.channel, and I am too lazy to rewrite all that.
%However, the LFP raster requires the original cluster code (in
%data.clusters).

%The list of clusters selected for analysis is stored in data.chanlist and data.channelsortorder(data.chanlist)
%data.chanlist contains the channel indices selected for plotting, in ascending order
%data.channelsortorder(data.chanlist) contains the actual channel numbers selected for plotting

        [data.chanlist, channelsortorder] = updatechans(app);
        data.channelsortorder = channelsortorder;
        [~, sorti] = sort(channelsortorder);
        data.channels = sorti(data.clusters); 
    else
        data.chanlist = updatechans(app);
    end

    pause(.1)
    data.nchans = max(data.chanlist);
    [spetfreq, spetlevel, stimspets, trials] = freqlevelspet(data.spets, data.stimons, data.frqs, data.lvls);
%spike times aligned to stimulus onsets are stored in stimspets
%The frequency and level of the last stimulus preceding each spike is stored in spetfreq and spetlevel
    binwindow = app.BinEdit.Value./1000; %PSTH bin size (sec)
    analwin = [app.OnsetEdit.Value app.OffsetEdit.Value]./1000; %ftc analysis window;
%tuning curve
    ax = app.FTCAxes; cla(app.FTCAxes, 'reset');
    if app.FTCPopupCheck.Value
        figure(2); clf; ax = axes();
    end
    if numel(unique(data.frqs))>1 || numel(unique(data.lvls)) > 1 %not LDS
        spontwin = [-0.2 0];
        spontrates = findspontrate(data.spets, data.stimons, data.channels, data.chanlist, spontwin);
        totaltime = (spontwin(2)-spontwin(1))*numel(data.stimons)/2;
        spontsig = poissinv(0.95, spontrates);
        spontsig = ((spontsig-spontrates)*2.5+1+spontrates)./totaltime;
        spontrates = spontrates./totaltime;
    end
    labeledfra = [];
    fratable = table([]);
    ftctable = table([]);
    if numel(unique(data.frqs))>1 && numel(unique(data.lvls)) > 1 %FRA
        fra = multifra(stimspets, spetfreq, spetlevel, analwin, data.channels, trials, data.chanlist, data.frqs, data.lvls, ax);
        for chan = 1:size(fra,1)
            chanfra = squeeze(fra(chan,:,:));
            [BLi(chan, 1), BFi(chan, 1)] = ind2sub(size(chanfra), find(chanfra == max(chanfra(:)), 1));
            peakrate(chan,1) = max(chanfra(:));
        end
        for freq = 1:size(fra, 3)
            ftcstats(:,freq,:) = ftcthreshold(fra(:,:,freq), sort(unique(data.lvls)), spontrates, spontsig);
        end
        realfra = permute(fra, [2 3 1]);
        Rate = realfra(:);
        frachannel = repmat(data.chanlist', numel(unique(data.frqs))*numel(unique(data.lvls)), 1);
        Channel = frachannel(:);
        frqlist = sort(unique(data.frqs));
        frafreq = repmat(frqlist', numel(unique(data.lvls)), numel(data.chanlist));
        Frequency = frafreq(:);
        lvllist = sort(unique(data.lvls));
        Level = repmat(lvllist, numel(unique(data.frqs)) * numel(data.chanlist), 1);
        fratable = table(Channel, Frequency, Level, Rate);
        for chan = 1:size(fra, 1)
            labeledfra(chan,:,:) = [[chan; lvllist], [frqlist'; squeeze(fra(chan,:,:))]];
        end

        
        
%         hold(ax, 'on')
%         plot(ax, 1:numel(unique(data.frqs)), 1+numel(unique(data.lvls)) - ftcstats(:,:,1)./10, 'wo-', 'LineWidth', 4)
%         hold(ax, 'off')
        [~,CFi] = min(ftcstats(:,:,1), [], 2);
        frqlist = sort(unique(data.frqs));
        lvllist = sort(unique(data.lvls));
        BL = lvllist(BLi);
        BF = frqlist(BFi);
        CF = frqlist(CFi);
        allthresh = ftcstats(:,:,1);
        thresh = diag(squeeze(ftcstats(:,CFi,1)));
        if app.ClustsCheck.Value
            ftctable = table(data.channelsortorder(data.chanlist), peakrate, BF./1000, BL, CF./1000, thresh, allthresh, spontsig, ...
                'VariableNames', {'Cluster', 'PeakRate_Hz', 'BestFreq_kHz', 'BestLevel_dB', 'CharFreq_kHz', 'Thresh_CF_dB', 'AllThresh', 'ThresholdRate_Hz'});
        else
            ftctable = table(data.chanlist, peakrate, BF./1000, BL, CF./1000, thresh, allthresh, spontsig, ...
                'VariableNames', {'Channel', 'PeakRate_Hz', 'BestFreq_kHz', 'BestLevel_dB', 'CharFreq_kHz', 'Thresh_CF_dB', 'AllThresh', 'ThresholdRate_Hz'});
        end
        for i = 1:numel(frqlist)
            threshnames{i} = ['Thresh_' num2str(round(frqlist(i)./1000)), 'kHz_dB'];
        end
        ftctable= splitvars(ftctable, 'AllThresh', 'NewVariableNames', threshnames);
    elseif numel(unique(data.frqs))>1 && numel(unique(data.lvls)) == 1 %ISO-I
        trials = reshape(trials, 1, numel(trials));
        ftc = synftc(stimspets, spetfreq, analwin, data.channels, trials(1:end), data.chanlist, data.frqs, ax);
        [oprate, obf] = max(ftc, [], 2);
        frqlist = sort(unique(data.frqs))./1000;
        ftcstats = ftcthreshold(ftc, sort(unique(data.frqs)), spontrates, spontsig);
        if app.ClustsCheck.Value
            ftctable = table(data.channelsortorder(data.chanlist), oprate, frqlist(obf), ftcstats, 'VariableNames', {'Cluster', 'PeakRate_Obs_Hz', 'BestFrequency_Obs_kHz', 'ftcstats'});
        else
            ftctable = table(data.chanlist, ftcstats, 'VariableNames', {'Channel', 'ftcstats'});
        end
        ftctable = splitvars(ftctable, 'ftcstats', 'NewVariableNames', {'PeakRate_Model_Hz', 'BestFrequency_Model_kHz', 'Bandwidth2STD_Octaves', 'LowerCutoff_kHz', 'UpperCutoff_kHz'});
        set(ax,'XTickLabel', sort(unique(data.frqs))./1000)
        xlabel(ax, 'Frequency (kHz)')
        if app.ClustsCheck.Value                
            set(ax, 'YTickLabel', channelsortorder);
        end
    elseif numel(unique(data.frqs))==1 && numel(unique(data.lvls)) > 1 %ISO-F
        trials = reshape(trials, 1, numel(trials));          
        ftc = synftc(stimspets, spetlevel, analwin, data.channels, trials(1:end), data.chanlist, data.lvls, ax);
        ftcstats = ftcthreshold(ftc, sort(unique(data.lvls)), spontrates, spontsig);
        if app.ClustsCheck.Value
            ftctable = table(data.channelsortorder(data.chanlist), ftcstats(:,4), 'VariableNames', {'Cluster', 'Threshold_dBSPL'});
        else
            ftctable = table(data.chanlist, ftcstats(:,4), 'VariableNames', {'Channel', 'Threshold_dBSPL'});
        end
        set(ax,'XTickLabel', sort(unique(data.lvls)))
        xlabel(ax, 'Level (dB)')
        title(ax, [app.TankEdit.Value newline num2str(app.OnsetEdit.Value) '-' num2str(app.OffsetEdit.Value) ' ms'])
        if app.ClustsCheck.Value                
            set(ax, 'YTickLabel', channelsortorder);
        end
    end
%running spike rate average
    bincount = NaN(round(max(data.spets)/binwindow), data.nchans);
    for chan = data.chanlist'
        [bincount(:, chan), bintime] = hist(data.spets(data.channels==chan), 0.5*binwindow:binwindow:max(data.spets));
    end
    if app.fratescaleEdit.Value == 0
        yscale = max(bincount(:));
    else
        yscale = app.fratescaleEdit.Value*binwindow;
    end
    [sx, sz] = size(bincount);
    frate = bincount./yscale - (1:sz);
    cla(app.rateAxes, 'reset');
    plotfrate(app, frate, bintime, data.stimons, data.stimoffs, data.chanlist);
    if app.ClustsCheck.Value
        set(app.rateAxes, 'YTick', -numel(app.data.cluster):-1);
        set(app.rateAxes, 'YTickLabel', flipud(channelsortorder));
    end
    data.ratetable = table(reshape(bintime, [numel(bintime), 1]), bincount, 'VariableNames', {'BinTime', 'BinCount'});
    app.frateTable.Data = [(1:max(data.chanlist))', max(bincount)'./binwindow];
    app.frateTable.ColumnName = {'Chan', 'MaxRate'};
%PSTH and LFPs
        psthdata = psthlfpplots(app, data);
%Cluster Waveforms
    if app.ClustsCheck.Value
        cluster = app.data.cluster;
        if app.ClustPopupCheck.Value
            figure(2); clf; ax = axes();
        else
            ax = app.Clust1Axes; cla(app.Clust1Axes, 'reset');
        end
        hold(ax, 'on');
        vvec = [1:-1/16:0 0:1/16:1 1:-1/16:0 0:1/16:1]';
        vv = 1:32;
        colmat = [vvec(vv+32) vvec(vv+22) vvec(vv+12)];

        clustampmat = NaN(numel(cluster), 1);
        for i = 1:numel(cluster)
            clustampmat(i,:) = cluster(i).peakChannel2(1);
        end
        peakChannel2 = clustampmat(channelsortorder); %1xn vector of nearest channel
        lfptemp = [zeros(26, size(data.LFP, 2)); data.LFP; zeros(26, size(data.LFP, 2))];
        msnip = NaN(numel(cluster), 51, size(data.LFP, 2));
        for i = 1:numel(data.chanlist)
            clus = channelsortorder(data.chanlist(i));
            snippets = NaN(numel(cluster(clus).spikes), 51, size(data.LFP, 2));
            for spike = 1:numel(cluster(clus).spikes)
                spikewin = round(cluster(clus).spikes(spike) * data.fs) + (1:51);
                snippets(spike,:,:) = lfptemp(spikewin, :);
            end
%             for chan = 1:size(data.LFP, 2)
%                 r = triu(corr(permute(snippets(:,:,chan), [2 1])));
%                 n = numel(cluster(clus).spikes);
%                 corscore(chan) = sum(sum(triu(r)))./(n*(n-1)/2);
%             end
%             corscore
            msnip(i,:,:) = squeeze(mean(snippets, 1));
            plot(ax,((1:size(data.LFP, 2)) + repmat((-25:25)', 1,size(data.LFP, 2))./60), squeeze(msnip(i,:,:))./range(msnip(:)) - i, 'k');
            plot(ax,(cluster(clus).peakChannel2(1) + (-25:25)./60), ...
            msnip(i,:,cluster(clus).peakChannel2(1))./range(msnip(:)) - i, 'Color', colmat(peakChannel2(i),:));
            waveforms(i,:) = squeeze(msnip(i,:,peakChannel2(i)));
        end
        save('waveforms.mat', 'waveforms');
        hold(ax, 'off')
        set(ax, 'YTick', -numel(cluster):-1);
        set(ax, 'YTickLabel', flipud(channelsortorder));
        ylabel(ax, 'Cluster'); xlabel(ax, 'Channel');

        if app.ClustPopupCheck.Value
            figure(3); clf; ax = axes();
        else
            ax = app.Clust2Axes; cla(app.Clust2Axes, 'reset');
        end
        hold(ax, 'on');
        nrows = ceil(sqrt(numel(data.chanlist)));
        for i = 1:numel(data.chanlist)
            j = i-1;
            clus = channelsortorder(data.chanlist(i));
            snippets = NaN(max([1 numel(cluster(clus).spikes)]), 51);
            for spike = 1:numel(cluster(clus).spikes)
                spikewin = round(cluster(clus).spikes(spike) * data.fs) + (1:51);
                snippets(spike,:) = lfptemp(spikewin, peakChannel2(i));
            end
            xoff = floor(j/nrows);
            yoff = nrows - mod(j, nrows);
            if numel(cluster(clus).spikes) > 20
                snipdens = NaN(60, 51);
                yrange = max(max(abs(snippets)));
                binvec = linspace(-yrange, yrange, 60);
                for samp = 1:51
                    snipdens(:,samp) = hist(snippets(:,samp), binvec);
                end
                imagesc(ax, xoff + [0 1], yoff + [-1 0], 1-snipdens./max(snipdens(:)))
                colormap(ax, gray)
            else
                yrange = max(max(abs(snippets)));
                plot(ax, xoff + (1:size(snippets, 2))./size(snippets, 2), yoff -.5 + snippets./yrange./2, 'k')
            end

                plot(ax, xoff + (1:size(msnip, 2))./size(msnip, 2), yoff -.5 + msnip(i,:,peakChannel2(i))./yrange./2, 'r:')
                text(ax, xoff, yoff, num2str(clus));
        end
        axis(ax, [0 nrows 0 nrows])
        hold(ax, 'off');
        set(ax, 'YTick', []); set(ax, 'XTick', []);
    end %clustcheck
    end %not ra4 data
    app.data = data;
    app.DialogueLabel.Text = 'Ready';
    app.data.psth1table = psthdata.psth1table;
    app.data.psth2table = psthdata.psth2table;
    app.data.psth3table = psthdata.psth3table;
    app.data.lfp1table = psthdata.lfp1table;
    app.data.lfp2table = psthdata.lfp2table;
    app.data.lfpstatstable = psthdata.lfpstatstable;
    app.data.statstable = psthdata.statstable;
    app.data.ftctable = ftctable;
    app.data.fra = labeledfra;
    app.data.fratable = fratable;