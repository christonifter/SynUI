function [pst, yrange2, averate2, stims] = synpstbin(app, data, i, spetfreq, targetfreq, spetlevel, targetlvl)
    psthwindow = app.PSTHBinEdit.Value./1000; 
    reffield = ['PSTH' num2str(i) 'RefDrop'];
    startfield = ['PSTH' num2str(i) 'StartEdit'];
    endfield = ['PSTH' num2str(i) 'EndEdit'];
% ref is an array with the times (sec) that will be used to align spike times to
% the stimulus onset or offset
    if strcmpi(app.(reffield).Value, 'onset')
        ref = data.stimons;
    else
        ref = data.stimoffs;
    end
% align the spikes to the references. We are temporarily shifting times
% back by (PSTH Start Window) amount, to allow us to include negative
% spike times (spike times preceding the reference).
    if isempty(data.spets)
        spettable = [];
    else
        spettable = repmat(data.spets, 1, numel(data.stimons)) - (ref + app.(startfield).Value)';
    end
% align the spikes to the last reference. We then shift times forward by
% the (PSTH Start Window) amount.
    spettable(spettable<0) = NaN;
    spettable(spettable>(app.(endfield).Value - app.(startfield).Value)) = NaN;
    [PSTHspets, trial] = min(spettable + app.(startfield).Value, [], 2, 'omitnan');

% We will save the subset of spikes that fall within the user-defined
% windows, and each spike's channel, cluster, last stim freq, last stim
% level, and last stim trial
    PSTHsel = ~isnan(PSTHspets) & ismember(data.channels, data.chanlist);
    trial = trial(PSTHsel);
    PSTHclusters = data.clusters(PSTHsel);
    PSTHchans = data.channels(PSTHsel);
    PSTHsubspets = PSTHspets(PSTHsel);
    PSTHfrqs = spetfreq(PSTHsel);
    PSTHlvls = spetlevel(PSTHsel);
    selection = ones(size(PSTHsubspets));
    if targetfreq >= 0
        selection = selection .* (PSTHfrqs==targetfreq); 
    end
    if targetlvl >= 0
        selection = selection .* (PSTHlvls == targetlvl);
    end
    pst.PSTHspets = PSTHsubspets(find(selection));
    pst.PSTHchans = PSTHchans(find(selection));
    pst.PSTHclusters = PSTHclusters(find(selection));
    frqtrials = trial(find(selection));

% create PSTH, and measure max and average bincounts.
    binvec = app.(startfield).Value:psthwindow:app.(endfield).Value;
    bincount = NaN(length(binvec), max(data.clusters), max(data.chanlist));
    bintime = NaN(length(binvec), max(data.clusters), max(data.chanlist));
    yrange = NaN(max(data.clusters), max(data.chanlist));
    averate = NaN(max(data.clusters), max(data.chanlist));
    for chan = data.chanlist'
        for cluster = 1:max(data.clusters)
            [bincount(:, cluster,chan), bintime(:, cluster,chan)] = ...
                hist(pst.PSTHspets(pst.PSTHchans==chan & pst.PSTHclusters==cluster), binvec);
            yrange(cluster, chan) = max(bincount(:, cluster, chan));
            averate(cluster, chan) = sum(bincount(:, cluster, chan))/ ...
                (bintime(end, cluster, chan)-bintime(1, cluster, chan));
        end
    end
    [sx, sy, sz] = size(bincount);
    pst.bincount = reshape(bincount./psthwindow, sx, sy*sz);
    pst.bintime = reshape(bintime, sx, sy*sz);
    stimdur = median(data.stimoffs - data.stimons);
    stims = [-1.*strcmp(app.(reffield).Value, 'Offset').*stimdur, stimdur];
    yrange2 = reshape(yrange./psthwindow, sy*sz, 1);
    averate2 = reshape(averate, sy*sz, 1);