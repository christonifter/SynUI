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
% windows, and each spike's channel, last stim freq, last stim
% level, and last stim trial
    PSTHsel = ~isnan(PSTHspets) & ismember(data.channels, data.chanlist);
    trial = trial(PSTHsel);
    PSTHchans = data.channels(PSTHsel);
    PSTHsubspets = PSTHspets(PSTHsel);
    PSTHfrqs = spetfreq(PSTHsel);
    PSTHlvls = spetlevel(PSTHsel);
    selection = ones(size(PSTHsubspets));
    if targetfreq >= 0
        selection = selection .* (round(PSTHfrqs)==round(targetfreq)); 
    end
    if targetlvl >= 0
        selection = selection .* (round(PSTHlvls) == round(targetlvl));
    end
    pst.PSTHspets = PSTHsubspets(find(selection));
    pst.PSTHchans = PSTHchans(find(selection));
    frqtrials = trial(find(selection));

% create PSTH, and measure max and average bincounts.
    binvec = app.(startfield).Value:psthwindow:app.(endfield).Value;
    bincount = NaN(length(binvec), numel(data.chanlist));
    bintime = NaN(length(binvec), numel(data.chanlist));
    yrange = NaN(numel(data.chanlist), 1);
    averate = NaN(numel(data.chanlist), 1);
    for chan = 1:numel(data.chanlist)
        [bincount(:,chan), bintime(:,chan)] = hist(pst.PSTHspets(pst.PSTHchans==data.chanlist(chan)), binvec);
        yrange(chan) = max(bincount(:, chan));
        averate(chan) = sum(bincount(:, chan))/ (bintime(end, chan)-bintime(1, chan));
    end
    pst.bincount = bincount./psthwindow;
    pst.bintime = bintime;
    stimdur = median(data.stimoffs - data.stimons);
    stims = [-1.*strcmp(app.(reffield).Value, 'Offset').*stimdur, stimdur];
    yrange2 = yrange./psthwindow;
    averate2 = averate;