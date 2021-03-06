function [pst, yrange2, averate2, stims] = ...
    synpstbin2(data, spetfreq, targetfreq, spetlevel, targetlvl, binsize, psthwindow, reference)

% ref is an array with the times (sec) that will be used to align spike times to
% the stimulus onset or offset
% if reference == 0
%     PSTHspets = data.spets;
%     PSTHspets(data.spets<app.(startfield).Value | data.spets>app.(endfield).Value) = NaN;
%     trial = ones(size(PSTHspets));
% else
% align the spikes to the references. We are temporarily shifting times
% back by (PSTH Start Window) amount, to allow us to include negative
% spike times (spike times preceding the reference).
    if isempty(data.spets)
        spettable = [];
    else
        spettable = repmat(data.spets, 1, numel(reference)) - reference';
    end
% align the spikes to the last reference. We then shift times forward by
% the (PSTH Start Window) amount.
    spettable(spettable<psthwindow(1) | spettable>psthwindow(2)) = NaN;
    [PSTHspets, trial] = min(spettable, [], 2, 'omitnan');
% end
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
    binvec = psthwindow(1):binsize:psthwindow(2);
    bincount = NaN(length(binvec), numel(data.chanlist));
    bintime = NaN(length(binvec), numel(data.chanlist));
    yrange = NaN(numel(data.chanlist), 1);
    averate = NaN(numel(data.chanlist), 1);
    for chan = 1:numel(data.chanlist)
        [bincount(:,chan), bintime(:,chan)] = hist(pst.PSTHspets(pst.PSTHchans==data.chanlist(chan)), binvec);
        yrange(chan) = max(bincount(:, chan));
        averate(chan) = sum(bincount(:, chan))/ (bintime(end, chan)-bintime(1, chan));
    end
    pst.bincount = bincount./binsize./numel(reference);
    pst.bintime = bintime;
    stimdur = median(data.stimoffs - data.stimons);
    stims = [data.stimons(1), data.stimoffs(1) - data.stimons(1)] - reference(1);
    yrange2 = yrange./binsize./numel(reference);
    averate2 = averate;
