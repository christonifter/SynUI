function spontrates = findspontrate(spets, stimons, channels, chanlist, analwin)

%find spontaneous rate from recording during pulsed sound presentation.
%The strategy is to count spikes from a time window before each sound onset.
%inputs: 
%spets = spike times in sec
%stimons = sound onsets in sec
%channels = channel (recording site or cluster) for each spike
%analwin = time window relative to sound onset to collect spikes from.


analwins = (stimons(:) + analwin)';
analwins = analwins(:);
if analwins(1)< 0
    analwins = analwins(3:end);
end
spikestimeindex = discretize(spets, analwins);
spikeisinanalwin = mod(spikestimeindex,2) == 1;
analchans = channels(spikeisinanalwin);
totaltime = (analwin(2)-analwin(1))*numel(analwins)/2;

spontrates = histcounts(analchans, .5+(0:numel(chanlist)))'./ totaltime;
spontrates = spontrates(chanlist);