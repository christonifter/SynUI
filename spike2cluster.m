function [clusters, spikes, experiment] = spike2cluster(rez, experiment)
spikeTimes     = rez.st3(:,1);
spikeClusters = rez.st3(:,2);
peakChannel = rez.iNeighPC';
peakChannel2 = NaN(numel(experimentsi), numel(unique(spikeClusters)),32);
for i = 1:numel(experiment)
    tankpath = [inpath subject char(folders(experimentsi(i)))];
    data = readtank(tankpath);
    spikei = (spikeTimes > segmentstarts(i) & spikeTimes < segmentends(i));
    spikes = (spikeTimes(spikei) - segmentstarts(i) + 1)./fs;
    clusters = spikeClusters(spikei);
    lfptemp = [zeros(26, 32); data.LFP; zeros(26, 32)];
    clear cluster
    for j = 1:max(spikeClusters)
        clind = spikei & (spikeClusters == j);
        experiment(i).cluster(j).spikes = (spikeTimes(clind) - segmentstarts(i) + 1)./fs;
        experiment(i).cluster(j).peakChannel = peakChannel(j,:);
        snippets = NaN(numel(experiment(i).cluster(j).spikes), 51 ,32); %snippets of spikes from this experiment, all clusters
        for spike = 1:numel(experiment(i).cluster(j).spikes)
            spikewin = round(experiment(i).cluster(j).spikes(spike) * fs) + (1:51);
            snippets(spike,:,:) = lfptemp(spikewin, :);
        end
        msnip = squeeze(mean(snippets, 1));
        peakamp = range(msnip, 1);
        [~, chord] = sort(peakamp, 'descend');
        peakChannel2(i,j,:) = chord;
        experiment(i).cluster(j).peakChannel2 = squeeze(peakChannel2(i,j,:));
        experiment(i).cluster(j).msnip = msnip';
        experiment(i).cluster(j).snips = snippets;
    end
end

