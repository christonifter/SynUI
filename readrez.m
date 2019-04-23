% st = readNPY('spike_times.npy'); % these are in samples, not seconds
% clu = readNPY('spike_clusters.npy');

figure(1); clf(1);
% plot spikeWaves
for i = 1:20
    j = i +20;
    subplot(4,5,i)
    clind = find(out.spikeTemplates == j);
    channel = out.peakChannel(j, 1);
    plot((-50:50).*1000./rez.ops.fs, squeeze(out.spikeWaves(channel,:, clind)) - squeeze(mean(out.spikeWaves(channel,1, clind), 2))');
    title(num2str(out.peakChannel(j)))
end


figure(2); clf(2); ax = gca;    hold(ax, 'on');

%plot rasters
for i = 1:max(out.spikeTemplates)
    clind = find(out.spikeTemplates == i);
    plot(ax, out.spikeTimes(clind)./rez.ops.fs, i, 'k*')
    xlim(ax, [0 60])
end
hold(ax, 'off');