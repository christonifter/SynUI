totalneurons = 4;
x = [ones(2,1).*200; ones(2,1).*600];
y = [(1:2).*400 (1:2).*400-100]';
z = ones(4,1).*450;

%%now they are 1000x1 vectors. put next to each other 
neuron_positions = [x,y,z];
electrode_z = [ones(32,1).*450];
electrode_x = [ones(16,1).*250;ones(16,1).*650];
electrode_y = repmat([0:100:1500]',2,1);
electrode_matrix = [electrode_x, electrode_y, electrode_z];

distance = NaN(32,totalneurons);
for electrode = 1:32
    for neuron = 1:totalneurons
        position_difference = (neuron_positions(neuron,:)- electrode_matrix(electrode,:));
        distance(electrode,neuron) = norm(position_difference,2);
    end
end
amplitude = 400./distance;


spike_train
halfnoise = normrnd(0,.2, [size(amplitude, 1) size(spiketrain, 1)/2]);
simtanks = amplitude * spiketrain' + [halfnoise; halfnoise];
figure(13); plot(simtanks(:,1:(samplingrate*2))'-(1:32).*2, 'k')

ops.fbinary = ['D:\Spikes\simIC\twentyneurons\mergedata-001-i16.dat'];
fid = fopen(ops.fbinary, 'w');
x = max(max(abs(simtanks)));
fwrite(fid, int16((2^15)./x .* simtanks), 'int16');
fclose(fid);
rez = tanksort(ops, 'D:\Spikes\simIC\twentyneurons\');
<<<<<<< Updated upstream

spikeTimes     = rez.st3(:,1);
spikeClusters = rez.st3(:,2);
=======
[~,sorti] = sort(rez.st3(:,1));
spikeTimes = rez.st3(sorti,1);
spikeClusters = rez.st3(sorti,2);
spikeAmp = rez.st3(sorti,4);
>>>>>>> Stashed changes
peakChannel = rez.iNeighPC';
peakChannel2 = NaN(numel(unique(spikeClusters)),32);

<<<<<<< Updated upstream
% [sr,sc] = find(triu(rez.simScore,1)>0.9)
%     spikeClusters2= spikeClusters
% for i = 1:numel(sr)
%     spikeCluster2(spikeClusters == sc) = sr;
% end
% 
% 
    clear cluster
    for j = 1:max(spikeClusters)
=======
channelnoise = std(simtanks, [], 2);
clear cluster
for j = 1:max(spikeClusters)
>>>>>>> Stashed changes
        clind = find(spikeClusters == j);
        cluster(j).spikes = (spikeTimes(clind))./rez.ops.fs;
        cluster(j).peakChannel = peakChannel(j,:);
    end

<<<<<<< Updated upstream
=======
% mergy;
spikeTimes = rez.st3(:,1);
spikeClusters = rez.st3(:,2);


clear cluster
for j = 1:max(spikeClusters)
    clind = find(spikeClusters == j);
    cluster(j).spikes = (spikeTimes(clind))./rez.ops.fs;
    cluster(j).peakChannel = peakChannel(j,:);
end
for j = 1:max(spikeClusters)
    snippets = NaN(numel(cluster(j).spikes), 51 ,32); %snippets of spikes from this experiment, all clusters
    for spike = 1:numel(cluster(j).spikes)
        spikewin = round(cluster(j).spikes(spike) * rez.ops.fs) + (1:51);
        snippets(spike,:,:) = simtanks(:, spikewin)';
    end
    msnip = squeeze(mean(snippets, 1));
    peakamp = range(msnip, 1);
    [~, chord] = sort(peakamp, 'descend');
    peakChannel2(j,:) = chord;
    cluster(j).peakChannel2 = squeeze(peakChannel2(j,:));
    cluster(j).msnip = msnip';
    cluster(j).snips = snippets;
end


>>>>>>> Stashed changes
size(spikeTimes)
max(spikeClusters)
figure(14); subplot(2,1,1); imagesc(rez.U(:,:,1)); 
xlabel('Cluster'); ylabel('Channel');colorbar;
subplot(2,1,2); imagesc(triu(rez.simScore,1)); colormap(jet); colorbar;
xlabel('Cluster'); ylabel('Cluster'); 
figure(12);subplot(1,2,2); plot(rez.W(1:51,:,1)-(1:max(spikeClusters)))

