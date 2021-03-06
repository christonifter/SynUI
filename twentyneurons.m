    totalneurons = 20;
    x = [ones(10,1).*200; ones(10,1).*600];
    y = [(1:10).*160 (1:10).*160-100]';
    z = ones(20,1).*450;

% distance test
%     totalneurons = 16;
%     x = [ones(16,1).*260];
%     y = (1:16)'.*100;
%     z = (1:16)'.*50;

%spatial ambiguity test
%     totalneurons = 16;
%     x = [ones(16,1).*600];
%     y = repmat([(1:2:8).*200 (2:2:8).*200], 1, 2)';
%     z = reshape(repmat([(1:4).*100+200], 4, 1), 16, 1);
    
% totalneurons = 4;
% x = [ones(2,1).*200; ones(2,1).*600];
% y = [(1:2).*400+400 (1:2).*400+300]';
% z = ones(4,1).*450;

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
amplitude = 400./distance; %spike amplitude = 1 at distance = 100 um. noise amplitude = 1.


spike_train
rng(1)
simtanks = [zeros(32, 26) amplitude * spiketrain' + normrnd(0,.2, [size(amplitude, 1) size(spiketrain, 1)]) zeros(32, 26)];
figure(13); plot(0.5*simtanks(:,1:(samplingrate*10))'-(1:32), 'k')
ops.fbinary = ['D:\Spikes\twentyneurons\run01\mergedata-001-i16.dat'];
fid = fopen(ops.fbinary, 'w');
x = max(max(abs(simtanks)));
fwrite(fid, int16((2^15)./x .* simtanks(:, 27:(end-26))), 'int16');
fclose(fid);
rez = tanksort(ops, 'D:\Spikes\twentyneurons\run01\');
spikeTimes     = rez.st3(:,1);
spikeClusters = rez.st3(:,2);
peakChannel = rez.iNeighPC';
peakChannel2 = NaN(numel(unique(spikeClusters)),32);
% 
% clear cluster
% for j = 1:max(spikeClusters)
%     clind = find(spikeClusters == j);
%     cluster(j).spikes = (spikeTimes(clind))./rez.ops.fs;
%     cluster(j).peakChannel = peakChannel(j,:);
% end
% 
% figure(14); subplot(2,1,1); imagesc(rez.U(:,:,1)); 
% xlabel('Cluster'); ylabel('Channel');colorbar;
% subplot(2,1,2); imagesc(triu(rez.simScore,1)); colormap(jet); colorbar;
% xlabel('Cluster'); ylabel('Cluster'); 
% figure(12);subplot(1,2,2); plot(rez.W(1:51,:,1)-(1:max(spikeClusters)))
% channelnoise = std(simtanks, [], 2);
% for j = 1:max(spikeClusters)
%     snippets = NaN(numel(cluster(j).spikes), 51 ,32); %snippets of spikes from this experiment, all clusters
%     for spike = 1:numel(cluster(j).spikes)
%         spikewin = round(cluster(j).spikes(spike) * rez.ops.fs) + (1:51);
%         snippets(spike,:,:) = simtanks(:, spikewin)';
%     end
%     msnip = squeeze(mean(snippets, 1));
%     cluster(j).msnip = msnip';
% end
% 
% figure(21);    
% c = 1;
% rez.merges2=[];
% for i = 1:size(rez.merges, 1)
%     clus1 = rez.merges(i,1);
%     clus2 = rez.merges(i,2);
%     gapmat1 = cluster(clus1).msnip - cluster(clus2).msnip - 2* channelnoise;
%     gapmat2 = cluster(clus2).msnip - cluster(clus1).msnip - 2* channelnoise;
%     [r1, c1] = find(gapmat1>0);
%     [r2, c2] = find(gapmat2>0);
%     subplot(4,4,c)
%     imagesc(cluster(clus1).msnip - cluster(clus2).msnip, [-1 1]); colorbar; colormap('jet');
%     title([num2str(clus1) '-' num2str(clus2)])
%     hold on;
%     plot([c1; c2], [r1; r2], 'wx')';
%     hold off;
%     if numel([r1; r2]) == 0
%         rez.merges2 = [rez.merges2; rez.merges(i,:)];
%     else
%         if c > 15
%             c = 1; figure();
%         else
%             c = c+1;
%         end
%     end
% end
% out = supergroups(rez.merges2);
% for i = 1:size(out,1)
%     rez.st3(rez.st3(:,2)==out(i,1), 2) = out(i,2);
% end
% spikeTimes     = rez.st3(:,1);
% spikeClusters = rez.st3(:,2);

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
save('D:\Spikes\twentyneurons\run01\rez.mat', 'rez', '-v7.3');
save('D:\Spikes\twentyneurons\run01\cluster.mat', 'cluster', '-v7.3');