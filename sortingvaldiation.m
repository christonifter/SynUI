inpath = 'D:\Tanks\';
outpath = 'D:\Spikes\';
subject = '19-66-190802-115555\';
site = 'Site02Run01\';
groupcode = 2;
[groups,folders,~]=xlsread([inpath subject 'Experiments.xlsx']);
experimentsi = find(groups == groupcode);

%% load data (slow)
load([outpath subject site 'rez.mat'])
load([outpath subject site 'peakChannel2.mat'])
spikeTimes     = rez.st3(:,1);
spikeClusters = rez.st3(:,2);
peakChannel = rez.yc(rez.iNeighPC)';
mergedata = [];
for i = 1:numel(experimentsi)
    tankpath = [inpath subject char(folders(experimentsi(i)))];
    data = readtank(tankpath);
    mergedata = [mergedata; data.LFP];
    segmentends(i) = size(mergedata,2);
end
mergedata = [zeros(26, 32); mergedata; zeros(26, 32)];
fs = data.fs;
clear data
segmentstarts = 1+ [0; segmentends(1:(end-1))];

%% save snippets (kinda slow)
for j = 1:max(spikeClusters)
    j
    clind = find(spikeClusters == j);
    spets = (spikeTimes(clind))./fs;
    snippets = NaN(numel(clind), 51 ,32); %snippets of spikes from this experiment, all clusters
    for spike = 1:numel(clind)
        spikewin = spikeTimes(clind(spike)) + (1:51);
        snippets(spike,:,:) = mergedata(spikewin, :);
    end
    msnip(j,:,:) = squeeze(mean(snippets, 1));
    peakamp(j,:) = squeeze(range(msnip(j,:,:), 2));
end

%% count spikes per experiment (fast)
spikecount = NaN(max(spikeClusters), numel(experimentsi));
for i = 1:numel(experimentsi)
    load([outpath subject site char(folders(experimentsi(i))) '.mat'])
    for j = 1:max(spikeClusters)
        spikecount(j,i) = numel(cluster(j).spikes);
    end
end

%% plots (fast)
figure(1); clf(1); hold on;
for j = 1:max(spikeClusters)
    pj = max(peakamp(j,:), [], 2);
    plot(repmat(((-25:25)./60)', [1, 32]) + (1:32), squeeze(msnip(j,:,:))./pj - j, 'k')
    [~, chord] = sort(peakamp(j,:), 'descend');
    peakChannelmerge(j,:) = chord;
    plot((-25:25)./60 + peakChannelmerge(j, 1), squeeze(msnip(j,:,peakChannelmerge(j, 1)))./pj - j, 'r')
end
hold off;

nearchan = peakChannel2(:,:,1)';

figure(2); clf; imagesc(nearchan); colormap('jet');
xlabel('Experiment')
ylabel('Cluster')
colorbar()
title('Nearest Channel')
figure(3); clf; subplot(2,2,1)
imagesc(log(spikecount)); colormap('jet');
xlabel('Experiment')
ylabel('Cluster')
title('Spike Count')
subplot(2,2,2)
plot(sum(spikecount,2), -1.*(1:max(spikeClusters)))
subplot(2,2,3)
plot(sum(spikecount,1))

nearchan(spikecount<100) = NaN;
figure(4); imagesc(nearchan); colormap('jet')
hosatchel = get(gca, 'colormap');
hosatchel(1,:) = [1 1 1];
set(gca, 'colormap', hosatchel);
nanstd(nearchan, [], 2)



% 
% 
% 
% 
% spikewin = -20:20;
% spksamp = cluster(11).spikes.*data.fs;
% spikewins = round(spksamp(spksamp>20 & spksamp<(length(data.LFP)-20)) + spikewin);
% 
% snippets = double(reshape(data.LFP(spikewins, 10), size(spikewins)));
% r = triu(corr(snippets'));
% n = size(snippets, 1);
% size(snippets)
% corscore = sum(sum(triu(r)))./(n*(n-1)/2);
% corscore
% plot(spikewin, snippets)
