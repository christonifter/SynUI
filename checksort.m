figure(1)
% function checksort(folderpath)
for j = 1:numel(clusterallexp)
    clspikecount(j,1) = clusterallexp(j).nsp;
%     clamp(j,1) = range(clusterallexp(j).msnip(:);
    clamp(j,1) = mean(rez.st3(rez.st3(:,2) == j, 4));
     clslope(j,1) = mean(clusterallexp(j).slopes);
    snr(j,1) = max(clusterallexp(j).snrs(:));
    snrsnr(j,1) = snr(j,1)./mean(max(clusterallexp(j).snrs));
    bimodality(j,1) = HartigansDipTest(rez.st3(rez.st3(:,2) == j, 4));
    [~,peaklat] = max(clusterallexp(j).snrs(:, peakChannel(j)));
    halfwin = min([peaklat-1, 51-peaklat]);
%     asym(j,1) = sqrt(sum((clusterallexp(j).msnip((peaklat-halfwin):peaklat, peakChannel(j)) - clusterallexp(j).msnip((peaklat+halfwin):-1:peaklat, peakChannel(j))).^2))./snr(j,1);
    asym(j,1) = corr(clusterallexp(j).msnip((peaklat-halfwin):peaklat, peakChannel(j)), clusterallexp(j).msnip((peaklat+halfwin):-1:peaklat, peakChannel(j)));
end
refractory = rez.est_contam_rate(:);
uselist = find(snr > 2 & snrsnr > 5 & asym > 0.2);
% uselist = 1:numel(clusterallexp);
%  uselist = find(rez.good);
% [~, uselist] = sort(asym, 'ascend','MissingPlacement','last');
sttemplate = NaN(32,51,numel(uselist));

for j = 1:numel(uselist)
    sttemplate(:,:,j) = permute(clusterallexp(uselist(j)).msnip', [1 3 2]);
end

nrows = min([10 numel(uselist)]);
ncols = 6;
figure(1); clf; 
subplot(1,ncols,1); plot(0,0); hold on;
for j = 1:nrows
    imagesc(0, (j-1)*33, clusterallexp(uselist(j)).snrs'./3); colormap('jet')
    imagesc(52, (j-1)*33, permute(rez.dWU(6:56,:,merged2orig(uselist(j))), [2 1 3])./15, [-1 1]); colormap('jet')
    plot([0 52*2], [peakChannel(uselist(j),1) peakChannel(uselist(j),1)]-1+(j-1)*33, 'k:')

end
axis([0 103 0 33*j])

subplot(1,ncols,2)
plot(0,0);
hold on;
for j = 1:nrows
     plot(2E3*squeeze(clusterallexp(uselist(j)).snips(1:ceil((1+clusterallexp(uselist(j)).nsp)/200):end,:))'+j, 'k')
%     plot(1E3*sttemplate(clusterallexp(uselist(j)).peakChannel2(1),:,j)'+j, 'r')
    plot(2E3*sttemplate(peakChannel(uselist(j), 1),:,j)'+j, 'r')

    sampnoise = std(clusterallexp(uselist(j)).snips(:));
    sampsnr = abs(clusterallexp(uselist(j)).msnip(:,clusterallexp(uselist(j)).peakChannel2(1)))./sampnoise;
    plot(sampsnr./10+j, 'g');
end
ylim([0.5 nrows+.5])


subplot(1,ncols,3)
plot(0,0);
hold on;
for j = 1:nrows
    st = rez.st3(rez.st3(:,2) == uselist(j), 1)./rez.ops.fs;
    [y, x] = hist(diff(sort(st)), 0:2E-4:1E-2);
    plot(x,y./max(y(1:end-1))+j)
end
set(gca, 'XTick', (0:10)./1E3);
grid on;
xlim([0 9E-3])
ylim([1 nrows+1])

subplot(1,ncols,4)
plot(0,0);
hold on;
for j = 1:nrows
    plot(rez.acg([401:500, 502:600],uselist(j))./max(rez.acg(401:600,uselist(j)))+j)
end
ylim([1 nrows+1])
subplot(1,ncols,5)
plot(0,0);
hold on;
for j = 1:nrows
    plot(rez.st3(rez.st3(:,2) == uselist(j), 1), rez.st3(rez.st3(:,2) == uselist(j), 4)./200+j, 'k.')
end
% for k = 2:(numel(segmentstarts))
%     plot([rez.segmentstarts(k), rez.segmentstarts(k)], [0 nrows+1], 'r')
% end
ylim([1 nrows+1])

subplot(1,ncols,6)
plot(0,0);
hold on;
for j = 1:nrows
    [y, x] = hist(rez.st3(rez.st3(:,2) == uselist(j), 4), 50);
    plot(x,y./max(y)+j)
end
ylim([1 nrows+1])

figure(2); clf; hold on;
[FILEPATH,~,~] = fileparts(rez.ops.fbinary);
exps = dir([FILEPATH, '\*.mat']);
for exp = 1:numel(exps)-2
    load([FILEPATH '\' exps(exp).name]);
    for clust = 1:nrows
        imagesc((exp-1)*52, (clust-1)*33, squeeze(cluster(uselist(clust)).msnip), [-1E-4 1E-4]); colormap('jet'); colorbar;
    end
end
for clust = 1:nrows
    plot([0 52*(numel(exps)-2)], [peakChannel(uselist(clust),1) peakChannel(uselist(clust),1)]-1+(clust-1)*33, 'k:')
end
axis([0 52*exp 0 clust*33])
% 
% figure(3); clf; ax = gca; hold on;
%     vvec = [1:-1/16:0 0:1/16:1 1:-1/16:0 0:1/16:1]';
%     vv = 1:32;
%     colmat = [vvec(vv+32) vvec(vv+22) vvec(vv+12)];
%     
%     
% for clust = 1:10
%     spikei = rez.st3(:,2) == uselist(clust);
%     spikePC = rez.cProjPC(spikei,:,peakChannel(clust,1));
%     scatter3(spikePC(:,1), spikePC(:,2), spikePC(:,3), .1, 'MarkerEdgeColor', colmat(mod(clust, 16).*2, :));
% end
% grid on;
% for i = 1:600
%     view(ax, i,30);
%     pause(.1)
% end
