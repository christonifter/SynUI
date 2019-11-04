figure(1)
% function checksort(folderpath)
for j = 1:numel(clusterallexp)
    clspikecount(j,1) = clusterallexp(j).nsp;
%     clamp(j,1) = range(clusterallexp(j).msnip(:);
    clamp(j,1) = mean(rez.st3(rez.st3(:,2) == j, 4));
    clslope(j,1) = mean(clusterallexp(j).slopes);
    sampnoise = std(clusterallexp(j).snips(:));
    sampsnr = abs(clusterallexp(j).msnip(:,peakChannel(j,1)))./sampnoise;
    clsnr(j,1) = max(sampsnr);

end

% uselist = SUlist;
[~, uselist] = sort(clsnr, 'descend','MissingPlacement','last');
sttemplate = NaN(32,51,numel(uselist));

for j = 1:numel(uselist)
    sttemplate(:,:,j) = permute(clusterallexp(uselist(j)).msnip', [1 3 2]);
end

nrows = min([10 numel(uselist)]);
ncols = 6;
figure(1); clf; 
subplot(1,ncols,1); plot(0,0); hold on;
for j = 1:nrows
    imagesc(0, (j-1)*33, sttemplate(:,:,j)./15E-5, [-1 1]); colormap('jet')
    imagesc(52, (j-1)*33, permute(rez.dWU(6:56,:,uselist(j)), [2 1 3])./15, [-1 1]); colormap('jet')
    plot([0 52*2], [peakChannel(uselist(j),1) peakChannel(uselist(j),1)]-1+(j-1)*33, 'k:')

end
axis([0 103 0 33*j])

subplot(1,ncols,2)
plot(0,0);
hold on;
for j = 1:nrows
     plot(1E3*squeeze(clusterallexp(uselist(j)).snips(1:ceil((1+clusterallexp(uselist(j)).nsp)/200):end,:))'+j, 'k')
    plot(1E3*sttemplate(peakChannel(uselist(j),1),:,j)'+j, 'r')

    sampnoise = std(clusterallexp(uselist(j)).snips(:));
    sampsnr = abs(clusterallexp(uselist(j)).msnip(:,peakChannel(uselist(j),1)))./sampnoise;
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
    plot(rez.acg(401:600,uselist(j))./max(rez.acg(401:600,uselist(j)))+j)
end

subplot(1,ncols,5)
plot(0,0);
hold on;
for j = 1:nrows
    plot(rez.st3(rez.st3(:,2) == uselist(j), 1), rez.st3(rez.st3(:,2) == uselist(j), 4)./200+j, 'k.')
end
for k = 2:(numel(segmentstarts))
    plot([segmentstarts(k), segmentstarts(k)], [0 nrows+1], 'r')
end
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
c = 1;
for clust = 1:nrows
    for exp = 1:numel(experiment)
        imagesc((exp-1)*52, (clust-1)*33, experiment(exp).cluster(uselist(clust)).msnip', [-3E-4 3E-4]); colormap('jet'); colorbar;
    end
    plot([0 52*exp], [peakChannel(uselist(clust),1) peakChannel(uselist(clust),1)]-1+(clust-1)*33, 'k:')
end
axis([0 52*exp 0 clust*33])
% st1 = rez.st3(rez.st3(:,2) == 18, 1)./rez.ops.fs;
% 
% [K, Qi, Q00, Q01, Ri] = ccg(st1, st1, 500, 1E-3);
