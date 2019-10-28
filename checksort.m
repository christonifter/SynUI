% function checksort(folderpath)
uselist = SUlist;
% uselist = 1:numel(cluster);
sttemplate = NaN(32,51,numel(uselist));

for j = 1:numel(uselist)
    sttemplate(:,:,j) = permute(cluster(uselist(j)).msnip, [1 3 2]);
end

nrows = min([10 numel(uselist)]);
ncols = 5;
figure(1); clf;

% nrows = ceil(sqrt(numel(uselist)));
for j = 1:nrows
    subplot(nrows, ncols, ncols*(j-1)+1)
    imagesc(sttemplate(:,:,j), [-2E-4 2E-4]); colormap('jet')
end

subplot(1,ncols,2)
plot(0,0);
hold on;
for j = 1:nrows
    plot(1E3*squeeze(cluster(uselist(j)).snips(1:10:end,:))'-j, 'k')
    plot(1E3*sttemplate(peakChannel2(1,uselist(j),1),:,j)'-j, 'r')
end
ylim([-10.5 -0.5])
subplot(1,ncols,3)
plot(0,0);
hold on;
for j = 1:nrows
    st = rez.st3(rez.st3(:,2) == uselist(j), 1)./rez.ops.fs;
    [y, x] = hist(diff(sort(st)), 0:2E-4:1E-2);
    plot(x,y./max(y(1:end-1))-j)
end
grid on;
xlim([0 9E-3])

subplot(1,ncols,4)
plot(0,0);
hold on;
for j = 1:nrows
    plot(rez.acg(401:600,uselist(j))./max(rez.acg(401:600,uselist(j)))-j)
end

st = rez.st3(rez.st3(:,2) == 13, 1)./rez.ops.fs;
[K, Qi, Q00, Q01, Ri] = ccg(st, st, 500, 2E-4);