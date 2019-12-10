figure(7);    
c = 1;
rez.merges2=[];
for i = 1:size(rez.merges, 1)
    clus1 = rez.merges(i,1);
    clus2 = rez.merges(i,2);
    gapmat1 = cluster(clus1).msnip - cluster(clus2).msnip - 2* channelnoise;
    gapmat2 = cluster(clus2).msnip - cluster(clus1).msnip - 2* channelnoise;
    [r1, c1] = find(gapmat1>0);
    [r2, c2] = find(gapmat2>0);
%     subplot(4,5,c)
%     imagesc((cluster(clus1).msnip - cluster(clus2).msnip)'./channelnoise', [-3 3]); colorbar; colormap('jet');
%     title([num2str(clus1) '-' num2str(clus2)])
%     hold on;
%     plot([r1; r2], [c1; c2], 'ko');
%     hold off;
    if numel([r1; r2]) < 3
        rez.merges2 = [rez.merges2; rez.merges(i,:)];
    else
%         if c > 19
%             c = 1; figure();
%         else
%             c = c+1;
%         end
    end
end
out = supergroups(rez.merges2);
%merge clusters and remove empty clusters
if numel(out)>1
    for clump = unique(out(:,2))'
        userows = find(out(:,2) == clump);
        nsp = zeros(numel(userows), 1);
        for i = 1:numel(userows)
            nsp(i) = sum(spikeClusters == out(userows(i), 1));
        end
        [~,bestclusti] = max(nsp);
        out(userows, 2) = out(userows(bestclusti), 1);
    end
    otherclusters = setdiff(unique(spikeClusters), out(:,1));
    clustlist = sortrows([out; [otherclusters(:) otherclusters(:)]]);
    spikeClusters2 = spikeClusters;
    for i = 1:size(out,1)
        spikeClusters2(spikeClusters==out(i,1)) = out(i,2);
    end
    newclusters = unique(spikeClusters2);
    ignitionremix = [newclusters, (1:numel(newclusters))'];
    peakChannel1 = rez.iNeighPC';
    for i = 1:numel(unique(spikeClusters2))
        spikeClusters3(spikeClusters2==newclusters(i), 1) = i;
        clustlist(find(clustlist(:,2)==newclusters(i)), 3) = i;
        merged2orig(i,1) = clustlist(find(clustlist(:,3) == i, 1), 2);
    end
        peakChannel = peakChannel1(merged2orig,:);
    rez.st3(:,2) = spikeClusters3;
else
    peakChannel = rez.iNeighPC';
    clustlist = sort(unique(spikeClusters));
    merged2orig = sort(unique(spikeClusters));
end