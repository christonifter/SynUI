for j = 1:numel(clusterallexp)
    nsp(j,1) = clusterallexp(j).nsp;
    amp(j,1) = mean(rez.st3(rez.st3(:,2) == j, 4));
    slope(j,1) = mean(clusterallexp(j).slopes);
    snr(j,1) = max(clusterallexp(j).snrs(:));
    snrsnr(j,1) = snr(j,1)./mean(max(clusterallexp(j).snrs));
    bimodality(j,1) = HartigansDipTest(rez.st3(rez.st3(:,2) == j, 4));
    asym(j,1) = corr(clusterallexp(j).msnip((peaklat-halfwin):peaklat, peakChannel(j)), clusterallexp(j).msnip((peaklat+halfwin):-1:peaklat, peakChannel(j)));
    for k = (j+1):numel(clusterallexp)
        stsim(j,k) = corr(clusterallexp(j).msnip(:), clusterallexp(k).msnip(:));
    end    
end
    refractory = rez.est_contam_rate(:);
    position = peakChannel(:,1);
    
    
sortingquality = table(nsp, snr, snrsnr, position, refractory, bimodality, asym);

figure(10); clf; ax = gca;
subplot(2,2,1)
plot(sortingquality.nsp, sortingquality.snr, 'k.');
xlabel('nsp')
ylabel('snr')
set(gca, 'xScale', 'log')
subplot(2,2,2)
plot(sortingquality.snrsnr, sortingquality.snr, 'k.');
xlabel('snrsnr')
ylabel('snr')
subplot(2,2,3)
plot(sortingquality.refractory, sortingquality.snr, 'k.');
xlabel('refractory')
ylabel('snr')
subplot(2,2,4)
plot(sortingquality.asym, sortingquality.snr, 'k.');
xlabel('asym')
ylabel('snr')

figure(11);
imagesc(stsim, [0 1]); colormap('jet'); colorbar;
figure(12); clf; ax = gca;
scatter3(sortingquality.snr, sortingquality.snrsnr, sortingquality.asym, 1, 'ko');
hold on;
uselist = find(snr > 2 & snrsnr > 5 & asym > 0.2);
scatter3(sortingquality.snr(uselist), sortingquality.snrsnr(uselist), sortingquality.asym(uselist), 2, 'ro');
xlabel('snr')
ylabel('snrsnr')
zlabel('asym')
hold off;
