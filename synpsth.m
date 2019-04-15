function synpsth(stimspets, binvec, channels, chanlist, ax)
nchans = numel(chanlist);
normbincount = NaN(nchans,numel(binvec));
for chan = chanlist'
    [bincount, bintime] = hist(stimspets(channels==chan), binvec);
    yrange = max(bincount) - min(bincount);
    normbincount(chan,:) = (bincount- min(bincount))./(1.1*yrange) - chan;
%     averate = sum(bincount)/(binvec(end)-binvec(1));
%     annotation(gcf, 'textbox', [.8, 0.98-chan/(1.2*nchans), 0, 0 ] , 'String', num2str(averate/235), 'FitBoxToText', 'on')
end
plot(ax, bintime, normbincount);
grid(ax, 'on');
ylabel(ax, 'Channel')
xlabel(ax, 'Time (sec)')