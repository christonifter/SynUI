
function synftc(stimspets, spetfreq, analwin, channels, trials, chanlist, frqs, ax)
if ~exist('ax', 'var')
    ax = axes();
end
nchans = numel(chanlist);
frqlist = sort(unique(frqs));
ftc = NaN(nchans, numel(frqlist));
scale = 1./(analwin(2)-analwin(1))./trials;
for chan = 1:nchans
    ftc(chan,:) = hist(spetfreq(stimspets>analwin(1)& stimspets<analwin(2) & ...
        channels ==chanlist(chan)), frqlist);
end
ftc = ftc .*scale;
imagesc(ax, ftc);
colormap(ax,'jet')
set(ax,'XTick', 1:numel(frqlist))
set(ax,'YTick', 1:nchans)
set(ax,'YTickLabel', chanlist)
ylabel(ax, 'Channel')

colorbar(ax)
axis(ax, [.5, numel(frqlist)+.5, .5, nchans+.5 ])
hosatchel = get(ax, 'colormap');
hosatchel(1,:) = [1 1 1];
set(ax, 'colormap', hosatchel);