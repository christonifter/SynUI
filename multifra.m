function fra = multifra(stimspets, spetfreq, spetlevel, analwin, channels, trials, chanlist, frqs, lvls, ax)
if ~exist('ax', 'var')
    ax = axes();
end
frqlist = sort(unique(frqs));
lvllist = sort(unique(lvls));
fra = NaN(numel(lvllist), numel(frqlist), numel(chanlist));
scale = 1./(analwin(2)-analwin(1))./trials./numel(chanlist);
for chanc = 1:numel(chanlist)
    chani = (channels == chanlist(chanc));
    for frqc = 1:numel(frqlist)
        frqi = (spetfreq == frqlist(frqc));
        for lvlc = 1:numel(lvllist)
            lvli = (spetlevel == lvllist(lvlc));
            fra(lvlc,frqc, chanc) = sum(stimspets>analwin(1)& stimspets<analwin(2) & ...
                chani & frqi & lvli);
        end
    end
end
fra = fra.*scale;
mfra = sum(fra, 3, 'omitnan');
imagesc(ax, flipud(mfra));
colormap(ax,'jet')
set(ax,'XTick', 1:numel(frqlist))
set(ax,'XTickLabel', frqlist./1000)
set(ax,'YTick', 1:numel(lvllist))
set(ax,'YTickLabel', lvllist(numel(lvllist):-1:1))
ylabel(ax, 'Level (dBSPL)')
xlabel(ax, 'Frequency (kHz)')
title(ax, 'Frequency Response Area')
colorbar(ax)
axis(ax, [.5, numel(frqlist)+.5, .5, numel(lvllist)+.5 ])

% v1 = [zeros(1,30) linspace(1/20, 1, 20) ones(1,20) linspace(1,1/20, 20) zeros(1,30)]
hosatchel = get(ax, 'colormap');
hosatchel(1,:) = [1 1 1];
set(ax, 'colormap', hosatchel);
