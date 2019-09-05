function fra = multifra(stimspets, spetfreq, spetlevel, analwin, channels, trials, chanlist, frqs, lvls, ax)
if ~exist('ax', 'var')
    ax = axes();
end
frqlist = sort(unique(frqs));
lvllist = sort(unique(lvls));
fra = NaN(numel(lvllist), numel(frqlist), numel(chanlist));
scale = 1./(analwin(2)-analwin(1))./trials;
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
fra = permute(fra.*scale, [3 1 2]);
mfra = squeeze(sum(fra, 1, 'omitnan'))./numel(chanlist);

cscale = [0 max(mfra(:))];

imagesc(ax, flipud(mfra));
colormap(ax,'jet')
set(ax,'XTick', 1:numel(frqlist))
set(ax,'XTickLabel', round(frqlist./100)./10)
set(ax,'YTick', 1:numel(lvllist))
set(ax,'YTickLabel', lvllist(numel(lvllist):-1:1))
ylabel(ax, 'Level (dBSPL)')
xlabel(ax, 'Frequency (kHz)')
title(ax, 'Frequency Response Area')
colorbar(ax)
axis(ax, [.5, numel(frqlist)+.5, .5, numel(lvllist)+.5 ])


% ngrads = 100;
% vvec = [zeros(1, ngrads) 0:1/ngrads:1 ones(1, ngrads) 1:-1/ngrads:0 zeros(1, ngrads)]';
% vv = 1:(ngrads*3.5);
% colmat = [vvec(vv) vvec(vv+round(ngrads*.75)) vvec(vv+round(ngrads*1.5))];
% threshi = round(ngrads*thresh/max(mfra(:)));
% 
% colmat(1,:) = [1 1 1];
%         
% hosatchel = get(ax, 'colormap');
% hosatchel(1,:) = [1 1 1];
% set(ax, 'colormap', colmat);
% % set(ax, 'colormap', hosatchel);
