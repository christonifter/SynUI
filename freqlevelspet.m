function [spetfreq, spetlevel, stimspets, trials] = freqlevelspet(spets, stimons, frqs, lvls)
frqlist = sort(unique(frqs));
lvllist = sort(unique(lvls));
trials = NaN(numel(lvllist), numel(frqlist));
for frq = 1:numel(frqlist)
    for lvl = 1:numel(lvllist)
        trials(lvl, frq) = sum(frqs == frqlist(frq) & lvls == lvllist(lvl));
    end
end

if isempty(spets)
    spettable = [];
else
    spettable = repmat(spets, 1, numel(stimons)) - stimons';
end
spettable(spettable<0) = NaN;

[stimspets, spetlaststimi] = min(spettable, [], 2, 'omitnan');
spetfreq = frqs(spetlaststimi);
spetlevel = lvls(spetlaststimi);
