function [spetfreq, stimspets, trials] = freqanal(spets, stimons, frqs)
frqlist = sort(unique(frqs));
trials = NaN(size(frqlist));
for frq = 1:numel(frqlist)
    trials(frq) = sum(frqs == frqlist(frq));
end

if isempty(spets)
    spettable = [];
else
    spettable = repmat(spets, 1, numel(stimons)) - stimons';
end
spettable(spettable<0) = NaN;

[stimspets, spetlaststimi] = min(spettable, [], 2, 'omitnan');
spetfreq = frqs(spetlaststimi);

