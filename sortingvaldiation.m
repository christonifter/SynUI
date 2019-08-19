
spikewin = -20:20;
spksamp = cluster(1).spikes.*data.fs;
spikewins = round(spksamp(spksamp>20 & spksamp<(length(data.LFP)-20)) + spikewin);

snippets = double(reshape(data.LFP(spikewins, 1), size(spikewins)));
r = triu(corr(snippets'));
n = size(snippets, 1);
corscore = sum(sum(triu(r)))./(n*(n-1)/2);
corscore
plot(spikewin, snippets)
