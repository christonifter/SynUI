subplot(2,1,1)
t = (1:size(data.raw, 1))./data.fs;
chan = 16;
figure(1)
plot(t, data.raw(:,chan))
xlim([0 1])
subplot(2,1,2)
plot(t, data.LFP(:,chan))
xlim([0 1])
period = 47.6/1000;
analwin = [0 20];
cyclet = analwin(1):period:analwin(2);
cycletsamp = round(cyclet*data.fs);
cycleperiod = min(diff(cycletsamp));
for cycle = 1:numel(cyclet)
    startsamp = 1+cycletsamp(cycle);
    endsamp = cycletsamp(cycle)+cycleperiod;
    lfpcycle(:,:,cycle) = data.LFP(startsamp:endsamp,:);
end
figure(2)
plot((1:cycleperiod)./data.fs, mean(lfpcycle, 3));