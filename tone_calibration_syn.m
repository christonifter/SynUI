data=TDTbin2mat('D:\Synapse\Tanks\19-19-200211-125706\Block3\');
micsensitivity = 2; %in mV/Pa
micgain = 40; %in dB

freqvec = data.epocs.Fre1.data;
stimtimes = data.epocs.Fre1.onset + 0.1;
leveldata = data.streams.Mic1.data;
levelfs = data.streams.Mic1.fs;
leveltimes = round(stimtimes.*levelfs)+1;
levelvec = NaN(size(freqvec));

stimdata = data.streams.Stim.data;
stimfs = data.streams.Stim.fs;

for i = 1:numel(freqvec)
    levelvec(i) = mean(leveldata(leveltimes(i, 1):leveltimes(i, 1) + round(0.8.*levelfs)));
end
% levelsplvec = 20*log10(levelvec/1000/100/0.001459)+93.9;
levelsplvec = 20*log10(levelvec/micsensitivity)+93.9 - micgain;


figure(1)
subplot(2,2,1)
plot((1:numel(stimdata))./stimfs, stimdata)
title('DAC signal')
subplot(2,2,3)
plot((1:numel(leveldata))./levelfs, leveldata)
title('Microphone recorded RMS level')
ylabel('mV')
xlabel('sec')

subplot(2,2,2)
plot(freqvec./1000, levelvec, 'o-')
title('amplitude response')
ylabel('mV')
ylim([0 300])
subplot(2,2,4)
plot(freqvec./1000, levelsplvec, 'o-')
title('amplitude response')
ylabel('dB SPL')
xlabel('kHz')
ylim([80 100])
