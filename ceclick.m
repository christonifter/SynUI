fs = 50000;
duration = .15;
t = (1:(duration*fs))./fs;
ncomps = 100;
chirprate = 20;
freqs = 2000+(chirprate:chirprate:(chirprate*ncomps));
comps = zeros(ncomps, numel(t));
kdon = 0.092;
ddon = 0.4356;
phase = (kdon/(1-ddon)) .* freqs .^ -ddon; %phase in sec
phaserad = phase.*freqs.*2.*pi;
for i = 1:ncomps
    comps(i,:) = cos(2*pi*freqs(i)*t + phaserad(i));
end
subplot(2,2,1)
plot(t, sum(comps))
xlabel('Time (sec)')
title(['CE Chirp Time Waveform, ' num2str(chirprate) ' chirps/sec, 2000-4000 Hz'])

[sgram, specf, spect] = pspectrum(sum(comps), fs, 'spectrogram', 'FrequencyLimits', [0 5000]);

subplot(2,2,3)
imagesc(flipud(sgram))
set(gca, 'YTick', 0:100:1000)
set(gca, 'YTickLabel', specf(1001:-100:1))
set(gca, 'XTick', 0:50:numel(spect)-1)
set(gca, 'XTickLabel', spect(1:50:numel(spect)))
colormap('jet')
xlabel('Time (sec)')
ylabel('Freq (Hz)')
title('Spectrogram')

subplot(2,2,2)
plot(t, comps./2+(1:ncomps)')
xlabel('Time (sec)')
ylabel('Component')
title('Component Time Waveforms')

subplot(2,2,4)
plot(phase.*1E3, freqs, 'ko-')
xlabel('Phase Delay (ms)')
ylabel('Frequency (Hz)')
ylim([0 5000])
title('Phase Delay')


