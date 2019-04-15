subplot(3,2,1)
plot(data.raw(:,1))
subplot(3,2,2)
[p,f] = pspectrum(data.raw, data.fs);
[p2,f2] = pspectrum(data.LFPraw, data.LFPfs);
[p3,f3] = pspectrum(data.LFP, data.fs);
semilogy(f, p(:,1)', 'r');
hold on;
semilogy(f, p2(:,1)', 'g');
semilogy(f, p3(:,1)', 'b');
hold off;
xlim([0 1000])

subplot(3,2,3)
plot(data.LFPraw(:,1))

subplot(3,2,5)
plot(data.LFP(:,1))

