function sttracesplot(app, ax, traces, fs)
if ~isempty(traces)
y = traces(:,:,1);
t = (1:size(y, 2))./fs;
cla(ax);
hold(ax, 'on');
plot(ax, t, y', 'k-');
plot(ax, t, mean(y)', 'r-');
hold(ax, 'off');
end