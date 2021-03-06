function [peaks, lats, slopes] =lfppeaks(y)
[peaks(1,:), lats(1,:)] = min(y);
[slopes(1,:), lats(4,:)] = min(diff(y)); %instantaneous falling slope
[slopes(2,:), lats(5,:)] = max(diff(y)); %instantaneous rising slope

for chan = 1:size(y,2)
    [peaks(2,chan), lats(2,chan)] = max(y(1:lats(1,chan), chan));
    [peaks(3,chan), i] = max(y(lats(1,chan):end, chan));
    lats(3,chan) = i + lats(1,chan);
     [u,d] = findcross(y(:,chan));
     x = find(d<lats(1,chan), 1, 'last');
     if ~isempty(x)
         lats(6,chan) = d(x);
     else
         lats(6,chan) = NaN;
     end
     x = find(u>lats(1,chan), 1, 'first');
     if ~isempty(x)
         lats(7,chan) = u(x);
     else
         lats(7,chan) = NaN;
     end
end
    slopes(3,:) = peaks(1,:)./(lats(1,:)-lats(6,:)); %average falling slope
    slopes(4,:) = -1.*peaks(1,:)./(lats(7,:)-lats(1,:)); %average rising slope