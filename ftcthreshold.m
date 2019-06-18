function ftcstats = ftcthreshold(ftc, frqlist, spontrate, threshcrit)
nchans = size(ftc, 1);
ftcstats = NaN(nchans, 5);
for chan = 1:nchans
     try
         if numel(frqlist) < threshcrit
            a = fit(frqlist, ftc(chan, :)'-spontrate(chan), 'gauss1');
            if a.a1 > threshcrit %threshold of 10 Hz above spontaneous works okay.
             frqthresh = [-a.c1 * sqrt(log(a.a1/10)) + a.b1, a.c1 * sqrt(log(a.a1/10)) + a.b1] ;
             ftcstats(chan, :) = [a.a1, a.b1, 2*a.c1, frqthresh];
            end
         else
            a = fit(log(frqlist), ftc(chan, :)'-spontrate(chan), 'gauss1');
            ftcstats(chan, :) = [a.a1, exp(a.b1), 2*a.c1/log(2), exp([a.b1 - a.c1, a.b1 + a.c1])];
         end
    catch
    end
end
