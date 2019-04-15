clear x y z
n = 9;
div = 10;
spacing =3;

x = repmat(mod(spacing.*(1:n), div)', 9, 1);

y = NaN(n,9);
for leveli = 1:9
    y(:,leveli) = x((1:n)+leveli);
end
z(:,1) = x;
z(:,2) = reshape(y, numel(y), 1);


n = 5;
div = 7;
spacing =3;
r = repmat(mod(spacing.*(1:n), div)', 9,1);
r(r==5) =4;
r(r==6) =5;
s = [r, x(1:45)];
check = s(:,1)*10+s(:,2);
sort(check);


freqs = 4000.* 2.^(0:0.5:4)';
levels = (10:10:90)';
freqhalfoctvec= freqs(z(:,1));
freqhalfoctvec = [freqhalfoctvec; flip(freqhalfoctvec(1:end-3)); flip(freqhalfoctvec(end-2:end))];
levelhalfoctvec= levels(z(:,2));
levelhalfoctvec = [levelhalfoctvec; flip(levelhalfoctvec(1:end-3)); flip(levelhalfoctvec(end-2:end))];

freqs = 4000.* 2.^(0:4)';

freqoctvec= freqs(s(:,1));
leveloctvec = levels(s(:,2));
freqoctvec = [freqoctvec; flip(freqoctvec(1:end-3)); flip(freqoctvec(end-2:end))];
leveloctvec = [leveloctvec; flip(leveloctvec (1:end-3)); flip(leveloctvec (end-2:end))];
