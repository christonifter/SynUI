function [ups, downs] = findcross(y)
t1=y(1:end-1);
t2=y(2:end);
tt=t1.*t2;
indx=find(tt<0);
if ~isempty(indx)
    if y(indx(1)) < 0
        ups = indx(1:2:end);
        downs = indx(2:2:end);
    else
        ups = indx(2:2:end);
        downs = indx(1:2:end);
    end
else
    ups = [];
    downs = [];
end