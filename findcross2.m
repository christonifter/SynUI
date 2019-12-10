function [ups, downs] = findcross2(y)
t1=y(1:end-1);
t2=y(2:end);
tt=t1.*t2;
indx=find(tt<0);

zos = find(y(2:end-1) == 0);
zz = y(zos).*y(zos+2);
indx2 = zos(find(zz<0))+1;
indx = sort([indx;indx2]);

if ~isempty(indx)
    if y(indx(1)) < 0
        ups = indx(1:2:end);
        downs = indx(2:2:end);
    elseif y(indx(1)) == 0
        if y(indx(1)-1) < 0
            ups = indx(1:2:end);
            downs = indx(2:2:end);
        else
            ups = indx(2:2:end);
            downs = indx(1:2:end);
        end
    else
        ups = indx(2:2:end);
        downs = indx(1:2:end);
    end
else
    ups = [];
    downs = [];
end
end
