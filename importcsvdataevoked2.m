% if ~exist('utbl', 'var')
    statspath = 'C:\Users\cmlem\Documents\ExperimentData\evokeStats\';
    paramspath = 'C:\Users\cmlem\Documents\ExperimentData\evokeParam\';
    outpath = 'C:\Users\cmlem\Documents\ExperimentData\summarystats\';
    y = dir(statspath);
    z = dir(paramspath);
    
    utbl = [];
    partbl = [];
    expseq = [];
        
    for j = 3:numel(y)
        tbl = readtable([statspath '\' y(j).name]);
        partblx = readtable([paramspath '\' z(j).name]);
        bn = char(partblx.BlockName);
        BlockName = {bn};
        yrn = str2double(y(j).name(1:2));
        subn = str2double(y(j).name(4:5));
        try
            subn = str2double(y(j).name(4:6));
        end
        siten = 0;
        exporder = j-2;
        expblock = str2double(bn(1:3));
        expyr = 2000+str2double(bn(end-10:end-9));
        expmo = str2double(bn(end-8:end-7));
        expday = str2double(bn(end-6:end-5));
        exphr = str2double(bn(end-3:end-2));
        expmin = str2double(bn(end-1:end));
        exptime = datetime(expyr, expmo, expday, exphr, expmin, 0);
        tbl.Channel = tbl.Channel + yrn*1E7 + subn *1E4 + siten * 1E3;
        partblx = addvars(partblx, expblock, exptime);
        partbl = [partbl; repmat(partblx, size(tbl, 1), 1)];
        utbl = [utbl; tbl];
    end

    isLDSdriven_discontinuous = utbl.LDSonsetrate_Hz > utbl.Last4secCI95_Hz;
    baselineCI95_Hz = poissinv(95/100, utbl.BaselineRate1_Hz * 2.5)/2.5; %assuming LDS onset rate was measured from a 2.5 sec window   
    isLDSdriven_continuous = utbl.LDSonsetrate_Hz > baselineCI95_Hz;
    utbl = [utbl partbl table(baselineCI95_Hz, isLDSdriven_continuous, isLDSdriven_discontinuous)];
%     utbl = addvars(utbl, isLDSdriven);
% end
writetable(utbl, [outpath 'allchannels.csv'])