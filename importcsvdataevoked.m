% if ~exist('utbl', 'var')
    statspath = 'C:\Users\cmlem\Documents\ExperimentData\evokeStats\';
    paramspath = 'C:\Users\cmlem\Documents\ExperimentData\evokeParam\';
    LDSpath = 'C:\Users\cmlem\Documents\ExperimentData\evokeLDS\';
    outpath = 'C:\Users\cmlem\Documents\ExperimentData\summarystats\';
    y = dir(statspath);
    z = dir(paramspath);
    v = dir(LDSpath);
    
    utbl = [];
    partbl = [];
    expseq = [];
        
    for j = 3:numel(y)
        tbl = readtable([statspath '\' y(j).name]);
        partblx = readtable([paramspath '\' z(j).name]);
        psthtbl = readtable([LDSpath '\' v(j).name]);
        bn = char(partblx.BlockName);
        BlockName = {bn};
        yrn = str2double(y(j).name(1:2));
        subn = str2double(y(j).name(4:5));
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
        onsetrate = table2array(psthtbl(1,2:size(psthtbl,2)))';
%         partblx = partblx(:, 1:11);
        if ~isempty(cell2mat(regexp(partblx.Properties.VariableNames, 'LDSPreMS')))
            Time = partblx.Time;
            LDSPreSEC = partblx.LDSPreMS./1000;
            LDSDurSEC = partblx.LDSDurationMS./1000;
            LDSISISec = partblx.LDSISIMS./1000;
            LDSNPulses = partblx.LDSRepeats;
            if ~isempty(regexpi(y(j).name, 'dB'))
                LDSLevel = str2double(y(j).name((1:2) + regexpi(y(j).name, 'dB')-3));
            else
                LDSLevel = NaN;
            end
            HPCF = partblx.HPCF;
            LPCF = partblx.LPCF;
            ModDepth = partblx.ModDepth;
            ModFrequency = partblx.ModFrequencyHz; 
            ModExponent = partblx.ModExponent;
            clear partblx
            partblx = table(Time, LDSPreSEC, LDSDurSEC, LDSISISec, LDSNPulses, LDSLevel, HPCF, LPCF, ModDepth,...
                ModFrequency, ModExponent, BlockName);
        end
        partblx = addvars(partblx, expblock, exptime);
        partbl = [partbl; repmat(partblx, size(tbl, 1), 1)];
        tbl = addvars(tbl, onsetrate);
        utbl = [utbl; tbl];
    end
    
    AD = (utbl.ADSpikeCount_Cont);
%     logAD(~isfinite(logAD)) = NaN;
    
    utbl = [utbl partbl table(AD)];
    isSoundEvoked = utbl.onsetrate > poissinv(.99, utbl.AveRate_PSTH1_Hz*2.5)/2.5;
    utbl = addvars(utbl, isSoundEvoked);
% end
ADtbl = utbl(utbl.EarlySpikeCount> 0,:);
nonADtbl = utbl(utbl.EarlySpikeCount == 0,:);
ADsoundtbl = ADtbl(ADtbl.onsetrate > poissinv(.99, ADtbl.AveRate_PSTH1_Hz*2.5)/2.5, :);    
nonADsoundtbl = nonADtbl(nonADtbl.onsetrate > poissinv(.99, nonADtbl.AveRate_PSTH1_Hz*2.5)/2.5, :);


writetable(utbl, [outpath 'allchannels.csv'])
writetable(ADtbl, [outpath 'ADchannels.csv'])
writetable(ADsoundtbl, [outpath 'ADsoundevokedchannels.csv'])
writetable(nonADsoundtbl, [outpath 'nonADsoundevokedchannels.csv'])
writetable(nonADtbl, [outpath 'nonADchannels.csv'])