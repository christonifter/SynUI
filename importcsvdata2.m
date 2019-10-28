if ~exist('utbl', 'var')
    spontfolders = dir('D:\ExperimentData\spontaneous\');
    spontparams = dir('D:\ExperimentData\spontparams\');
    LDSfolders = dir('D:\ExperimentData\LDSstats\');
    LDSPSTHfolders = dir('D:\ExperimentData\LDSPSTH\');
    utbl = [];
    partbl = [];
    LDSpeak = [];
    expseq = [];
    for i = 3:numel(spontfolders)
        y = dir(['D:\ExperimentData\spontaneous\' spontfolders(i).name]);
        z = dir(['D:\ExperimentData\spontparams\' spontparams(i).name]);
        w = dir(['D:\ExperimentData\LDSstats\' LDSfolders(i).name]);
        v = dir(['D:\ExperimentData\LDSPSTH\' LDSPSTHfolders(i).name]);

        yrn = str2double(spontfolders(i).name(1:2));
        subn = str2double(spontfolders(i).name(4:5));
        siten = str2double(spontfolders(i).name(8));

        for j = 3:numel(y)
            tbl = readtable(['D:\ExperimentData\spontaneous\' spontfolders(i).name '\' y(j).name]);
            partblx = readtable(['D:\ExperimentData\spontparams\' spontparams(i).name '\' z(j).name]);
            LDStbl = readtable(['D:\ExperimentData\LDSstats\' LDSfolders(i).name '\' w(j).name]);
            psthtbl = readtable(['D:\ExperimentData\LDSPSTH\' LDSPSTHfolders(i).name '\' v(j).name]);
            exporder = j-2;
            expyr = str2double(y(j).name((end-19):(end-18)));
            expmo = str2double(y(j).name((end-17):(end-16)));
            expday = str2double(y(j).name((end-15):(end-14)));
            exphr = str2double(y(j).name((end-12):(end-11)));
            expmin = str2double(y(j).name((end-10):(end-9)));
            exptime = datetime(expyr, expmo, expday, exphr, expmin, 0);
            partblx = partblx(:, 1:11);
            partblx = addvars(partblx, exporder, exptime);
            
            onsetrate = table2array(psthtbl(1,2:size(psthtbl,2)))';
            LDSpeak = [LDSpeak; LDStbl.MaxRate_PSTH2_Hz];
            partbl = [partbl; repmat(partblx, size(tbl, 1), 1)];
            tbl.Cluster = tbl.Cluster + yrn*1E7 + subn *1E4 + siten * 1E3;
            tbl = addvars(tbl, onsetrate);
            utbl = [utbl; tbl];
        end
    end
    AD = (utbl.ADSpikeCount_Cont);
%     logAD(~isfinite(logAD)) = NaN;

    utbl = [utbl partbl table(LDSpeak) table(AD)];

end
usoundtbl = utbl(utbl.onsetrate > poissinv(.99, utbl.AveRate_PSTH1_Hz*2.5)/2.5, :);
includeclusts = unique(utbl.Cluster(utbl.ADSpikeCount_Cont > 0));
ADtbl = utbl(ismember(utbl.Cluster, includeclusts), :);
ADsoundtbl = ADtbl(ADtbl.onsetrate > poissinv(.99, ADtbl.AveRate_PSTH1_Hz*2.5)/2.5, :);

umod = fitlm(ADsoundtbl, 'AD ~ LDSLevel')
figure(1); plot(ADsoundtbl.LDSLevel, ADsoundtbl.AD, 'k.')
hold on;
plot(1:100, (table2array(umod.Coefficients(1,1)) + table2array(umod.Coefficients(2,1)).*(1:100)), 'r--') 
hold off;
xlim([45 100])
box off;
% set(gca, 'YTick', [1 10 100 1000]);
% set(gca, 'YTickLabel', [1 10 100 1000]);
% set(gcf, 'Position', [100 100 300 300]);
xlabel('Sound Level (dB SPL)')
ylabel('LSA spikes')

umod = fitlm(ADsoundtbl, 'AD ~ LDSaverate')
figure(2); plot(ADsoundtbl.LDSaverate, ADsoundtbl.AD, 'k.')
hold on;
plot(1:100, (table2array(umod.Coefficients(1,1)) + table2array(umod.Coefficients(2,1)).*(1:100)), 'r--') 
hold off;
box off;
% set(gca, 'YTick', [1 10 100 1000]);
% set(gca, 'YTickLabel', [1 10 100 1000]);
% set(gcf, 'Position', [100 100 300 300]);
xlabel('Sound-evoked spike rate')
ylabel('LSA spikes')

[r, p] = corr(ADsoundtbl.LDSaverate, ADsoundtbl.AD, 'rows', 'complete')
[r, p] = corr(ADsoundtbl.LDSLevel, ADsoundtbl.AD, 'rows', 'complete');



uclusters = sort(unique(ADsoundtbl.Cluster));
[repcount, m] = hist(ADsoundtbl.Cluster, uclusters);
hirepclusts = uclusters(repcount>1);
hireps = repcount(repcount>1);

regtbl = ADsoundtbl(ismember(ADsoundtbl.Cluster, hirepclusts), :);
regtbl.Clusterc = categorical(regtbl.Cluster);
lmrate = fitlm(regtbl,'AD~LDSaverate*Clusterc');

figure(2); hold on;
for i = 2:numel(hirepclusts)
    slopei = table2array(lmrate.Coefficients(end/2+i,1)) + table2array(lmrate.Coefficients(2,1));
    intercepti = table2array(lmrate.Coefficients(1+i,1)) + table2array(lmrate.Coefficients(1,1));
     LDSratei = regtbl.LDSaverate(regtbl.Cluster == hirepclusts(i));
     ADi = regtbl.AD(regtbl.Cluster == hirepclusts(i));
     
     [~, sorti] = sort(LDSratei);
     plot([min(LDSratei) max(LDSratei)], [intercepti+slopei*min(LDSratei) intercepti+slopei*max(LDSratei)], '-')
     LDSratem(i,1) = mean(LDSratei);
     ADm(i,1) = mean(ADi);
end
hold off;
slopes = table2array(lmrate.Coefficients((end/2+1):end,1)) + table2array(lmrate.Coefficients(2,1));
figure(4);
hist(slopes, -5:.05:5)
xlim([-5 5]);
figure(5)
plot(LDSratem, ADm, 'k.')
[r p] = corr(LDSratem, ADm)


% ADmat = NaN(numel(hirepclusts), max(repcount));
% LDSratemat = NaN(numel(hirepclusts), max(repcount));
% LDSlevelmat = NaN(numel(hirepclusts), max(repcount));
% for i = 1:numel(hirepclusts)
%     subind = find(ADsoundtbl.Cluster == hirepclusts(i));
%     ADmat(i, 1:hireps(i)) = ADsoundtbl.logAD(subind);
%     LDSratemat(i, 1:hireps(i)) = ADsoundtbl.LDSaverate(subind);
%     LDSlevelmat(i, 1:hireps(i)) = ADsoundtbl.LDSLevel(subind);
% end
% regtbl = table(ADmat, LDSratemat, LDSlevelmat);
% regtbl = splitvars(regtbl);
% fitrm(regtbl, 'ADmat_1-ADmat_8 ~ LDSratemat_1-LDSratemat_8')



uclusters = sort(unique(utbl.Cluster));
usubs = sort(unique(floor(utbl.Cluster./1000)));

for i = 1:max(utbl.exporder)
    baselinem(i,1) = mean(utbl.AveRate_PSTH1_Hz(utbl.exporder== i));
    evokem(i,1) = mean(utbl.LDSaverate(utbl.exporder== i));
end
fatiguemodel = fitlm(utbl,'AveRate_PSTH1_Hz~exporder');

exptimeelapsed = NaN(size(utbl, 1), 1);

for i = 1:numel(usubs)
    thissubject = floor(utbl.Cluster./1000) == usubs(i);
    exptimeelapsed(thissubject) = (datenum(utbl.exptime(thissubject)) - min(datenum(utbl.exptime(thissubject)))).*24.*60;
end
utbl = addvars(utbl, exptimeelapsed);
plot(exptimeelapsed, utbl.AveRate_PSTH1_Hz, 'k.')
fatiguemodel2 = fitlm(utbl,'AveRate_PSTH1_Hz~exptimeelapsed');