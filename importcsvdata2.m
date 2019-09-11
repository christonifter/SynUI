% spontfolders = dir('D:\ExperimentData\spontaneous\');
% spontparams = dir('D:\ExperimentData\spontparams\');
% LDSfolders = dir('D:\ExperimentData\LDSstats\');
% LDSPSTHfolders = dir('D:\ExperimentData\LDSPSTH\');
% utbl = [];
% partbl = [];
% LDSpeak = [];
% for i = 3:numel(spontfolders)
%     y = dir(['D:\ExperimentData\spontaneous\' spontfolders(i).name]);
%     z = dir(['D:\ExperimentData\spontparams\' spontparams(i).name]);
%     w = dir(['D:\ExperimentData\LDSstats\' LDSfolders(i).name]);
%     v = dir(['D:\ExperimentData\LDSPSTH\' LDSPSTHfolders(i).name]);
%     
%     yrn = str2double(spontfolders(i).name(1:2));
%     subn = str2double(spontfolders(i).name(4:5));
%     siten = str2double(spontfolders(i).name(8));
%     
%     for j = 3:numel(y)
%         tbl = readtable(['D:\ExperimentData\spontaneous\' spontfolders(i).name '\' y(j).name]);
%         partblx = readtable(['D:\ExperimentData\spontparams\' spontparams(i).name '\' z(j).name]);
%         LDStbl = readtable(['D:\ExperimentData\LDSstats\' LDSfolders(i).name '\' w(j).name]);
%         psthtbl = readtable(['D:\ExperimentData\LDSPSTH\' LDSPSTHfolders(i).name '\' v(j).name]);
%         onsetrate = table2array(psthtbl(1,2:size(psthtbl,2)))';
%         LDSpeak = [LDSpeak; LDStbl.MaxRate_PSTH2_Hz];
%         partbl = [partbl; repmat(partblx(:, 1:11), size(tbl, 1), 1)];
%         tbl.Cluster = tbl.Cluster + yrn*1E7 + subn *1E4 + siten * 1E3;
%         tbl = addvars(tbl, onsetrate);
%         utbl = [utbl; tbl];
%     end
% end
% 
% utbl = [utbl partbl table(LDSpeak)];
includeclusts = unique(utbl.Cluster(utbl.ADSpikeCount_Cont > 0));
ADtbl = utbl(ismember(utbl.Cluster, includeclusts), :);
ADsoundtbl = ADtbl(ADtbl.onsetrate > poissinv(.99, ADtbl.AveRate_PSTH1_Hz*2.5)/2.5, :);
umod = fitlm(ADsoundtbl, 'ADSpikeCount_Cont ~ LDSaverate+LDSLevel +Cluster')
figure(1); plot(ADsoundtbl.LDSLevel, ADsoundtbl.ADSpikeCount_Cont, 'k.')
% hold on;
% semilogy(1:120, table2array(umod.Coefficients(1,1)) + table2array(umod.Coefficients(3,1)).*(1:120), 'r--') 
% hold off;
xlim([-5 100])
umod = fitlm(ADsoundtbl, 'ADSpikeCount_Cont ~ LDSLevel + Cluster')
figure(2); plot(ADsoundtbl.LDSaverate, ADsoundtbl.ADSpikeCount_Cont, 'k.')


[r, p] = corr(ADsoundtbl.LDSaverate, ADsoundtbl.ADSpikeCount_Cont)

[r, p] = corr(ADsoundtbl.LDSLevel, ADsoundtbl.ADSpikeCount_Cont)