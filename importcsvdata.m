y = dir('D:\ExperimentData\19-32-S2LDSsn\');
z = dir('D:\ExperimentData\19-32-S2LDSsn-params\');
utbl = [];
partbl = [];
for i = 3:numel(y)
    tbl = readtable(['D:\ExperimentData\19-32-S2LDSsn\' y(i).name]);
    partblx = readtable(['D:\ExperimentData\19-32-S2LDSsn-params\' z(i).name]);
    partbl = [partbl; repmat(partblx, size(tbl, 1), 1)];
    utbl = [utbl; tbl];
%     if rem(i,2)
%         LDSstats(:,1:3,ceil((i-2)/2)) = [tbl.AveRate_PSTH2_Hz, tbl.ADSpikeCount_Cont, tbl.ADOnset_cont_sec];
%     else
%         LDSstats(:,4:5,ceil((i-2)/2)) = [tbl.AveRate_PSTH1_Hz, tbl.ADSpikeCount_Cont];
%     end
end

utbl = [utbl partbl];
rndv = normrnd (10,2,size(utbl,1), 1);
utbl = addvars(utbl, rndv);
includeclusts = unique(utbl.Cluster(utbl.ADSpikeCount_Cont > 0));
ADtbl = utbl(ismember(utbl.Cluster, includeclusts), :);
umod = fitlm(ADtbl, 'ADSpikeCount_Cont ~ LDSaverate+LDSLevel +Cluster')

plot(ADtbl.LDSLevel, ADtbl.ADSpikeCount_Cont, 'k.')


% LDSrates = squeeze(LDSstats(:,1,:));
% LDSevoke = squeeze(LDSstats(:,2,:));
% LDSsig = squeeze(LDSstats(:,3,:));
% 
% prerates = squeeze(LDSstats(:,4,:));
% AUC = squeeze(LDSstats(:,5,:));
% 
% 
% % AUC(AUC==0) = NaN; %include only AD trials
% % AUC(isnan(LDSsig)) = NaN; %include only sound-driven trials 
% 
% 
% figure(1); clf;
% xvar = LDSrates(51,:)';
% yvar = AUC(51,:)';
% plot(xvar, yvar, 'k.')
% [r,p] = corr(xvar, yvar, 'rows', 'complete')
ylabel('AUC')
xlabel('LDS rate')