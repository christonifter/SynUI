y = dir('D:\ExperimentData\19-66-S2LDSsn\');
for i = 3:numel(y)
    table = readtable(['D:\ExperimentData\19-66-S2LDSsn\' y(i).name]);
    
    if rem(i,2)
        LDSstats(:,1:3,ceil((i-2)/2)) = [table.AveRate_PSTH2_Hz, table.ADSpikeCount_Cont, table.ADOnset_cont_sec];
    else
        LDSstats(:,4:5,ceil((i-2)/2)) = [table.AveRate_PSTH1_Hz, table.ADSpikeCount_Cont];
    end
end


LDSrates = squeeze(LDSstats(:,1,:));
LDSevoke = squeeze(LDSstats(:,2,:));
LDSsig = squeeze(LDSstats(:,3,:));

prerates = squeeze(LDSstats(:,4,:));
AUC = squeeze(LDSstats(:,5,:));
AUC(AUC==0) = NaN; %include only AD trials
% AUC(isnan(LDSsig)) = NaN; %include only sound-driven trials 


figure(1); clf;
semilogy(LDSrates(:), AUC(:), 'k.')
[r,p] = corr(LDSrates(:), log(AUC(:)), 'rows', 'complete')
ylabel('AUC')
xlabel('LDS rate')