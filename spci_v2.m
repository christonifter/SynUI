% rules:
% 1. This script will take all files in the current working directory and
% average certain columns together.

% 2. If you want input files ordered correctly, use the same number of digits in
% each input filename. Example file1: 'xx_2_7s' -> 'xx_02_07s'
% This script will read the files in the current working folder, and sort 
% them alphabetically. and 'xx_2_7s' will be put behind 'xx_12_17s'

% 3. Use only characters, digits and underscores in your input filenames. This
%script will try to label the columns with the filenames the columns are
%taken from. MATLAB requires certain rules when naming column names. If you
%break this rule, the script will work, but the column names won't include the
%source filenames.


clear
outpath = 'D:\ExperimentData\tone_evoked\';
outname = '19-99_007_19kHz_50-60dB_blockave';
y = dir(pwd);
for filei = 1:(numel(y)-2)
    ss = readtable(y(filei+2).name);
    sp1(:, filei) = ss.SpectralPowerMod1;
    sp2(:, filei) = ss.SpectralPowerMod2;
    [~, fni, ~] = fileparts(y(filei+2).name);
    fn(1, filei) = {['post_sp' num2str(filei) '_' fni]};
    cfn(1, filei) = {['change' num2str(filei) '_' fni]};
end
Channel = ss.Channel;
for row = 1:size(sp2, 1)
    pre_mean(row, 1) = mean(sp1(row,:));
    df = size(sp2, 2)-1;
    pre_SEM(row, 1) = std(sp1(row,:))/sqrt(df+1); % Standard Error
    for col = 1:size(sp2, 2)
        [change(row,col), p(row, col), CI2(row,:)] = ttest(sp1(row, :), sp2(row, col));
        if pre_mean(row,1)>sp2(row, col)
            change(row,col) = change(row,col)*-1;
        end
    end
end
out = table(Channel, pre_mean, pre_SEM, CI2, sp2, change);
out2 = out;
try
    out2 =splitvars(out, {'sp2', 'change'}, 'NewVariableNames', {fn, cfn});
catch ME
    ME
end
writetable(out2, [outpath outname '.csv']);
