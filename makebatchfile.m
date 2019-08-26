% tankpath = uigetdir('D:\Synapse\Tanks');
clustpath = 'D:\Spikes\';
outpath = 'D:\ExperimentData\';
batchpath = 'D:\Batches\';

[filepath, sub, ~] = fileparts(tankpath);
y = dir(tankpath);
c = 0;
clear flist
for i = 3:numel(y)
    if y(i).isdir
        c = c +1;
        tanklist{c,1} = [tankpath '\' y(i).name '\'];
        clusterlist{c,1} = [clustpath sub '\' y(i).name '.mat'];
        outlist{c,1} = [outpath sub '\' y(i).name];
    end
end

xlswrite([batchpath sub], [tanklist clusterlist outlist]);