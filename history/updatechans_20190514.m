function [chanlist, channelsortorder] = updatechans(app)
    if app.ClustsCheck.Value == 1;
        cluster = app.data.cluster;
        chanlist = (1:numel(cluster))';
        channelsortorder = chanlist;
        clustampmat = NaN(numel(cluster), 1);
        spikecount = NaN(numel(cluster), 1);
        AD = NaN(numel(cluster), 1);
        for i = 1:numel(cluster)
            clustampmat(i,:) = cluster(i).peakChannel2(1);
            spikecount(i) = numel(cluster(i).spikes);
            postrate = (1+sum(cluster(i).spikes>data.stimoffs(1)))/(1+max(data.spets)-data.stimoffs(1));
            prerate = (1+sum(cluster(i).spikes<data.stimons(1)))/(1+data.stimons(1));
            AD(i) = postrate/prerate;
        end
        switch app.SortOrderDropDown.Value
            case 'Nearest Channel'
                [~, channelsortorder] = sort(clustampmat);

            case 'Spike Count'
                [~, channelsortorder] = sort(spikecount, 'descend');
            case 'AD'
                [~, channelsortorder] = sort(AD, 'descend');
        end
        
        if app.AllRadio.Value
        elseif app.TopRadio.Value
            chanlist = channelsortorder(1:str2double(app.ClustTopEdit.Value));
        elseif app.BottomRadio.Value
            chanlist = channelsortorder((end-str2double(app.ClustBottomEdit.Value)):end);
        elseif app.ClustRadio.Value
            chanlist = eval(['[' app.ClustClustEdit.Value ']']);
        end
    else
        chanvec = zeros(32,1);
        for chan = 1:32
            chanfield = ['CheckBox_' num2str(chan)];
            if app.(chanfield).Value
                chanvec(chan) = 1;
            end
        end
        chanlist = find(chanvec);
        channelsortorder = 1:numel(chanlist);
    end
end
