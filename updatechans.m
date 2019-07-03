function [chanlist, channelsortorder] = updatechans(app)
    if app.ClustsCheck.Value == 1
        cluster = app.data.cluster;
        data = app.data;
        chanlist = (1:numel(cluster))';
        channelsortorder = chanlist;
        clustampmat = NaN(numel(cluster), 1);
        spikecount = NaN(numel(cluster), 1);
        AD = NaN(numel(cluster), 1);

        analwin(1,:) = [0 data.stimons(1)];
        analwin(2,:) = [data.stimoffs(1) size(data.LFP, 1)/data.fs];
        for i = 1:numel(cluster)
            clustampmat(i,:) = cluster(i).peakChannel2(1);
            spikecount(i) = numel(cluster(i).spikes);
            postcount = sum(cluster(i).spikes>analwin(2,1) & cluster(i).spikes<analwin(2,2));
            posttime = analwin(2,2) - analwin(2,1);
            postrate = postcount/posttime;
            precount = sum(cluster(i).spikes>analwin(1,1) & cluster(i).spikes<analwin(1,2));
            pretime = analwin(1,2) - analwin(1,1);
            prerate = precount/pretime;
            if postrate + prerate == 0
                den = 1;
            else
                den = sqrt(postrate + prerate);
            end
            AD(i) = (postrate - prerate) / den;
        end
        switch app.SortOrderDropDown.Value
            case 'Nearest Channel'
                [~, channelsortorder] = sort(clustampmat);
            case 'Spike Count'
                [~, channelsortorder] = sort(spikecount, 'descend');
            case 'AD'
                [ADsort, channelsortorder] = sort(AD, 'descend');
        end
        [~, sorti] = sort(channelsortorder);
        if app.AllRadio.Value
        elseif app.TopRadio.Value
            chanlist = (1:str2double(app.ClustTopEdit.Value))';
        elseif app.BottomRadio.Value
            chanlist = ((numel(cluster)-str2double(app.ClustBottomEdit.Value)+1):numel(cluster))';
        elseif app.ClustRadio.Value
%             chanlist = eval(['[' app.ClustClustEdit.Value ']']);
%             chanlist = sorti(chanlist);
            usechannels = eval(['[' app.ClustClustEdit.Value ']']);
            nousechannels = setdiff(channelsortorder, usechannels);
            
            chanlist = 1:numel(usechannels);
            channelsortorder= [usechannels'; nousechannels];
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
    chanlist = sort(chanlist(:));
end
