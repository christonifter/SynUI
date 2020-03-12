%Streams data from synapse to MATLAB during preview and recording
function livestreamapp4(app, mode) 
app.DialogueLabel.Text = 'Beginning Experiment';
readsynpars(app);
chanlist = updatechans(app);
spikevar = 'eSpk';
wavevar = 'RAW1';
stimvar = 'Stim';
stimvar2 = 'Stmb';
stimvarf = 'Freq';
offvar = 'StOn';
scalpvar = 'RA4W';
evvar = 'EvOn';

app.DialogueLabel.Text = 'Connecting to SynapseLive';

if app.syn.getMode() == 0
    c = clock();
    d = [num2str(mod(c(1), 100)) sprintf('%02d', c(2)), sprintf('%02d', c(3))];
    t = [sprintf('%02d', c(4)) sprintf('%02d', c(5))];
    blockname = [sprintf('%03d', app.BlockEdit.Value) '-' app.syn.getCurrentExperiment() '-' d '-' t];
    app.syn.setCurrentBlock(blockname);
end

exp = SynapseLive('MODE', mode, 'EXPERIMENT', app.ExperimentDrop.Value);
updaterate = app.UpdateEdit.Value;
binwindow = app.BinEdit.Value./1000; %PSTH bin size (sec)
updatemode(app);
if updaterate<2*binwindow; pause(2*binwindow); end

spets = [];
channels = [];
stimons = [];
stimoffs = [];
evons = [];
evoffs = [];
frqs = [];
lvls = [];
currtime = [];
parmat = [];
sttrace = [];
rawdata = [];
app.DialogueLabel.Text = 'Streaming Data';
loopc = 1;
loopcon = 1;
while loopcon
    pause(updaterate)

    if app.syn.getMode() < 2 %stops if synapse is set to idle or standby
        break
    end
    new = exp.update;
    try 
    if isempty(new)
        break
    end

%read UI data  
    updaterate = app.UpdateEdit.Value;
    binwindow = app.BinEdit.Value./1000; %PSTH bin size (sec)
    ratewin = app.RangeEdit.Value;
    analwin = [app.OnsetEdit.Value app.OffsetEdit.Value]./1000; %ftc analysis window;
    chanlist = updatechans(app);


%read continuously updated data
    currtime = [currtime; app.syn.getParameterValue('ReadTime', 'Time')];
    currpars = NaN(size(app.pars));
    for par = 1:numel(app.pars)
        currpars(par) = app.syn.getParameterValue(char(app.gizmos(par)), char(app.pars(par)));
    end
    parmat = [parmat; currpars'];
    
    if isfield(new.data.streams, scalpvar)
        rawdata = [rawdata; new.data.streams.(scalpvar).data'];
        fs = new.data.streams.(scalpvar).fs;
    end
%     if isfield(new.data.streams, scalpvar)
%         chans = new.data.streams.(scalpvar).chan;
%         nchans = numel(unique(chans));
%         
%         if rem(numel(chans), nchans)>0
%             disp('data collection problem: unequal trials among channels');
%         end
%         for chan = 1:nchans
%             chani = find(chans == chan);
%             newchandata(:,:,chan) = new.data.streams.(scalpvar).data(chani, :);
%         end
%         sttrace = [sttrace; newchandata];
%         app.data.sttrace = sttrace;
%         app.data.fs = new.data.streams.(scalpvar).fs;
%         clear newchandata
%     end
%     
    if isfield(new.data.snips, spikevar)
        spets = [spets; new.data.snips.(spikevar).ts];
        channels = [channels; new.data.snips.(spikevar).chan];
    end
    
    gizlist = app.syn.getGizmoNames();
    threshlist = [];
    if ~isempty(cell2mat(regexp(gizlist, 'Neu1')))
        nthreshs = sum(cell2mat(regexp(app.syn.getParameterNames('Neu1'), 'Threshold')));
        threshlist = zeros(nthreshs,1);
        for i = 1:nthreshs
            threshlist(i,1) = app.syn.getParameterValue('Neu1', ['Threshold' num2str(i)]);
        end
    end

%read triggered data
    nstims = 0;
    if ~sum(isnan(new.data.epocs.(offvar).data))
        stimons = [stimons; new.data.epocs.(offvar).onset(new.data.epocs.(offvar).data == 1)];
        stimoffs = [stimoffs; new.data.epocs.(offvar).onset(new.data.epocs.(offvar).data == 0)];
        nstims = numel(new.data.epocs.(offvar).onset(new.data.epocs.(offvar).data == 0));
    end
    if isfield(new.data.streams, scalpvar)
        if ~sum(isnan(new.data.epocs.(evvar).data))
            evons = [evons; new.data.epocs.(evvar).onset(new.data.epocs.(evvar).data == 1)];
            evoffs = [evoffs; new.data.epocs.(evvar).onset(new.data.epocs.(evvar).data == 0)];
            nevs = numel(new.data.epocs.(evvar).onset(new.data.epocs.(evvar).data == 0));
        end
    end
    if isfield(new.data.streams, scalpvar)
        usechan = 1;
        startsamp = round(evons.*fs);
        endsamp = startsamp + round(.02*fs);
        if size(rawdata, 1)<endsamp(end)
            startsamp = startsamp(1:(end-1));
            endsamp = endsamp(1:(end-1));
        end
         for cycle = 1:numel(startsamp)
            sttrace(cycle,:) = rawdata(startsamp(cycle):endsamp(cycle), usechan);
         end
%         app.data.fs = fs;
    end
    
    if isfield(new.data.epocs, stimvar) %FRA
        frqs = [frqs; new.data.epocs.(stimvar).data(~isnan(new.data.epocs.(stimvar).data))];
    elseif isfield(new.data.epocs, stimvarf) %Iso-I
        frqs = [frqs; new.data.epocs.(stimvarf).data(~isnan(new.data.epocs.(stimvarf).data))];
    elseif ismember('Frequency', app.pars)
        frqparind = find(strcmpi(app.pars, 'Frequency'));
        frqs = [frqs; app.syn.getParameterValue(char(app.gizmos(frqparind)), 'Frequency').* ones(nstims, 1)];
    else
        frqs = [frqs; zeros(nstims, 1)];
    end
    
    if isfield(new.data.epocs, stimvar2) %FRA
        lvls = [lvls; new.data.epocs.(stimvar2).data(~isnan(new.data.epocs.(stimvar2).data))];
    end
    if ismember('LDSAmplitude', app.pars)
        lvlparind = find(strcmpi(app.pars, 'LDSAmplitude'));
        lvls = [lvls; app.syn.getParameterValue(char(app.gizmos(lvlparind)), 'LDSAmplitude').* ones(nstims, 1)];
    end

    if ismember('Level', app.pars)
        lvlparind = find(strcmpi(app.pars, 'Level'));
        lvls = [lvls; app.syn.getParameterValue(char(app.gizmos(lvlparind)), 'Level').* ones(nstims, 1)];
    end

    endtext = '';
    if app.LastoffsetnsecButton.Value
        x = app.StimDurLabel.Text;
        y = datestr(seconds(0), 'MM:SS');
        stimdur = seconds(datetime(x, 'InputFormat', 'mm:ss')-datetime(y, 'InputFormat', 'mm:ss'));
        endtime = stimdur;
        endtext = ['(' datestr(seconds(endtime-currtime(end)),'MM:SS') ' remaining) '];
    end
    if app.SecondsButton.Value
        endtime = app.ExpDurationEdit.Value;
        endtext = ['(' datestr(seconds(endtime-currtime(end)),'MM:SS') ' remaining)'];
    end
    
    if isfield(new.data.streams, scalpvar)
        reporttrials = evons;
    else
        reporttrials = stimons;
    end

    app.DialogueLabel.Text = ['Streaming Data', newline, datestr(seconds(currtime(end)),'MM:SS'), ' ' ...
        endtext, num2str(loopc) ' polls, ' num2str(numel(reporttrials)) ' trials' ];
if isfield(new.data.streams, scalpvar)
    if isempty(stimoffs)
        ax = app.LFP1Axes;
        traces = sttrace;
    else
        ax = app.LFP2Axes;
        traces = sttrace(startsamp > stimons(1)*fs, :);
    end
    sttracesplot(app, ax, traces, fs)
else
%running spike rate average
    if numel(chanlist)>0
        bincount = NaN(round(currtime(end)/binwindow), max(chanlist));
        for chan = chanlist'
            [bincount(:,chan), bintime] = hist(spets(channels==chan), 0.5*binwindow:binwindow:currtime(end));
        end
        if app.fratescaleEdit.Value == 0
            yscale = max(max(bincount));
        else
            yscale = app.fratescaleEdit.Value*binwindow;
        end
        frate = bincount./yscale - (1:max(chanlist));
    app.frateTable.Data = [(1:max(chanlist))', max(bincount)'./binwindow];
    app.frateTable.ColumnName = {'Channel', 'Max rate'};
    else
        frate = NaN;
        bintime = NaN;
    end

    plotfrate(app, frate, bintime, stimons, stimoffs, chanlist)
     xlim(app.rateAxes, [max([0 currtime(end)-ratewin]), max([1 currtime(end)])]);                

%tuning curve
     if numel(stimoffs)>0
        plottrials = min([numel(stimons) numel(frqs) numel(lvls)]);
        [spetfreq, spetlevel, stimspets, trials] = freqlevelspet(spets, stimons(1:plottrials), frqs(1:plottrials), lvls(1:plottrials));
        if (numel(unique(lvls)) == 1) && (numel(unique(frqs))>1)
            trials = reshape(trials, 1, numel(trials));
            synftc(stimspets, spetfreq, analwin, channels, trials, chanlist, frqs, app.FTCAxes);
            set(app.FTCAxes,'XTickLabel', round(sort(unique(frqs))./100)./10)
            xlabel(app.FTCAxes, 'Frequency (kHz)')
            title(app.FTCAxes, 'Iso Intensity Function')
        end
        if (numel(unique(frqs)) == 1) && (numel(unique(lvls))>1)
            trials = reshape(trials, 1, numel(trials));
            synftc(stimspets, spetlevel, analwin, channels, trials, chanlist, lvls, app.FTCAxes);
            set(app.FTCAxes,'XTickLabel', sort(unique(lvls)))
            xlabel(app.FTCAxes, 'Level (dB)')
            title(app.FTCAxes, 'Rate Level Function')
        end
        if (numel(unique(lvls))>1) && (numel(unique(frqs))>1)
%             plottrials = min([numel(stimons) numel(frqs) numel(lvls)]);
%             [spetfreq, spetlevel, stimspets, trials] = freqlevelspet(spets, stimons(1:plottrials), frqs(1:plottrials), lvls(1:plottrials));
            multifra(stimspets, spetfreq, spetlevel, analwin, channels, trials, chanlist, frqs(1:plottrials), lvls(1:plottrials), app.FTCAxes);
        end
     end
end
loopc = loopc + 1;
if (~app.InfiniteButton.Value) && (currtime(end) > endtime)

        loopcon = 0;
        app.syn.setMode(0);
         app.StimPanel.Visible = 0;
        updatemode(app);
end
    catch ME
        ME
        ME.stack(1)
        ME.stack(2)
         [numel(stimons) numel(frqs) numel(lvls)]
    end
end %while

metext = 'Ready.';
if exist('ME', 'var')
    metext = 'A MATLAB Error was thrown :(';
end
app.DialogueLabel.Text = ['Data streaming ended. ' metext];

if strcmpi(mode, 'Record')
    try
        stimpanelvalues = [app.StimOnEdit.Value; app.StimOffEdit.Value; app.StimOnEdit.Value; ...
            app.TrainOffEdit.Value; app.TrainRepeatsEdit.Value; app.StimFreqEdit.Value; ...
            app.StimLevelEdit.Value; app.LDSPreEdit.Value; app.LDSDurEdit.Value; ...
            app.LDSISIEdit.Value; app.LDSRepeatsEdit.Value; app.LDSPostGapEdit.Value; ...
            app.LDSCFreqEdit.Value; app.LDSBWEdit.Value; app.LDSLevelEdit.Value; ...
            app.LDSModDepthEdit.Value; app.LDSModExpEdit.Value; app.LDSModFreqEdit.Value];
        stimpanelnames = [{'StimOn'}; {'StimOff'}; {'TrainOn'}; {'TrainOff'}; {'TrainRepeats'}; {'StimFreq'}; ...
            {'StimLevel'}; {'LDSPre'}; {'LDSDur'}; {'LDSISI'}; {'LDSRepeats'}; {'LDSPostGap'}; ...
            {'LDSCFreq'}; {'LDSBW'}; {'LDSLevel'}; {'LDSModDepth'}; {'LDSModExp'}; {'LDSModFreq'}];
        stimpanelvis = [strcmp(app.StimOnEdit.Visible, 'on'); strcmp(app.StimOffEdit.Visible, 'on'); strcmp(app.TrainOnEdit.Visible, 'on'); ...
            strcmp(app.TrainOffEdit.Visible, 'on'); strcmp(app.TrainRepeatsEdit.Visible, 'on'); strcmp(app.StimFreqEdit.Visible, 'on'); ...
            strcmp(app.StimLevelEdit.Visible, 'on'); strcmp(app.LDSPreEdit.Visible, 'on'); strcmp(app.LDSDurEdit.Visible, 'on'); ...
            strcmp(app.LDSISIEdit.Visible, 'on'); strcmp(app.LDSRepeatsEdit.Visible, 'on'); strcmp(app.LDSPostGapEdit.Visible, 'on'); ...
            strcmp(app.LDSCFreqEdit.Visible, 'on'); strcmp(app.LDSBWEdit.Visible, 'on'); strcmp(app.LDSLevelEdit.Visible, 'on'); ...
            strcmp(app.LDSModDepthEdit.Visible, 'on'); strcmp(app.LDSModExpEdit.Visible, 'on'); strcmp(app.LDSModFreqEdit.Visible, 'on')];
        if strcmp(app.StimulusPanel.Visible, 'off')
            stimpanelvis(1:7) = 0;
        end
        if strcmp(app.LDSPanel.Visible, 'off')
            stimpanelvis(8:end) = 0;
        end
        stimval = stimpanelvalues(find(stimpanelvis));
        stimname = stimpanelnames(find(stimpanelvis));
        stimtable = table(stimname, stimval);
    catch ME
        stimtable = table([]);
        ME
        disp('unable to save synui stimulus values');
    end

    
    [changesamp, ~]  = find(diff(parmat));
    changes = sort(unique(changesamp))+1;
    fulltable = array2table([currtime parmat], 'VariableNames', [{'Time'}; app.pars]);
    changearray = [[currtime(1), parmat(1,:)]; [currtime(changes),  parmat(changes,:)]];
    changetable = array2table(changearray, 'VariableNames', [{'Time'}; app.pars]);

    flist = dir(app.syn.getCurrentTank);
    folders = zeros(size(flist));
    dates = zeros(size(flist));
    for file = 3:numel(flist)
        folders(file) = flist(file).isdir;
        dates(file) = flist(file).datenum;
    end
    subind= find(folders);
     [~, i]=max(dates(subind));
     currfolder = [flist(subind(i)).folder '\' flist(subind(i)).name];
     save([currfolder '\params.mat'], 'fulltable', 'changetable', 'stimtable', 'threshlist')
     app.BlockEdit.Value = app.BlockEdit.Value + 1;
end


end %function livestreamapp