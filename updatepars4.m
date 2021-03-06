%Reads parameters stored in synapse        
function updatepars3(app) 
    %Reset things
    if app.syn.getMode == 0
        app.syn.setMode(1);
    end
    for i = 1:18
        synvaluename = ['Param' num2str(i) 'Label'];
        editname = ['Param' num2str(i) 'Edit'];
        paramname = ['Par' num2str(i) 'Label'];
        app.(synvaluename).Visible = 0;
        app.(paramname).Visible = 0;
        app.(editname).Visible = 0;
    end
    
    app.TonePanel.Visible = 0;
    app.LDSPanel.Visible = 0;
    gizlist = app.syn.getGizmoNames();
    c = regexp(gizlist, '^x');
    pargizs = gizlist(~cellfun('isempty', c));
    c = 1;
    app.pars = {};
    app.gizmos = {};
    
%     app.ToneOnEdit.Value = 0;
%     app.ToneOffEdit.Value = 0;
%     app.TrainOnEdit.Value = 0;
%     app.TrainOffEdit.Value = 0;
%     app.TrainRepeatsEdit.Value = 0;
%     app.ToneFreqEdit.Value = 0;
%     app.ToneLevelEdit.Value = 0;
%     app.LDSPreEdit.Value = 0;
%     app.LDSDurEdit.Value = 0;
%     app.LDSISIEdit.Value = 0;
%     app.LDSRepeatsEdit.Value = 0;
%     app.PostGapEdit.Value = 0;
%     app.CenterFreqEdit.Value = 0;
%     app.BandwidthEdit.Value = 0;
%     app.TRMSEdit.Value = 0;
%     app.ModDepthEdit.Value = 0;
%     app.ModExpEdit.Value = 0;
%     app.ModFreqEdit.Value = 0;
    
    %collect gizmo and param names
    for gizmo = 1:numel(pargizs)
       params = app.syn.getParameterNames(char(pargizs(gizmo)));
       for param = 1:numel(params)
           inf = app.syn.getParameterInfo(char(pargizs(gizmo)), char(params(param)));
           if ~strcmpi(inf.Type, 'Logic')
               app.gizmos{c,1} = char(pargizs(gizmo));
               app.pars{c,1} = char(params(param));
               c = c+1;
           end
       end
    end
    %Hide things
    if sum([ismember(app.pars, 'StimHiMS'); ismember(app.pars, 'TrainDurationMS')])
        app.TonePanel.Visible = 1;
        if sum(ismember(app.pars, 'StimHiMS'))
            app.ToneOnEdit.Visible = 1;
            app.ToneOnLabel.Visible = 1;
            app.ToneOffEdit.Visible = 1;
            app.ToneOffLabel.Visible = 1;
        else
            app.ToneOnEdit.Visible = 0;
            app.ToneOnLabel.Visible = 0;
            app.ToneOffEdit.Visible = 0;
            app.ToneOffLabel.Visible = 0;            
        end
        if sum(ismember(app.pars, 'TrainDurationMS'))
            app.TrainOnEdit.Visible = 1;
            app.TrainOnLabel.Visible = 1;
            app.TrainOffEdit.Visible = 1;
            app.TrainOffLabel.Visible = 1;
        else
            app.TrainOnEdit.Visible = 0;
            app.TrainOnLabel.Visible = 0;
            app.TrainOffEdit.Visible = 0;
            app.TrainOffLabel.Visible = 0;            
        end
        if sum(ismember(app.pars, 'StimFrequencyHz'))
            app.ToneFreqEdit.Visible = 1;
            app.ToneFreqLabel.Visible = 1;
        elseif sum(ismember(app.pars, 'AMHPCF'))
                app.ToneFreqEdit.Visible = 1;
                app.ToneFreqLabel.Visible = 1;
                app.ToneBWEdit.Visible = 1;
                app.ToneBWLabel.Visible = 1;
        else
            app.ToneFreqEdit.Visible = 0;
            app.ToneFreqLabel.Visible = 0;
            app.ToneBWEdit.Visible = 0;
            app.ToneBWLabel.Visible = 0;
        end
        if sum(ismember(app.pars, 'StimAmplitude'))
            app.ToneLevelEdit.Visible = 1;
            app.ToneLevelLabel.Visible = 1;
        else
            app.ToneLevelEdit.Visible = 0;
            app.ToneLevelLabel.Visible = 0;
        end
        if sum(ismember(app.pars, 'Stim2FrequencyHz'))
            app.Tone2FreqEdit.Visible = 1;
            app.Tone2FreqLabel.Visible = 1;
            app.Tone2LevelEdit.Visible = 1;
            app.Tone2LevelLabel.Visible = 1;
        else
            app.Tone2FreqEdit.Visible = 0;
            app.Tone2FreqLabel.Visible = 0;
            app.Tone2LevelEdit.Visible = 0;
            app.Tone2LevelLabel.Visible = 0;
        end
        if sum(ismember(app.pars, 'TrainRepeats'))
            app.TrainRepeatsEdit.Visible = 1;
            app.TrainRepeatsLabel.Visible = 1;
        else
            app.TrainRepeatsEdit.Visible = 0;
            app.TrainRepeatsLabel.Visible = 0;
        end

    end
    if sum(ismember(app.pars, 'LDSDurationMS'))
        app.LDSPanel.Visible = 1;
        if sum(ismember(app.pars, 'LDSRepeats'))
            app.LDSISIEdit.Visible = 1;
            app.LDSRepeatsEdit.Visible = 1;
            app.LDSRepeatsLabel.Visible = 1;
            app.PostGapEdit.Visible = 0;
            app.LDSISILabel.Visible = 1;
            app.PostGapLabel.Visible = 0;
        else
            app.LDSISIEdit.Visible = 0;
            app.LDSRepeatsEdit.Visible = 0;
            app.LDSRepeatsLabel.Visible = 0;
            app.PostGapEdit.Visible = 1;
            app.LDSISILabel.Visible = 0;
            app.PostGapLabel.Visible = 1;
        end
    end
    if sum(ismember(app.pars, 'LDSPreMS'))
        app.LDSPreEdit.Visible = 1;
        app.LDSPreLabel.Visible = 1;
    else
        app.LDSPreEdit.Visible = 0;
        app.LDSPreLabel.Visible = 0;
    end        

    if sum(ismember(app.pars, 'LDSAmplitude'))
        app.TRMSEdit.Visible = 1;
        app.TRMSLabel.Visible = 1;
    else
        app.TRMSEdit.Visible = 0;
        app.TRMSLabel.Visible = 0;
    end

    if sum(ismember(app.pars, 'buffer2'))
        app.DialogueLabel.Text = 'PC -> Audio -> RZ6 (~1 min)';
        pause(.1)
        bufferduration = 2000; %in ms, use an even number of seconds.
        bufferlength = round(bufferduration*195.3125);
        fid = fopen('C:\TDT\Synapse\DMR500ic120s30dB.f32');
        fsignal = fread(fid, bufferlength, 'float');
        nsignal = fsignal./max(fsignal);
        t = (1:(10*192000))./192000;
        nsignal = sin(2*pi*400*t);
        app.syn.setParameterValue('xDMR_audio','BufferLengthMS',bufferduration); 
        app.syn.setParameterValue('xDMR_audio','TrainRepeats',6); 
        app.syn.setParameterValues('xDMR_audio','buffer1',nsignal); 
        app.syn.setParameterValues('xDMR_audio','buffer2',nsignal);
        app.DialogueLabel.Text = 'Loading';
        pause(.1)
    end

    
    % build stim preview    
    hitime = 0;
    lotime = 0;
    ncyc = 1;
    nrepeats = 1;
    nldsrepeats = 1;
    t1 = 0;
    t2 = 0;
    gap = 0;
    gap2 = 0;
    ldsgap = 0;
    for par = 1:numel(app.pars)
        paramname = ['Par' num2str(par) 'Label'];
        synvaluename = ['Param' num2str(par) 'Label'];
        editname = ['Param' num2str(par) 'Edit'];
        app.(synvaluename).Visible = 1;
        app.(editname).Visible = 1;
        app.(paramname).Visible = 1;
        app.(paramname).Text = char(app.pars(par));
        app.(synvaluename).Text = num2str(app.syn.getParameterValue(char(app.gizmos(par)), char(app.pars(par))));
        parval = app.syn.getParameterValue(char(app.gizmos(par)), char(app.pars(par)));
        switch char(app.pars(par))
            case 'StimHiMS'
                hitime = parval/1000;
            case 'StimLoMS'
                lotime = parval/1000;
%             case 'StimNPulses'
%                 ncyc = parval;
            case 'TrainDurationMS'
                t1 = parval/1000;
            case 'GateLoMS'
                gap = parval/1000;
            case 'LDSISIMS'
                ldsgap = parval/1000;
            case 'TrainGapMS'
                gap = parval/1000;
            case 'PostLDSGapMS'
                gap2 = parval/1000;
            case 'TrainRepeats'
                nrepeats = parval;
                ncyc = parval;
            case 'GateNPulses'
                nrepeats = parval;
            case 'LDSPreMS' 
                t1 = parval/1E3;
            case 'LDSDurationMS'
                t2 = parval/1000;
            case 'LDSRepeats'
                nldsrepeats = parval;
            case 'BufferLengthMS'
                hitime = parval/1000;
                lotime = parval/1000;
        end
    end
    if ~isempty(regexp(app.ExperimentDrop.Value, 'FRAho'))
        ncyc = 81;
    end
    if ~isempty(regexp(app.ExperimentDrop.Value, 'FRAqo'))
        ncyc = 153;
    end
    if ~isempty(regexpi(app.ExperimentDrop.Value, 'rln4db'))
        nldsrepeats = nldsrepeats*21;
    end
    period = hitime + lotime;
    if t1 == 0
        t1 = ncyc*period;
    end
    cycleperiod = t1+gap;
    stimons = 0:period:(t1 - period);
    stimoffs = stimons + hitime;
    timevec = reshape([stimons; stimoffs], 1, numel([stimons; stimoffs]));
    timevec2 = reshape([timevec; timevec], numel(timevec)*2, 1);
    
    trainons = (0:(nrepeats-1)).*cycleperiod;
    trainoffs = trainons + t1;
    
    envelope = reshape([trainons; trainoffs], 1, numel([trainons; trainoffs]));
    timevec3 = repmat(timevec2, 1, nrepeats)+ cycleperiod .* ((1:nrepeats)-1);
    timevec3
    LDSon = (t1+gap)*nrepeats;
    phase1 = reshape(timevec3, 1, numel(timevec3));
    phase2 = repmat([0 0 t2 t2], nldsrepeats, 1) + (0:(nldsrepeats-1))'.*(t2+ldsgap);
    timevec4 = [phase1, reshape(phase2', 1, numel(phase2))+LDSon];
    if ~isempty(regexp(app.ExperimentDrop.Value, 'LDS'))
        timevec4 = [timevec4, phase1 + max(timevec4) + gap2];
    end
    yval = repmat([0 1 1 0], 1, numel(timevec4)/4);

    plot(app.StimAxes, timevec4, yval, 'k-')
    
        
    maxdur = max(timevec4);
    if app.LastoffsetnsecButton.Value
        maxdur = maxdur + app.ExpDurationEdit.Value;
    end
    if app.SecondsButton.Value
        maxdur = app.ExpDurationEdit.Value;
    end
    xlim(app.StimAxes, [0 maxdur]);
    app.StimPanel.Visible = 1;
    app.ExpDurationRadio.Visible = 1;
     app.StimDurLabel.Text = datestr(seconds(maxdur), 'MM:SS');
     for par = 1:numel(app.pars)
        mvarname = ['Param' num2str(par) 'Edit'];
        mvarname2 = ['Param' num2str(par) 'Label'];
        app.(mvarname).Value= str2double(app.(mvarname2).Text);
     end
     for par = 1:numel(app.pars)
        paramupdate(app, par)
     end

end