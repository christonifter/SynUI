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
    
    app.StimulusPanel.Visible = 0;
    app.LDSPanel.Visible = 0;
    gizlist = app.syn.getGizmoNames();
    c = regexp(gizlist, '^x');
    pargizs = gizlist(~cellfun('isempty', c));
    c = 1;
    app.pars = {};
    app.gizmos = {};
    
    
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
    
    
    hitime = 0;
    lotime = 0;
    ncyc = 1;
    nrepeats = 0;
    nldsrepeats = 1;
    t1 = 0;
    t2 = 0;
    gap = 0;
    gap2 = 0;
    ldsgap = 0;

    
    
    %Hide things
    if sum([ismember(app.pars, 'StimHiMS'); ismember(app.pars, 'TrainDurationMS')])
        app.StimulusPanel.Visible = 1;
        showstim = sum(ismember(app.pars, 'StimHiMS'));
            app.PulseDurEdit.Visible = showstim;
            app.PulseDurationmsLabel.Visible = showstim;
            app.PulseRateEdit.Visible = showstim;
            app.PulserateHzLabel.Visible = showstim;
            app.RampEdit.Visible = showstim;
            app.RiseFallmsLabel.Visible = showstim;
        showtrain = sum(ismember(app.pars, 'TrainDurationMS')); 
            app.TrainOnEdit.Visible = showtrain;
            app.TrainOnLabel.Visible = showtrain;
            app.TrainOffEdit.Visible = showtrain;
            app.TrainOffLabel.Visible = showtrain;

        if sum(ismember(app.pars, 'StimFrequencyHz'))
            app.StimFreqEdit.Visible = 1;
            app.StimFreqLabel.Visible = 1;
            app.StimBWEdit.Visible = 0;
            app.StimBWLabel.Visible = 0;
        elseif sum(ismember(app.pars, 'StimHPCF')) || sum(ismember(app.pars, 'StimOctave1'))
            app.StimFreqEdit.Visible = 1;
            app.StimFreqLabel.Visible = 1;
            app.StimBWEdit.Visible = 1;
            app.StimBWLabel.Visible = 1;
        else
            app.StimFreqEdit.Visible = 0;
            app.StimFreqLabel.Visible = 0;
            app.StimBWEdit.Visible = 0;
            app.StimBWLabel.Visible = 0;
        end
        
        if sum(ismember(app.pars, 'StimOctave1'))
            app.PulseDurEdit.Visible = 0;
            app.PulseDurationmsLabel.Visible = 0;
            sweepspeed = 1.5; % oct/ms
            
            Oct1 = app.syn.getParameterValue('xAud', 'StimOctave1');
            Oct2 = app.syn.getParameterValue('xAud', 'StimOctave2');
            RampMS = app.syn.getParameterValue('xAud', 'RiseFallMS');
            BW = Oct2 - Oct1;
            StimHiMS = BW * sweepspeed - RampMS * 1.25;
            app.syn.setParameterValue('xLDS', 'StimHiMS', StimHiMS);
            hitime = BW * sweepspeed/1E3;
        end
        if sum(ismember(app.pars, 'StimAmplitude'))
            app.StimLevelEdit.Visible = 1;
            app.StimLevelLabel.Visible = 1;
        else
            app.StimLevelEdit.Visible = 0;
            app.StimLevelLabel.Visible = 0;
        end
        if sum(ismember(app.pars, 'StimModFrequencyHz'))
            app.StimModFreqEdit.Visible = 1;
            app.StimModFreqLabel.Visible = 1;
        else
            app.StimModFreqEdit.Visible = 0;
            app.StimModFreqLabel.Visible = 0;
        end
        if sum(ismember(app.pars, 'Stim2FrequencyHz'))
            app.Stim2FreqEdit.Visible = 1;
            app.Stim2FreqLabel.Visible = 1;
        else
            app.Stim2FreqEdit.Visible = 0;
            app.Stim2FreqLabel.Visible = 0;
        end
        if sum(ismember(app.pars, 'Stim2Amplitude'))
            app.Stim2LevelEdit.Visible = 1;
            app.Stim2LevelLabel.Visible = 1;
        else
            app.Stim2LevelEdit.Visible = 0;
            app.Stim2LevelLabel.Visible = 0;
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
            app.LDSPostGapEdit.Visible = 0;
            app.LDSISILabel.Visible = 1;
            app.LDSPostGapLabel.Visible = 0;
        else
            app.LDSISIEdit.Visible = 0;
            app.LDSRepeatsEdit.Visible = 0;
            app.LDSRepeatsLabel.Visible = 0;
            app.LDSPostGapEdit.Visible = 1;
            app.LDSISILabel.Visible = 0;
            app.LDSPostGapLabel.Visible = 1;
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
        app.LDSLevelEdit.Visible = 1;
        app.LDSLevelLabel.Visible = 1;
    else
        app.LDSLevelEdit.Visible = 0;
        app.LDSLevelLabel.Visible = 0;
    end

    if sum(ismember(app.pars, 'buffer2'))
        app.DialogueLabel.Text = 'PC -> Audio -> RZ6 (~1 min)';
        pause(.1)
        bufferduration = 2000; %in ms, use an even number of seconds.
        bufferlength = round(bufferduration*195.3125);
        fid = fopen('C:/TDT/Synapse/DMR500ic120s30dB.f32');
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
    timevec2a = [0; 0; timevec2(3:(end-2), 1); t1; t1];
    timevec3 = repmat(timevec2a, 1, nrepeats)+ cycleperiod .* ((1:nrepeats)-1);
     
    LDSon = (t1+gap)*max([nrepeats 1]);
    phase1 = reshape(timevec3, 1, numel(timevec3));
    phase2 = repmat([0 0 t2 t2], nldsrepeats, 1) + (0:(nldsrepeats-1))'.*(t2+ldsgap);
    timevec4 = [phase1, reshape(phase2', 1, numel(phase2))+LDSon];
    if ~isempty(regexp(app.ExperimentDrop.Value, 'LDS'))
        timevec4 = [timevec4, phase1 + max(timevec4) + gap2];
    end
%     if strcmp(app.ExperimentDrop.Value(1), 'h')
%         timevec4 = 10+[timevec4, phase1 + max(timevec4) + gap2];
%     end
    yval = repmat([0 1 1 0], 1, numel(timevec4)/4);

    plot(app.StimAxes, timevec4, yval, 'k-')
    
        
    maxdur = max([1 timevec4]);
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
        make_pars_userfriendly(app, par)
     end

end