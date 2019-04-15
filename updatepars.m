%Reads parameters stored in synapse        
function updatepars(app) 
    if app.syn.getMode == 0
        app.syn.setMode(1);
    end
    app.IsoStimulusPanel.Visible = 0;
    app.NBNStimulusPanel.Visible = 0;
    app.TONEStimulusPanel.Visible = 0;
    app.FRAStimulusPanel.Visible = 0;
    app.EVOStimulusPanel.Visible = 0;

    switch app.experiment
        case 1
            app.pars = {'StimHiTime'; 'StimLoTime'; 'StimNPulses'; 'GateLoTime'; 'GateNPulses'; 'Level'};
            app.gizmos = {'xStimulus'; 'xStimulus'; 'xStimulus'; 'xStimulus'; 'xStimulus'; 'xStimulus'};
            prefix = 'Iso';
        case 2
            app.pars = {'TotLDSPulses'; 'LDSOnsetDelay'; 'LDSDuration'; 'LDSTimeBtwnPulses'; 'Level'; 'HPCF'; 'LPCF'; ...
                'ModPercentage'; 'ModFrequency'; 'SineWaveExponent'};
            app.gizmos = {'LDS_gauss'; 'LDS_gauss'; 'LDS_gauss'; 'LDS_gauss'; 'LDS_gauss'; 'LDS_gauss'; 'LDS_gauss'; ...
                'LDS_gauss'; 'LDS_gauss'; 'LDS_gauss'};
            prefix = 'NBN';
        case 3
            app.pars = {'TotLDSPulses'; 'LDSOnsetDelay'; 'LDSDuration'; 'LDSTimeBtwnPulses'; 'Level'; 'ToneFreq'; ...
                'ModPercentage'; 'ModFrequency'; 'SineWaveExponent'};
            app.gizmos = {'LDS_tones'; 'LDS_tones'; 'LDS_tones'; 'LDS_tones'; 'LDS_tones'; 'LDS_tones'; 'LDS_tones'; ...
                'LDS_tones'; 'LDS_tones'};
            prefix = 'TONE';
        case 4
            app.pars = {'GateHiTime'; 'GateLoTime'; 'GateNPulses'};
            app.gizmos = {'xFRA'; 'xFRA'; 'xFRA'};
            prefix = 'FRA';
        case 5
            app.pars = {'EvokeStimHi'; 'EvokeStimLo'; 'EvokeNPulses'; 'EvokeFrequency'; 'EvokeLevel'; 'LDSDuration';...
                'LDSFrequency'; 'LDSLevel'; 'NCycleRepeats'; 'ModulationFrequency'; 'ModulationDepth'};
            app.gizmos = {'xLDS'; 'xLDS'; 'xLDS'; 'xLDS'; 'xLDS'; 'xLDS'; 'xLDS'; 'xLDS'; 'xLDS'; 'xLDS'; 'xLDS'};
            prefix = 'EVO';
        case 6
            app.pars = {'EvokeStimHi'; 'EvokeStimLo'; 'PreLDSDuration'; 'EvokeFrequency'; 'EvokeLevel'; 'LDSDuration';...
                'LDSFrequency'; 'LDSLevel'; 'ModulationDepth'; 'ModulationFrequency'; 'SineExponent'};
            app.gizmos = {'LDSe'; 'LDSe'; 'LDSe'; 'LDSe'; 'LDSe'; 'LDSe'; 'LDSe'; 'LDSe'; 'LDSe'; 'LDSe'; 'LDSe'};
            prefix = 'LE1';
    end
    for par = 1:numel(app.pars)
        mvarname = [prefix char(app.pars(par)) 'Label'];
        app.(mvarname).Text = num2str(app.syn.getParameterValue(char(app.gizmos(par)), char(app.pars(par))));
    end
    switch app.experiment
        case 1
            hitime = str2double(app.IsoStimHiTimeLabel.Text);
            lotime = str2double(app.IsoStimLoTimeLabel.Text);
            period = hitime + lotime;
            stimons = 0:period:period*(str2double(app.IsoStimNPulsesLabel.Text)-1);
            stimoffs = stimons + hitime;
            trainperiod =  stimoffs(end)+ lotime + str2double(app.IsoGateLoTimeLabel.Text);
            train = reshape([stimons; stimoffs], 1, numel([stimons; stimoffs]));
            timemat = repmat(train, str2double(app.IsoGateNPulsesLabel.Text), 1) + (0:trainperiod:trainperiod*(str2double(app.IsoGateNPulsesLabel.Text)-1))';
            timevec = reshape(timemat', 1, numel(timemat));
            timevec2 = reshape([timevec; timevec], 1, numel(timevec)*2);
            yval = repmat([0 1 1 0], 1, numel(timevec2)/4);
            plot(app.IsoStimAxes, timevec2/1000, yval, 'k-')
            maxdur = max(timevec2/1000);
        case 2
            period = str2double(app.NBNLDSDurationLabel.Text) + str2double(app.NBNLDSTimeBtwnPulsesLabel.Text);
            del = str2double(app.NBNLDSOnsetDelayLabel.Text);
            stimons = del:period:(del + period*(str2double(app.NBNTotLDSPulsesLabel.Text))-1);
            stimoffs = stimons + str2double(app.NBNLDSDurationLabel.Text);
            timevec = reshape([stimons; stimoffs], 1, numel([stimons; stimoffs]));
            timevec2 = reshape([timevec; timevec], 1, numel(timevec)*2);
            yval = repmat([0 1 1 0], 1, numel(timevec2)/4);
            plot(app.NBNStimAxes, timevec2, yval, 'k-')
            maxdur = max(timevec2);
        case 3
            period = str2double(app.TONELDSDurationLabel.Text) + str2double(app.TONELDSTimeBtwnPulsesLabel.Text);
            del = str2double(app.TONELDSOnsetDelayLabel.Text);
            stimons = del:period:(del + period*(str2double(app.TONETotLDSPulsesLabel.Text))-1);
            stimoffs = stimons + str2double(app.TONELDSDurationLabel.Text);
            timevec = reshape([stimons; stimoffs], 1, numel([stimons; stimoffs]));
            timevec2 = reshape([timevec; timevec], 1, numel(timevec)*2);
            yval = repmat([0 1 1 0], 1, numel(timevec2)/4);
            plot(app.TONEStimAxes, timevec2/60, yval, 'k-')
            maxdur = max(timevec2);
        case 4
            period = str2double(app.FRAGateHiTimeLabel.Text) + str2double(app.FRAGateLoTimeLabel.Text);
            stimons = 0:period:(period*(str2double(app.FRAGateNPulsesLabel.Text))-1);
            stimoffs = stimons + str2double(app.FRAGateHiTimeLabel.Text);
            timevec = reshape([stimons; stimoffs], 1, numel([stimons; stimoffs]));
            timevec2 = reshape([timevec; timevec], 1, numel(timevec)*2)./1000;
            yval = repmat([0 1 1 0], 1, numel(timevec2)/4);
            plot(app.FRAStimAxes, timevec2/60, yval, 'k-')
            maxdur = max(timevec2);
        case 5
            period = str2double(app.EVOEvokeStimHiLabel.Text) + str2double(app.EVOEvokeStimLoLabel.Text);
            stimons = 0:period:(period*(str2double(app.EVOEvokeNPulsesLabel.Text))-1);
            stimoffs = stimons + str2double(app.EVOEvokeStimHiLabel.Text);
            timevec = reshape([stimons; stimoffs], 1, numel([stimons; stimoffs]));
            t1 = period*str2double(app.EVOEvokeNPulsesLabel.Text)./1000;
            t2 = str2double(app.EVOLDSDurationLabel.Text);
            nrepeats = str2double(app.EVONCycleRepeatsLabel.Text);
            timevec2 = reshape([timevec; timevec], numel(timevec)*2, 1)./1000;
            timevec3 = repmat([timevec2; t1; t1; t1+t2; t1+t2], 1, nrepeats);
            cycleperiod = max(max(timevec3));
            timevec4 = [reshape(timevec3 + cycleperiod .* ((1:nrepeats)-1), 1, numel(timevec3)), timevec2' + cycleperiod*nrepeats];
            yval = repmat([0 1 1 0], 1, numel(timevec4)/4);
            plot(app.EVOStimAxes, timevec4, yval, 'k-')
            maxdur = max(timevec4);
        case 6
            period = str2double(app.LE1EvokeStimHiLabel.Text)/1000 + str2double(app.LE1EvokeStimLoLabel.Text)/1000;
            stimons = 0:period:str2double(app.LE1PreLDSDurationLabel.Text);
            stimoffs = stimons + str2double(app.LE1EvokeStimHiLabel.Text)/1000;
            timevec = reshape([stimons; stimoffs], 1, numel([stimons; stimoffs]));
            t1 = str2double(app.LE1PreLDSDurationLabel.Text);
            t2 = str2double(app.LE1LDSDurationLabel.Text);
            nrepeats = 1;
            timevec2 = reshape([timevec; timevec], numel(timevec)*2, 1);
            timevec3 = repmat([timevec2; t1; t1; t1+t2; t1+t2], 1, nrepeats);
            cycleperiod = max(max(timevec3));
            timevec4 = [reshape(timevec3 + cycleperiod .* ((1:nrepeats)-1), 1, numel(timevec3)), timevec2' + cycleperiod*nrepeats];
            yval = repmat([0 1 1 0], 1, numel(timevec4)/4);
            plot(app.LE1StimAxes, timevec4, yval, 'k-')
            maxdur = max(timevec4);
    end
    stimpanel = [prefix 'StimulusPanel'];
    app.(stimpanel).Visible = 1;
    app.ExpDurationRadio.Visible = 1;
    app.StimDurLabel.Text = datestr(seconds(maxdur), 'MM:SS');
end