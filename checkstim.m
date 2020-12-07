function checkstim(app)
if ~isempty(app.data.stimsig)
    switch lower(app.PSTH1RefDrop.Value)
        case '0'
            offset = 0;
        case 'onset'
            offset = app.data.stimons(1);
        case 'offset'
            offset = app.data.stimoffs(1);
    end
    
    t1 = 1+round((app.PSTH1StartEdit.Value+offset)*app.data.stimfs);
    t2 = 1+round((app.PSTH1EndEdit.Value+offset)*app.data.stimfs);
    
    switch lower(app.PSTH2RefDrop.Value)
        case '0'
            offset = 0;
        case 'onset'
            offset = app.data.stimons(1);
        case 'offset'
            offset = app.data.stimoffs(1);
    end
    
    t3 = 1+round((app.PSTH2StartEdit.Value+offset)*app.data.stimfs);
    t4 = 1+round((app.PSTH2EndEdit.Value+offset)*app.data.stimfs);
    AmpRMS1 = rms(app.data.stimsig(t1:t2));
    AmpRMS2 = rms(app.data.stimsig(t3:t4));
    AmpPP1 = range(app.data.stimsig(t1:t2));
    AmpPP2 = range(app.data.stimsig(t3:t4));
    SynStimAmp = NaN;
    SynLDSAmp = NaN;
    SynStimLevel = NaN;
    SynLDSLevel = NaN;
    SynStimFreq = NaN;
    SynLDSCFreq = NaN;
    load([app.TankEdit.Value '\params.mat']');
    if ismember('StimAmplitude', changetable.Properties.VariableNames)
        SynStimAmp = changetable.StimAmplitude(1);
    end
    if ismember('LDSAmplitude', changetable.Properties.VariableNames)
        SynLDSAmp = changetable.LDSAmplitude(1);
    end
    if sum(strcmpi(stimtable.stimname, 'Stimlevel'))
        SynStimLevel = stimtable.stimval(find(strcmpi(stimtable.stimname, 'Stimlevel')));
    end
    if sum(strcmpi(stimtable.stimname, 'LDSlevel'))
        SynLDSLevel = stimtable.stimval(find(strcmpi(stimtable.stimname, 'LDSlevel')));
    end
    if sum(strcmpi(stimtable.stimname, 'StimFreq'))
        SynStimFreq = stimtable.stimval(find(strcmpi(stimtable.stimname, 'StimFreq')));
    end
    if sum(strcmpi(stimtable.stimname, 'LDSCFreq'))
        SynLDSCFreq = stimtable.stimval(find(strcmpi(stimtable.stimname, 'LDSCFreq')));
    end
    try
        load([app.TankEdit.Value '\params.mat']', 'changetable', 'stimtable');
    end
    app.data.stimstable = table(SynStimFreq, SynStimLevel, SynStimAmp,  AmpRMS1, AmpPP1, SynLDSCFreq, SynLDSLevel, SynLDSAmp, AmpRMS2, AmpPP2);
end
end