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
    SynLDSMD = NaN;
    SynLDSME = NaN;
    SynLDSBW = NaN;
    try
        load([app.TankEdit.Value '\params.mat']');
        if ismember('StimAmplitude', changetable.Properties.VariableNames)
            SynStimAmp = changetable.StimAmplitude(1);
        end
        if ismember('LDSAmplitude', changetable.Properties.VariableNames)
            SynLDSAmp = changetable.LDSAmplitude(1);
        end
        if ismember('StimLevel', changetable.Properties.VariableNames)
            SynStimLevel = changetable.StimLevel(1);
        end
        if ismember('LDSLevel', changetable.Properties.VariableNames)
            SynLDSLevel = changetable.LDSLevel(1);
        end
        if ismember('StimFrequency', changetable.Properties.VariableNames)
            SynStimFreq = changetable.StimFrequency(1);
        end
        
        if ismember('StimFrequencyHz', changetable.Properties.VariableNames)
            SynStimFreq = changetable.StimFrequencyHz(1);
        end
        if ismember('HPCF', changetable.Properties.VariableNames)
            SynLDSCFreq = sqrt(changetable.HPCF(1) * changetable.LPCF(1));
            SynLDSBW = log2(changetable.LPCF(1)./changetable.HPCF(1));
        end
        if ismember('ModDepth', changetable.Properties.VariableNames)
            SynLDSMD = changetable.ModDepth(1);
        end
        if ismember('ModExponent', changetable.Properties.VariableNames)
            SynLDSME = changetable.ModExponent(1);
        end
        
        if sum(strcmpi(stimtable.stimname, 'Tonelevel'))
            SynStimLevel = stimtable.stimval(find(strcmpi(stimtable.stimname, 'Tonelevel')));
        end
        if sum(strcmpi(stimtable.stimname, 'Stimlevel'))
            SynStimLevel = stimtable.stimval(find(strcmpi(stimtable.stimname, 'Stimlevel')));
        end
        if sum(strcmpi(stimtable.stimname, 'LDSlevel'))
            SynLDSLevel = stimtable.stimval(find(strcmpi(stimtable.stimname, 'LDSlevel')));
        end
        if sum(strcmpi(stimtable.stimname, 'TRMSEdit'))
            SynLDSLevel = stimtable.stimval(find(strcmpi(stimtable.stimname, 'TRMSEdit')));
        end
%         if sum(strcmpi(stimtable.stimname, 'ToneFreq'))
%             SynStimFreq = stimtable.stimval(find(strcmpi(stimtable.stimname, 'ToneFreq'.*1000)));
%         end
%         if sum(strcmpi(stimtable.stimname, 'StimFreq'))
%             SynStimFreq = stimtable.stimval(find(strcmpi(stimtable.stimname, 'StimFreq')));
%         end
%         if sum(strcmpi(stimtable.stimname, 'CenterFreq'))
%             SynLDSCFreq = stimtable.stimval(find(strcmpi(stimtable.stimname, 'CenterFreq'.*1000)));
%         end
%         if sum(strcmpi(stimtable.stimname, 'LDSCFreq'))
%             SynLDSCFreq = stimtable.stimval(find(strcmpi(stimtable.stimname, 'LDSCFreq')));
%         end
%         
%         if sum(strcmpi(stimtable.stimname, 'Bandwidth'))
%             SynLDSBW = stimtable.stimval(find(strcmpi(stimtable.stimname, 'Bandwidth')));
%         end
%         if sum(strcmpi(stimtable.stimname, 'LDSBW'))
%             SynLDSBW = stimtable.stimval(find(strcmpi(stimtable.stimname, 'LDSBW')));
%         end
%         if sum(strcmpi(stimtable.stimname, 'ModDepth'))
%             SynLDSMD = stimtable.stimval(find(strcmpi(stimtable.stimname, 'ModDepth')));
%         end
%         if sum(strcmpi(stimtable.stimname, 'LDSModDepth'))
%             SynLDSMD = stimtable.stimval(find(strcmpi(stimtable.stimname, 'LDSModDepth')));
%         end

        
    end
    app.data.stimstable = table(SynStimFreq, SynStimLevel, SynStimAmp,  AmpRMS1, AmpPP1, SynLDSCFreq, SynLDSLevel, SynLDSAmp, AmpRMS2, AmpPP2, SynLDSBW, SynLDSMD, SynLDSME);
end
end