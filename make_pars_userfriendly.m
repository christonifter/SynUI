function paramupdate(app, ind)
    LPCFi = find(ismember(app.pars, 'LPCF'));
    HPCFi = find(ismember(app.pars, 'HPCF'));
    LPCFfield = ['Param' num2str(LPCFi) 'Edit'];
    HPCFfield = ['Param' num2str(HPCFi) 'Edit'];
    Oct1i = find(ismember(app.pars, 'StimOctave1'));
    Oct2i = find(ismember(app.pars, 'StimOctave2'));
    Oct1field = ['Param' num2str(Oct1i) 'Edit'];
    Oct2field = ['Param' num2str(Oct2i) 'Edit'];
    StimLPCFi = find(ismember(app.pars, 'StimLPCF'));
    StimHPCFi = find(ismember(app.pars, 'StimHPCF'));
    StimLPCFfield = ['Param' num2str(StimLPCFi) 'Edit'];
    StimHPCFfield = ['Param' num2str(StimHPCFi) 'Edit'];
    StimHii = find(ismember(app.pars, 'StimHiMS'));
    StimLoi = find(ismember(app.pars, 'StimLoMS'));
    StimHifield = ['Param' num2str(StimHii) 'Edit'];
    StimLofield = ['Param' num2str(StimLoi) 'Edit'];
    parfield = ['Param' num2str(ind) 'Edit'];

    switch char(app.pars(ind))
        case 'StimHiMS'
            app.PulseDurEdit.Value = app.(StimHifield).Value + app.RampEdit.Value*1.25;
            app.PulseRateEdit.Value = 1E3/(app.(StimHifield).Value + app.(StimLofield).Value);
        case 'StimLoMS'
            app.PulseRateEdit.Value = 1E3/(app.(StimHifield).Value + app.(StimLofield).Value);
        case 'RiseFallMS'
            app.RampEdit.Value = app.(parfield).Value;
            app.PulseDurEdit.Value = app.(StimHifield).Value + app.RampEdit.Value*1.25;
        case 'TrainDurationMS'
            app.TrainOnEdit.Value = round(app.(parfield).Value)/1000;
        case 'TrainGapMS'
            app.TrainOffEdit.Value = round(app.(parfield).Value)/1000;
        case 'TrainRepeats'
            app.TrainRepeatsEdit.Value = app.(parfield).Value;
        case 'LDSPreMS'
            app.LDSPreEdit.Value = round(app.(parfield).Value)/1000;
        case 'LDSDurationMS'
            app.LDSDurEdit.Value = round(app.(parfield).Value)/1000;
        case 'LDSISIMS'
            app.LDSISIEdit.Value = round(app.(parfield).Value)/1000;
        case 'LDSRepeats'
            app.LDSRepeatsEdit.Value = round(app.(parfield).Value);
        case 'PostLDSGapMS'
            app.LDSPostGapEdit.Value = round(app.(parfield).Value)/1000;
        case 'StimAmplitude'
            app.StimLevelEdit.Value = round(20*log10(app.(parfield).Value)*1E3)/1E3;
        case 'StimFrequencyHz'
            app.StimFreqEdit.Value = round(app.(parfield).Value)/1000;
        case 'Stim2Amplitude'
            app.Stim2LevelEdit.Value = round(20*log10(app.(parfield).Value)*1E3)/1E3;
        case 'Stim2FrequencyHz'
            app.Stim2FreqEdit.Value = round(app.(parfield).Value)/1000;
        case 'LDSAmplitude'
             app.LDSLevelEdit.Value = round(20*log10(app.(parfield).Value)*1E3)/1E3;
            amp2rms2(app)
        case 'ModFrequencyHz'
            app.LDSModFreqEdit.Value = app.(parfield).Value;
        case 'HPCF'
            app.LDSCFreqEdit.Value = round(sqrt(app.(HPCFfield).Value * app.(LPCFfield).Value))/1E3;
            app.LDSBWEdit.Value = round(log2(app.(LPCFfield).Value / app.(HPCFfield).Value)*1E3)/1E3;  
        case 'LPCF'
            app.LDSCFreqEdit.Value = round(sqrt(app.(HPCFfield).Value * app.(LPCFfield).Value))/1E3;
            app.LDSBWEdit.Value = round(log2(app.(LPCFfield).Value / app.(HPCFfield).Value)*1E3)/1E3;  
        case 'StimOctave1'
            app.StimFreqEdit.Value = round(2^(0.5*(app.(Oct1field).Value + app.(Oct2field).Value)))/1E3;
            app.StimBWEdit.Value = app.(Oct2field).Value - app.(Oct1field).Value;  
        case 'StimOctave2'
            app.StimFreqEdit.Value = round(2^(0.5*(app.(Oct1field).Value + app.(Oct2field).Value)))/1E3;
            app.StimBWEdit.Value = app.(Oct2field).Value - app.(Oct1field).Value;
        case 'StimHPCF'
            app.StimFreqEdit.Value = round(sqrt(app.(StimHPCFfield).Value * app.(StimLPCFfield).Value))/1E3;
            app.StimBWEdit.Value = round(log2(app.(StimLPCFfield).Value / app.(StimHPCFfield).Value)*1E3)/1E3;  
        case 'StimLPCF'
            app.StimFreqEdit.Value = round(sqrt(app.(StimHPCFfield).Value * app.(StimLPCFfield).Value))/1E3;
            app.StimBWEdit.Value = round(log2(app.(StimLPCFfield).Value / app.(StimHPCFfield).Value)*1E3)/1E3;  
        case 'ModDepth'
            app.LDSModDepthEdit.Value = app.(parfield).Value*100;
%             amp2rms2(app)
        case 'ModExponent'
            app.LDSModExpEdit.Value = app.(parfield).Value ;
%             amp2rms2(app)
        case 'StimModFrequencyHz'
            app.StimModFreqEdit.Value = app.(parfield).Value;
    end
end