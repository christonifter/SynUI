function paramupdate(app, ind)
    LPCFi = find(ismember(app.pars, 'LPCF'));
    HPCFi = find(ismember(app.pars, 'HPCF'));
    LPCFfield = ['Param' num2str(LPCFi) 'Edit'];
    HPCFfield = ['Param' num2str(HPCFi) 'Edit'];
    StimLPCFi = find(ismember(app.pars, 'StimLPCF'));
    StimHPCFi = find(ismember(app.pars, 'StimHPCF'));
    StimLPCFfield = ['Param' num2str(StimLPCFi) 'Edit'];
    StimHPCFfield = ['Param' num2str(StimHPCFi) 'Edit'];
    parfield = ['Param' num2str(ind) 'Edit'];

    switch char(app.pars(ind))
        case 'StimHiMS'
            app.StimOnEdit.Value = round(app.(parfield).Value*1E3)/1E3;
        case 'StimLoMS'
            app.StimOffEdit.Value = round(app.(parfield).Value*1E3)/1E3;
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
            app.PostGapEdit.Value = round(app.(parfield).Value)/1000;
        case 'StimAmplitude'
            app.StimLevelEdit.Value = round(20*log10(app.(parfield).Value)*1E3)/1E3;
        case 'StimFrequencyHz'
            app.StimFreqEdit.Value = round(app.(parfield).Value)/1000;
        case 'Stim2Amplitude'
            app.Stim2LevelEdit.Value = round(20*log10(app.(parfield).Value)*1E3)/1E3;
        case 'Stim2FrequencyHz'
            app.Stim2FreqEdit.Value = round(app.(parfield).Value)/1000;
        case 'LDSAmplitude'
%             amp2rms(app)
        case 'ModFrequencyHz'
            app.ModFreqEdit.Value = app.(parfield).Value;
        case 'HPCF'
            app.CenterFreqEdit.Value = round(sqrt(app.(HPCFfield).Value * app.(LPCFfield).Value))/1E3;
            app.BandwidthEdit.Value = round(log2(app.(LPCFfield).Value / app.(HPCFfield).Value)*1E3)/1E3;  
        case 'LPCF'
            app.CenterFreqEdit.Value = round(sqrt(app.(HPCFfield).Value * app.(LPCFfield).Value))/1E3;
            app.BandwidthEdit.Value = round(log2(app.(LPCFfield).Value / app.(HPCFfield).Value)*1E3)/1E3;  
        case 'StimHPCF'
            app.StimFreqEdit.Value = round(sqrt(app.(StimHPCFfield).Value * app.(StimLPCFfield).Value))/1E3;
            app.StimBWEdit.Value = round(log2(app.(StimLPCFfield).Value / app.(StimHPCFfield).Value)*1E3)/1E3;  
        case 'StimLPCF'
            app.StimFreqEdit.Value = round(sqrt(app.(StimHPCFfield).Value * app.(StimLPCFfield).Value))/1E3;
            app.StimBWEdit.Value = round(log2(app.(StimLPCFfield).Value / app.(StimHPCFfield).Value)*1E3)/1E3;  
        case 'ModDepth'
            app.ModDepthEdit.Value = app.(parfield).Value*100;
%             amp2rms(app)
        case 'ModExponent'
            app.ModExpEdit.Value = app.(parfield).Value ;
%             amp2rms(app)
        case 'StimModFrequencyHz'
            app.StimModFreqEdit.Value = app.(parfield).Value;
    end
end