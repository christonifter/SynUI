function rms2amp2(app)
    dbshift = [0, -1.765, -.337, .466, 1.111, 1.535, 1.981, 2.246, 2.589]; %correction for modulation exponent

    HPi = find(ismember(app.pars, 'HPCF'));
    HPfield = ['Param' num2str(HPi) 'Edit'];
    LPi = find(ismember(app.pars, 'LPCF'));
    LPfield = ['Param' num2str(LPi) 'Edit'];
    BW =  (app.(LPfield).Value - app.(HPfield).Value)./1000;
     scaling = 18/(BW+3.75)+0.6;
     dbshift2 = 20*log10(scaling); %correction for bandwidth
%     dbshift2 = 12.5*BW;
        
    MEi = find(ismember(app.pars, 'ModExponent'));
    MEfield = ['Param' num2str(MEi) 'Edit'];
    MDi = find(ismember(app.pars, 'ModDepth'));
    MDfield = ['Param' num2str(MDi) 'Edit'];
    exp = round(app.(MDfield).Value) * app.(MEfield).Value+1;
    
    if sum(ismember(app.pars, 'LDSAmplitude'))
        Ampi = find(ismember(app.pars, 'LDSAmplitude'));
        Ampfield = ['Param' num2str(Ampi) 'Edit'];
        app.(Ampfield).Value = 10^((dbshift(exp) + dbshift2 + app.LDSLevelEdit.Value)/20);
    else
        Shifti = find(ismember(app.pars, 'dBShift'));
        Shiftfield = ['Param' num2str(Shifti) 'Edit'];
        app.(Shiftfield).Value = dbshift(exp) + dbshift2;
    end
    EAmpi = find(ismember(app.pars, 'StimAmplitude'));
    EAmpfield = ['Param' num2str(EAmpi) 'Edit'];
    if sum(ismember(app.pars, 'StimModFrequencyHz'))
        
            HPi = find(ismember(app.pars, 'StimHPCF'));
            HPfield = ['Param' num2str(HPi) 'Edit'];
            LPi = find(ismember(app.pars, 'StimLPCF'));
            LPfield = ['Param' num2str(LPi) 'Edit'];
            BW =  (app.(LPfield).Value - app.(HPfield).Value)./1000;
             scaling = 18/(BW+3.75)+0.6;
             dbshift2 = 20*log10(scaling);
        
        app.(EAmpfield).Value = 10^((dbshift(2) + dbshift2 + app.StimLevelEdit.Value)/20);
    elseif sum(ismember(app.pars, 'StimAmplitude'))
        app.(EAmpfield).Value = 10^(app.StimLevelEdit.Value/20);
    end
    
end