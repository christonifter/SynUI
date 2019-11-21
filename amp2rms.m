function amp2rms(app)
    dbshift = [0, -1.765, -.337, .466, 1.111, 1.535, 1.981, 2.246, 2.589];
    MEi = find(ismember(app.pars, 'ModExponent'));
    MEfield = ['Param' num2str(MEi) 'Edit'];
    MDi = find(ismember(app.pars, 'ModDepth'));
    MDfield = ['Param' num2str(MDi) 'Edit'];
    exp = round(app.(MDfield).Value) * app.(MEfield).Value+1;
    if sum(ismember(app.pars, 'LDSAmplitude'))
        Ampi = find(ismember(app.pars, 'LDSAmplitude'));
        Ampfield = ['Param' num2str(Ampi) 'Edit'];
        app.TRMSEdit.Value = 20*log10(app.(Ampfield).Value) - dbshift(exp);
    end

end