function findrms(app)
    dbshift = [0, -1.765, -.337, .466, 1.111, 1.535, 1.981, 2.246, 2.589];
    MEi = find(ismember(app.pars, 'ModExponent'));
    MEfield = ['Param' num2str(MEi) 'Edit'];
    find(ismember(app.pars, 'ModDepth'))
    MDi = find(ismember(app.pars, 'ModDepth'));
    MDfield = ['Param' num2str(MDi) 'Edit'];
    Leveli = find(ismember(app.pars, 'LDSamplitude'));
    Levelfield = ['Param' num2str(Leveli) 'Edit'];
    exp = round(app.(MDfield).Value/100) * app.(MEfield).Value+1;
    app.(Levelfield).Value = dbshift(exp) + app.TRMSEdit.Value;
end