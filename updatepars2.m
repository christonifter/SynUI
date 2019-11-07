%Reads parameters stored in synapse        
function updatepars2(app) 
    %Reset things
    if app.syn.getMode == 0
        app.syn.setMode(1);
    end
    for i = 1:12
        synvaluename = ['Param' num2str(i) 'Label'];
        editname = ['Param' num2str(i) 'Edit'];
        paramname = ['Par' num2str(i) 'Label'];
        app.(synvaluename).Visible = 0;
        app.(paramname).Visible = 0;
        app.(editname).Visible = 0;
    end
    app.CenterFreqEditFieldLabel.Visible = 0;
    app.CenterFreqEdit.Visible = 0;
    app.BandwidthEditFieldLabel.Visible = 0;
    app.BandwidthEdit.Visible = 0;
    app.TRMSEdit.Visible = 0;

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
    if sum(ismember(app.pars, 'HPCF'))
        app.CenterFreqEdit.Visible = 1;
        app.BandwidthEdit.Visible = 1;
        app.CenterFreqEditFieldLabel.Visible = 1;
        app.BandwidthEditFieldLabel.Visible = 1;
    end
    if sum(ismember(app.pars, 'ModDepth'))
        app.TRMSEdit.Visible = 1;
    end
        
    
    % build stim preview    
    hitime = 0;
    lotime = 0;
    ncyc = 1;
    nrepeats = 1;
    t1 = 0;
    t2 = 0;
    gap = 0;
    gap2 = 0;
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
            case 'StimNPulses'
                ncyc = parval;
            case 'TrainDurationMS'
                t1 = parval/1000;
            case 'GateLoMS'
                gap = parval/1000;
            case 'LDSISISec'
                gap = parval;
            case 'TrainGapMS'
                gap = parval/1000;
            case 'PostLDSGapMS'
                gap2 = parval/1000;
            case 'NumberofTrains'
                nrepeats = parval;
            case 'GateNPulses'
                nrepeats = parval;
            case 'LDSPreSEC' 
                t1 = parval;
            case 'LDSDurSEC'
                t2 = parval;
            case 'LDSDurationMS'
                t2 = parval/1000;
            case 'LDSNPulses'
                nrepeats = parval;
        end
    end
    period = hitime + lotime;
    if t1 == 0
        t1 = ncyc*period;
    end
    stimons = 0:period:(t1 - period);
    stimoffs = stimons + hitime;
    timevec = reshape([stimons; stimoffs], 1, numel([stimons; stimoffs]));
    timevec2 = reshape([timevec; timevec], numel(timevec)*2, 1);
    timevec3 = repmat(timevec2, 1, nrepeats);
    cycleperiod = max(max(timevec3))+gap;
    LDSon = (t1+gap)*nrepeats;
    LDSoff = LDSon+t2;
    timevec3b = [reshape(timevec3 + cycleperiod .* ((1:nrepeats)-1), 1, numel(timevec3)), LDSon, LDSon, LDSoff, LDSoff];
    timevec4 = [timevec3b, LDSoff + reshape(timevec3 + gap2 + cycleperiod .* ((1:nrepeats)-1), 1, numel(timevec3))];
    
%     timevec3 = repmat([timevec2; t1; t1; t1+t2; t1+t2], 1, nrepeats);
%     cycleperiod = max(max(timevec3))+gap;
%     timevec4 = [reshape(timevec3 + cycleperiod .* ((1:nrepeats)-1), 1, numel(timevec3))];
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

end