function out = readtankfast(path)
    data=TDTbin2mat(path);
    spikevar = 'eSpk';
    wavevar = 'RAW1';
    stimvar = 'Freq';
    stimvar2 = 'Stmb';
    stimvar3 = 'Stim';
    offvar = 'StOn';
    evvar = 'EvOn';
    scalpvar = 'RA4W';
    stimsigvar = 'Sti1';
    
    out.frqs = 0;
    out.stimons = data.epocs.(offvar).onset(data.epocs.(offvar).data == 1);
    out.stimoffs = data.epocs.(offvar).onset(data.epocs.(offvar).data == 0);
    out.evons = [];
    out.evoffs = [];
    out.stimsig = [];
    out.stimfs = [];
    if isfield(data.epocs, evvar)
        evons = data.epocs.(evvar).onset(data.epocs.(evvar).data == 1);
        out.evoffs = data.epocs.(evvar).onset(data.epocs.(evvar).data == 0);
        out.evons = evons(1:numel(out.evoffs));
    end
    try
        out.stimsig = data.streams.(stimsigvar).data;
        out.stimfs = data.streams.(stimsigvar).fs;
    end
    try
        out.stimsig = data.streams.(stimvar3).data.*1E-3; %LDS-SN measured in mV
        out.stimfs = data.streams.(stimvar3).fs;
    end
    
if isfield(data.streams, scalpvar)
    out.raw = data.streams.(scalpvar).data'./20; %in V, 20x gain on the RA4LI
    out.fs = data.streams.(scalpvar).fs;
    
    out.evons = out.evons + 2.05*1E-3; %RA4 headstage delay
    out.evoffs = out.evoffs + 2.05*1E-3;
    out.stimons = out.stimons + 2.05*1E-3;
    out.stimoffs = out.stimoffs + 2.05*1E-3;
else
%           out.spets = data.snips.(spikevar).ts+.3; %FRAs collected before 1/7/2019 skipped the first stimulation
    out.spets = data.snips.(spikevar).ts;
    out.channels = double(data.snips.(spikevar).chan);
    clusters = double(data.snips.(spikevar).sortcode);
    clusterlist = sort(unique(clusters));
    for clusteri = 1:numel(clusterlist)
        clusters(clusters == clusterlist(clusteri)) = clusteri;
    end
    out.clusters = double(clusters);
    out.stimons = data.epocs.(offvar).onset(data.epocs.(offvar).data == 1);
    out.stimoffs = data.epocs.(offvar).onset(data.epocs.(offvar).data == 0);
    out.lvls = zeros(size(out.stimons));


    if isfield(data.epocs, stimvar3) %FRA
        out.lvls = data.epocs.(stimvar2).data;
        out.frqs = data.epocs.(stimvar3).data(2:end);
    elseif isfield(data.epocs, stimvar2) %Iso-F
        out.lvls = data.epocs.(stimvar2).data;
        out.frqs = zeros(numel(out.stimons), 1); %RLN
        try
            load([path '/params.mat']');
            if ismember('Frequency', changetable.Properties.VariableNames)
                out.frqs = ones(size(out.stimons)).*changetable.Frequency(1);
            end
        end
    else
        if isfield(data.epocs, stimvar) %Iso-intensity
            out.frqs = data.epocs.(stimvar).data;
        else %LDS
            out.frqs = zeros(numel(out.stimons), 1);
        end
        try
            load([path '/params.mat']');
            if ismember('Level', changetable.Properties.VariableNames)
                out.lvls = ones(size(out.stimons)).*changetable.Level(1);
            end
        end
    end
    out.evons = [];
    out.evoffs = [];
    if isfield(data.epocs, evvar)
        evons = data.epocs.(evvar).onset(data.epocs.(evvar).data == 1);
        out.evoffs = data.epocs.(evvar).onset(data.epocs.(evvar).data == 0);
        out.evons = evons(1:numel(out.evoffs));
        out.frqs = zeros(numel(out.evons), 1); %RLN
        out.lvls = zeros(numel(out.evons), 1); %RLN
        try
            load([path '/params.mat']');
            if ismember('Frequency', changetable.Properties.VariableNames)
                out.frqs = ones(size(out.evons)).*changetable.Frequency(1);
            end
            if ismember('Levels', changetable.Properties.VariableNames)
                out.lvls = ones(size(out.evons)).*changetable.Levels(1);
            end
        end

    end
    if 1 == 0
        fny = data.streams.RAW1.fs/2; 
        fcl = 300;
        fcu = 3000;
        out.fs = data.streams.RAW1.fs;
         out.raw = data.streams.RAW1.data';
%          out.LFPraw = data.streams.LFPw.data';
        out.LFPfs = data.streams.RAW1.fs;
        
        [b,a] = butter(2, [fcl/fny fcu/fny]);
        out.LFP= squeeze(filter(b,a, out.raw));
    else
        out.fs = data.streams.RAW1.fs;
        out.raw = [];
        out.LFPfs = [];
        out.LFP= [];
        
    end
     
end
end