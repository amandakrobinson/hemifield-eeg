function run_preprocessing(partid)

    % preprocess hemifield eeg data:
    % interpolate bad channels
    % average rereference
    % high pass 0.1hz
    % low pass 100hz
    % epoch -200 to 800ms for each image
    % perform CSD/laplacian transformation on data
    % convert data to cosmomvpa format
    
    %% eeglab
    if isempty(which('eeglab'))
        addpath('~/Dropbox (Personal)/MATLAB/eeglab2021.1')
    end
    eeglab
    
    %% cosmomvpa
    if isempty(which('cosmo_wtf'))
         addpath('~/Dropbox/MATLAB/CoSMoMVPA/mvpa')
    end

    addpath('~/Dropbox/MATLAB/CSDtoolbox/func')

    %% get files
    datapath = 'data';
    
    contfn = sprintf('%s/derivatives/eeglab/sub-%02i_task-rsvp_continuous.set',datapath,partid);
    
    % load laplacian montage
    load(sprintf('%s/derivatives/CSD_montage.mat',datapath));
    
    %% run preprocessing
    if isfile(contfn)
        fprintf('Using %s\n',contfn)
    	EEG_cont = pop_loadset(contfn);
    else
        % load EEG file
        EEG_cont = pop_loadbv(sprintf('%s/sub-%02i/eeg/',datapath,partid), sprintf('sub-%02i_task-rsvp_eeg.vhdr',partid));
        EEG_cont = eeg_checkset(EEG_cont);
        EEG_cont.setname = partid;
        EEG_cont = eeg_checkset(EEG_cont);
        
        % get/interpolate bad channels - needed for CSD
        [~, badidx] = pop_rejchan(EEG_cont, 'elec',1:63 ,'threshold',5,'norm','on','measure','kurt');
        EEG_cont = eeg_interp(EEG_cont,badidx);
        EEG_cont = eeg_checkset(EEG_cont);
        
        % average re-reference
        % add Cz channel and re-reference, adding Cz back into dataset
        EEG_cont=pop_chanedit(EEG_cont, 'append',63,'changefield',{64 'labels' 'Cz'},'setref',{'' 'Cz'});
        Czloc = struct('labels',{'Cz'},'type',{''},'theta',{0},'radius',{0},'X',{5.2047e-15},'Y',{0},'Z',{85},'sph_theta',{0},'sph_phi',{90},'sph_radius',{85},'urchan',{64},'ref',{''},'datachan',{0});
        EEG_cont = pop_reref( EEG_cont, [],'refloc',Czloc);
        EEG_cont = eeg_checkset(EEG_cont);
        
        % high pass filter
        EEG_cont = pop_eegfiltnew(EEG_cont, 0.1,[]);

        % low pass filter
        EEG_cont = pop_eegfiltnew(EEG_cont, [],100);

        pop_saveset(EEG_cont,contfn);
    end

    %% add eventinfo to events
    eventsfncsv = sprintf('%s/sub-%02i/eeg/sub-%02i_task-rsvp_events.csv',datapath,partid,partid);
    eventsfntsv = sprintf('%s/sub-%02i/eeg/sub-%02i_task-rsvp_events.tsv',datapath,partid,partid);
    eventlist = readtable(eventsfncsv);

    % add in onset and duration
    idx = strcmp({EEG_cont.event.type},'E  1');
    onset = vertcat(EEG_cont.event(idx).latency);
    duration = 100*ones(size(onset));

    neweventlist = [table(onset,duration,'VariableNames',{'onset','duration'}) eventlist];

    % edit eventlist to change "italics" trials to non-italics
    % order of images is 0-11 animate (fishx4,birdx4,facex4), 12-23 inanimate
    % (boats, tree, tools), 24-47 words (fish,bird,face,boat,tree,tool)
    % but python indices start at 0
    % word order is lower,lower-italics,upper,upper-italics
    its = 25:2:47; % italics indices
    
    fixedeventlist = neweventlist;
    
    % sub for stimC,stimL and stimR  (minus 1 to get to non-italic same word)
    central=ismember(fixedeventlist.stimC,its);
    fixedeventlist.stimC(central) = fixedeventlist.stimC(central)-1;
    
    left=ismember(fixedeventlist.stimL,its);
    fixedeventlist.stimL(left) = fixedeventlist.stimL(left)-1;
    
    right=ismember(fixedeventlist.stimR,its);
    fixedeventlist.stimR(right) = fixedeventlist.stimR(right)-1;
    
    % save to tsv file
    writetable(neweventlist,'tmp.csv','Delimiter','\t')
    movefile('tmp.csv',eventsfntsv)
    
    %% create individual trial epochs
    EEG_epoch = pop_epoch(EEG_cont, {'E  1'}, [-.1 .8]);
    EEG_epoch = eeg_checkset(EEG_epoch);
    
    %% CSD/laplacian transformation
    fprintf('Performing laplacian transform on sub-%02i data\n',partid);
    data = CSD(single(EEG_epoch.data),G,H);          % compute CSD for <channels-by-samples-by-epochs> 3-D data matrix
    fprintf('Laplacian transform done\n');
    CSD_3D_final = double(data);               % final CSD data

    EEG_CSD = EEG_epoch;
    EEG_CSD.data = CSD_3D_final;
    
    %% convert to cosmo
    fprintf('Converting to cosmo format\n');

    ds = cosmo_flatten(permute(EEG_CSD.data,[3 1 2]),{'chan','time'},{{EEG_CSD.chanlocs.labels},EEG_CSD.times},2);
    ds.a.meeg=struct(); %or cosmo thinks it's not a meeg ds 
    ds.sa = table2struct(fixedeventlist,'ToScalar',true);
    cosmo_check_dataset(ds,'meeg');

    %% save epochs
    fprintf('Saving...\n');

    save(sprintf('%s/derivatives/cosmomvpa/sub-%02i_task-rsvp_cosmomvpa_CSD.mat',datapath,partid),'ds','-v7.3')
    
end
