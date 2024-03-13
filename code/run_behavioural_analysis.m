
%% analyse data from online similarity experiments using triplet task

% read and re-map stimuli
% reorders stimuli so images then words, in same object order

fprintf('re-map\n')
stimfns = dir('stimuli/*.png');
stimnames = {stimfns.name};
stims={};
stimorder = {'fish','bird','face','boats','tree','tools'};
wordorder = {'lowercase','uppercase'};
nimstims = 24;
for s=1:length(stimfns)
    name = stimfns(s).name;
    namespl = strsplit(name,{'_','.'});
    
    if ~contains(name,'stim')
        stimtype = stimnames(contains(stimnames,'stim0')&contains(stimnames,name(1:4)));
        stimparts = strsplit(stimtype{1},{'/','_'});
        
        wordnum = find(contains(wordorder,namespl{2}));
        objnum = find(contains(stimorder,namespl{1}));
        imnum(s) = nimstims+(objnum-1)*2+wordnum;
        %         name = sprintf('stim%03i_%s_%s_%s.png',imnum,stimparts{3},stimorder{objnum},wordorder{wordnum});
    else
        imnum(s) = str2num(namespl{1}(5:7));
    end
    stims{imnum(s),1} = fullfile('stimuli/',name);
end


%% now read similarity data from different tasks and collate into RDM

f=figure(1);clf;
f.Position = [f.Position(1:2) 1000 400];

f=figure(2);clf;
f.Position = [f.Position(1:2) 1000 800];

tasks = {'image_similarity' 'concept_similarity'};

for t = 1:length(tasks)
    fns = dir(sprintf('%s/data/hemiraptor-triplets*.csv',tasks{t}));
    
    T = [];subnr = 0;
    Ta=[];
    for f=1:numel(fns)
        fn = fullfile(fns(f).folder,fns(f).name);
        TS = readtable(fn);
        if size(TS,1)>1
            subnr = subnr+1;
            TS.subjectnr(:) = subnr;
            TS.stimulus = [];
            TS = TS(strcmp(TS.test_part,'triplet'),:);
            %TS.button_pressed = str2double(TS.button_pressed);
            T = [T; TS];
        end
        fprintf('%s, file %i/%i\n',tasks{t},f,numel(fns));
    end
    
    %% reorder similarity data
    
    [~,T.stim0number]=ismember(T.stim0,stims);
    [~,T.stim1number]=ismember(T.stim1,stims);
    [~,T.stim2number]=ismember(T.stim2,stims);
    
    %% create RDM
    fprintf('create RDM\n')
    allcombs = [T.stim0number T.stim1number T.stim2number];
    choice = T.button_pressed;
    RDMsum = zeros(length(stims));
    RDMcounts = zeros(length(stims));
    for i=1:numel(choice)
        % for every choice, add 1 to the similarity of the two items that were
        % not the odd-one-out (three choices were coded as 0, 1, 2)
        v = allcombs(i,(0:2)~=choice(i));
        RDMsum(v(1),v(2)) = RDMsum(v(1),v(2))+1;
        RDMsum(v(2),v(1)) = RDMsum(v(2),v(1))+1;
        % add 1 to the counts of all items compared, to compute the mean later
        for v = combnk(allcombs(i,:),2)'
            RDMcounts(v(1),v(2)) = RDMcounts(v(1),v(2))+1;
            RDMcounts(v(2),v(1)) = RDMcounts(v(2),v(1))+1;
        end
    end
    RDMmean = 1 - RDMsum./RDMcounts;
    
    save(sprintf('results/stats_%s.mat',tasks{t}),'RDMmean')
    
end
