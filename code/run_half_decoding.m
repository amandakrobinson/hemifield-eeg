function run_half_decoding(subjectnr)

%% decode stimuli from EEG data for hemifield experiment
% separately for clusters over left and right hemispheres
% decoding is performed separately on odd and even trial sequences,
% for later comparison across hemispheres without using the same trials

%% cosmomvpa
if ismac
    addpath('~/Dropbox (Personal)/MATLAB/CoSMoMVPA/mvpa')

else
    addpath('../CoSMoMVPA/mvpa');
    fprintf('paths loaded, getting nprocs\n');

    % start cluster, give it a unique directory
    % starting a pool can fail when 2 procs are requesting simultaneous
    % thus try again after a second until success
        pool=[];
        while isempty(pool)
            try
                pc = parcluster('local');
                pc.JobStorageLocation=tempdir;
                pool=parpool(pc);
            catch err
                disp(err)
                delete(gcp('nocreate'));
                pause(1)
            end
        end

end
nproc=cosmo_parallel_get_nproc_available();
fprintf('%d procs available\n',nproc);
cosmo_warning('off')

%% data details
fprintf('loading data...\n')

fn = sprintf('data/derivatives/cosmomvpa/sub-%02i_task-rsvp_cosmomvpa_CSD.mat',subjectnr);
outfn = sprintf('results/sub-%02i_half_decoding.mat',subjectnr);
fprintf('loading %s\n',fn);tic
load(fn,'ds')
fprintf('loading data finished in %i seconds\n',ceil(toc))
timevect=ds.a.fdim.values{2};

%% get stim combinations
ds.sa.targets = ds.sa.stimC; %temp
conds = unique(ds.sa.targets(ds.sa.targets~=-1));
condcomb = combnk(conds,2);

% get stim names/combos
imnums = (ds.sa.stimL(ds.sa.stimL~=-1));
nims = unique(imnums);
ims = ds.sa.stimnameL(ds.sa.stimL~=-1);
imnames={};
for i = 1:length(nims)
    n = unique(ims(imnums==nims(i)));
    imnames(i) = n(1);
end

%% factors
condnames = {'Central' 'Peripheral' 'MultiplePeripheral'};
vfs = {'LVF' 'RVF'};

%% elecs
elecsel = [14:16 44:46; 18:20 48:50];
elecs = {'L_elecs' 'R_elecs'};

%% classifier details
ma={};
ma.classifier = @cosmo_classify_lda;
ma.nproc = nproc;
ma.progress=0;
ma.normalization='demean';
ma.output = 'accuracy';

%% decode
res_cell={};

%% for each presentation condition
for cond = 1:length(condnames) % central, peripheral, multiple peripheral
    
    %% pairwise decoding for 24 images and 24 words
    % order of images is 1-12 animate (fishx4,birdx4,facex4), 13-24 inanimate
    % (boats, tree, tools), 25-36 words (fish, bird, face, boats, tree, tools)
    
    dsb = cosmo_slice(ds,ds.sa.condition==cond); % central objs

    %% for each visual field
    for lrvf=1:length(vfs) % LVF then RVF
        
        if cond==1 % only run central through once with stimC (not L/R VF)
            dsb.sa.targets = dsb.sa.stimC;
            if lrvf == 2
                break
            end
        else % peripheral conditions
            if lrvf == 1
                dsb.sa.targets = dsb.sa.stimL;
            else
                dsb.sa.targets = dsb.sa.stimR;
            end
        end
        
        dsbv = cosmo_slice(dsb,dsb.sa.targets~=-1);
        
        %% decode by hemisphere
        for clust = 1:2 % left cluster then right cluster
            
            % subset for specific electrodes
            dsh = cosmo_slice(dsbv,ismember(dsbv.fa.chan,elecsel(clust,:)),2);
            
            % get stream numbers
            streams = unique(dsh.sa.streamnumber);
            
            for h = 1:2 % decode half then other half of sequences
                
                dsbe = cosmo_slice(dsh,ismember(dsh.sa.streamnumber,streams(h:2:end)));
                
                % make chunks
                dsbe.sa.chunks = dsbe.sa.sequencenumber;
                
                nh = cosmo_interval_neighborhood(dsbe,'time','radius',0);
                
                %% parallel decode
                fprintf('p%d decoding %s, half %d using %s\n',subjectnr,condnames{cond},h,elecs{clust})

                for t=1:length(condcomb)
                    
                    % decode cv correct trials
                    ds_prop = cosmo_slice(dsbe,ismember(dsbe.sa.targets,condcomb(t,:)));
                    
                    partitions = cosmo_nchoosek_partitioner(ds_prop,1);
                    partitions = cosmo_balance_partitions(partitions,ds_prop,'balance_test',false);
                    
                    ma.partitions=partitions;
                    
                    r = cosmo_searchlight(ds_prop,nh,@cosmo_crossvalidation_measure,ma);
                    r.sa.subject = subjectnr;
                    if cond == 1
                        r.sa.vf = {'central'};
                    else
                        r.sa.vf = vfs(lrvf);
                    end
                    r.sa.clust = elecs(clust);
                    r.sa.cond = condnames(cond);
                    r.sa.half = h;
                    r.sa.cond1 = condcomb(t,1);
                    r.sa.cond2 = condcomb(t,2);
                    r.sa.im1 = imnames(nims==condcomb(t,1));
                    r.sa.im2 = imnames(nims==condcomb(t,2));
                    
                    res_cell{end+1} = r;
       
                end
            end
        end 
    end
end

res = cosmo_stack(res_cell);
save(outfn,'res','timevect','-v7.3')

end