function run_ch_searchlight_multiclass(subjectnr)


if ismac
    addpath('~/Dropbox/MATLAB/CoSMoMVPA/mvpa')
    ftpath = '~/Dropbox/MATLAB/fieldtrip-20171104';
    addpath(ftpath)
else
    addpath('../CoSMoMVPA/mvpa');
    addpath('../fieldtrip')
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
ft_defaults;

cosmo_warning('off')

%%
fn = sprintf('data/derivatives/cosmomvpa/sub-%02i_task-rsvp_cosmomvpa_CSD.mat',subjectnr);
outfn = sprintf('results/sub-%02i_ch_searchlight_multiclass.mat',subjectnr);
fprintf('loading %s\n',fn);tic
load(fn,'ds')
fprintf('loading data finished in %i seconds\n',ceil(toc))

%%
ma={};
ma.classifier = @cosmo_classify_lda;
ma.nproc = nproc;

% get channel neighbours
nh1 = cosmo_meeg_chan_neighborhood(ds,'count',4,'label','dataset','label_threshold',.99);
nh2 = cosmo_interval_neighborhood(ds,'time','radius',0);
nh = cosmo_cross_neighborhood(ds,{nh1,nh2});

%% factors
conditionlabels = {'central' 'LVF/RVF' 'LVF+RVF'};
vfs = {'LVF' 'RVF'};

%% for each presentation condition
res_cell = {};
for cond = 1:3 % central, peripheral, multiple peripheral
    
    %% pairwise decoding for 24 images and 12 words
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
        
        dsbv.sa.chunks = dsbv.sa.sequencenumber;
        dsbv.sa.chunks = cosmo_chunkize(dsbv,6);
        
        
        fprintf('subject %i condition: %s \n',subjectnr,conditionlabels{cond})
        
        partitions = cosmo_nchoosek_partitioner(dsbv,1);
        ma.partitions = cosmo_balance_partitions(partitions,dsbv);
        
        r = cosmo_searchlight(dsbv,nh,@cosmo_crossvalidation_measure,ma);
        if cond == 1
            r.sa.vf = {'central'};
            r.sa.vfnum = 0;
        else
            r.sa.vf = vfs(lrvf);
            r.sa.vfnum = lrvf;
        end
        r.sa.cond = cond;
        r.sa.condlabel = conditionlabels(cond);
        r.sa.subject = subjectnr;
        
        res_cell{end+1} = r;
    end
end

%% save
res = cosmo_stack(res_cell);
save(outfn,'res','-v7.3')
