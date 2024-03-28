function stats_neuralbehaviouralcorr()

%% correlations between neural responses and behaviour
% run RSA from two behavioural tasks to RDMs of left and right hemispheres
% then do contra and ipsi conditions to two behavioural tasks

if isempty(which('cosmo_wtf'))
    addpath('~/Dropbox/MATLAB/CoSMoMVPA/mvpa')
end

plotmodels = 1;

%% get behavioural models

conceptsim = load('./results/stats_concept_similarity.mat');
imagesim = load('./results/stats_image_similarity.mat');

mask = ones(size(conceptsim.RDMmean));
mask = tril(mask,-1);

conc = conceptsim.RDMmean(mask==1);
imsim = imagesim.RDMmean(mask==1);

%% construct models for correlation

mods = [imsim conc]; % take behaviour RDMs and stimulus model RDMs
modelnames = {'ImageTask' 'ConceptTask'};

% plot models
if plotmodels
    figure
    for a = 1:size(mods,2)
        subplot(1,size(mods,2),a)
        imagesc(squareform(mods(:,a)))
        axis square
        axis off
        title(modelnames{a})
    end
end

%% stack results from image decoding
fprintf('Loading data\n')
hfiles = dir('./results/sub-*_half_decoding.mat'); % two halves
names = {hfiles(:).name};
res_cell={};
cc = clock();mm='';
for f=1:length(names)
    fn = sprintf('results/%s',names{f});
    load(fn,'res');
    res_cell{f} = res;
    mm = cosmo_show_progress(cc,f/length(names),sprintf('%i/%i',f,length(names)),mm);
end
res_all = cosmo_stack(res_cell);

timevect = res_all.a.fdim.values{1};

%% correlate neural data per half with models, then take mean
fprintf('Computing neural-behaviour analyses\n')

conds = {'Peripheral' 'MultiplePeripheral'};
clust = {'L_elecs' 'R_elecs'};
ci = {'contra' 'ipsi'};

stats = struct();
stats.format = 'stats.condition.electrodecluster.vf';
stats.timevect = res_all.a.fdim.values{1};
stats.models.names = modelnames;
stats.models.models = mods;

for c = 1:length(conds)

    res_cond = cosmo_slice(res_all,strcmp(res_all.sa.cond,conds{c}));

    %% do neural-behaviour correlations per vf/hemi combination
    vfs = unique(res_cond.sa.vf); % LVF/RVF
    for v = 1:length(vfs)
        for cl = 1:length(clust)
            fprintf('Computing task correlations for %s, %s, %s \n',conds{c},vfs{v},clust{cl})

            s=struct();
            for p = 1:length(names)

                res_p = cosmo_slice(res_cond,res_cond.sa.subject==p);

                %% correlate models with neural data

                % get data per half
                dat1 = cosmo_slice(res_p,strcmp(res_p.sa.vf,vfs{v})&strcmp(res_p.sa.clust,clust{cl})&res_p.sa.half==1);
                dat2 = cosmo_slice(res_p,strcmp(res_p.sa.vf,vfs{v})&strcmp(res_p.sa.clust,clust{cl})&res_p.sa.half==2);
                x1 = dat1.samples;
                x2 = dat2.samples;

                % correlate and neural data with behavioural data
                r1 = corr(x1,mods,'type','Spearman');
                r2 = corr(x2,mods,'type','Spearman');
                rho = (r1 + r2)/2;

                s.rho(:,:,p) = rho;

            end

            %% run stats

            fprintf('Computing stats for %s, %s, %s\n',conds{c},vfs{v},clust{cl})

            dat = s.rho;
            s.condition = conds{c};
            s.eleccluster = clust{cl};
            s.vf = vfs{v};
            s.n = size(dat,3);
            s.mu = mean(dat,3);
            s.mu_all = dat;
            s.se = std(dat,[],3)./sqrt(s.n);
            s.tv = res_cond.a.fdim.values{1};

            % calculate bayesfactors
            s.bf = zeros(length(timevect),length(modelnames));
            s.onset = NaN(2,1);
            s.peak = NaN(2,1);
            for m = 1:length(modelnames)
                s.bf(:,m) = bayesfactor_R_wrapper(squeeze(dat(:,m,:)),...
                    'returnindex',2,'verbose',false,'args','mu=0,rscale="medium",nullInterval=c(-0.5,0.5)');

                % group mean onset & peak
                [mo,ot] = max(movmean(s.bf(:,m)>10,[0 9])==1); % onset 10 consecutive tp BF>10
                [~,peak] = max(s.mu(:,m)); % peak
                ot(mo==0) = []; % remove "onsets" with max=0
                if ~isempty(ot)
                    s.onset(m) = s.tv(ot);
                end
                s.peak(m) = s.tv(peak);

            end

            %% now do difference between model correlations
            mu_all = squeeze(s.mu_all(:,1,:)-s.mu_all(:,2,:));
            s.moddiff.mu = mean(mu_all,2);
            s.moddiff.se = std(mu_all,[],2)./size(mu_all,2);
            s.moddiff.bf = bayesfactor_R_wrapper(s.mu_all,...
                'returnindex',2,'verbose',false,'args','mu=0,rscale="medium",nullInterval=c(-0.5,0.5)');

            %% collate and save
            stats.(conds{c}).(clust{cl}).(vfs{v}) = s;
            fprintf('BFs done for %s, %s, %s\n',conds{c},clust{cl},vfs{v})
            fprintf('Saving\n')
            save('./results/stats_behavcorr_byhalf.mat','stats','timevect', '-v7.3');
        end
    end

    %% now do neural-behaviour correlation for contra/ipsi conditions
    for co = 1:length(ci) % contra/ipsi

        fprintf('Computing task correlations for %s, %s\n',conds{c},ci{co})
        s=struct();
        if co == 1 % contra
            x1=stats.(conds{c}).L_elecs.RVF.mu_all;
            x2=stats.(conds{c}).R_elecs.LVF.mu_all;
        else
            x1=stats.(conds{c}).L_elecs.LVF.mu_all;
            x2=stats.(conds{c}).R_elecs.RVF.mu_all;
        end

        s.mu_all = (x1+x2)/2;

        s.condition = conds{c};
        s.ci = ci{co};
        s.n = size(s.mu_all,3);
        s.mu = mean(s.mu_all,3);
        s.se = std(s.mu_all,[],3)./sqrt(s.n);

        % calculate bayesfactors
        s.bf = zeros(length(timevect),length(modelnames));
        s.onset = NaN(2,1);
        s.peak = NaN(2,1);
        for m = 1:length(modelnames)
            s.bf(:,m) = bayesfactor_R_wrapper(squeeze(s.mu_all(:,m,:)),...
                'returnindex',2,'verbose',false,'args','mu=0,rscale="medium",nullInterval=c(-0.5,0.5)');

            % group mean onset & peak
            [mo,ot] = max(movmean(s.bf(:,m)>10,[0 9])==1); % onset 10 consecutive tp BF>10
            [~,peak] = max(s.mu(:,m)); % peak
            ot(mo==0) = []; % remove "onsets" with max=0
            if ~isempty(ot)
                s.onset(m) = timevect(ot);
            end
            s.peak(m) = timevect(peak);
        end

        %% now do difference between model correlations
        mu_all = squeeze(s.mu_all(:,1,:)-s.mu_all(:,2,:));
        s.moddiff.mu = mean(mu_all,2);
        s.moddiff.se = std(mu_all,[],2)./size(mu_all,2);
        s.moddiff.bf = bayesfactor_R_wrapper(s.mu_all,...
            'returnindex',2,'verbose',false,'args','mu=0,rscale="medium",nullInterval=c(-0.5,0.5)');

        %% now do difference between RH and LH
        mu_all = x2-x1;
        s.rldiff.mu = mean(mu_all,3);
        s.rldiff.n = size(mu_all,3);
        s.rldiff.se = std(mu_all,[],3)./sqrt(s.n);
       % calculate bayesfactors
        s.rldiff.bf = zeros(length(timevect),length(modelnames));
        s.rldiff.onset = NaN(2,1);
        s.rldiff.peak = NaN(2,1);
        for m = 1:length(modelnames)
            s.rldiff.bf(:,m) = bayesfactor_R_wrapper(squeeze(mu_all(:,m,:)),...
                'returnindex',2,'verbose',false,'args','mu=0,rscale="medium",nullInterval=c(-0.5,0.5)');

            % group mean onset & peak
            [mo,ot] = max(movmean(s.rldiff.bf(:,m)>10,[0 9])==1); % onset 10 consecutive tp BF>10
            [~,peak] = max(s.rldiff.mu(:,m)); % peak
            ot(mo==0) = []; % remove "onsets" with max=0
            if ~isempty(ot)
                s.rldiff.onset(m) = timevect(ot);
            end
            s.rldiff.peak(m) = timevect(peak);
        end

        %% collate and save
        stats.(conds{c}).(ci{co}) = s;
        fprintf('BFs done for %s, %s\n',conds{c},ci{co})
        fprintf('Saving\n')
        save('./results/stats_neuralbehavcorr.mat','stats','timevect', '-v7.3');
    end
end

end