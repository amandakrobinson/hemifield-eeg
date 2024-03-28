function stats_rsa_within()

% function to correlate data halves for a given hemisphere
% eg RDM from odd sequences correlated with RDM from even
% measure of data consistency for contra vs ipsi

% correlate using spearman correlation

if isempty(which('cosmo_wtf'))
    addpath('~/Dropbox (Personal)/MATLAB/CoSMoMVPA/mvpa')
end

%% stack results
fprintf('Loading data\n')
files = dir('results/sub-*_half_decoding.mat');
res_cell={};
cc = clock();mm='';
for f=1:length(files)
    fn = sprintf('results/%s',files(f).name);
    load(fn,'res');
    res_cell{f} = res;
    mm = cosmo_show_progress(cc,f/length(files),sprintf('%i/%i',f,length(files)),mm);
end
res_all = cosmo_stack(res_cell);

timevect = res_all.a.fdim.values{1};

%% computing reliability: half-half rsa correlations within hemispheres/VF conditions
fprintf('Computing stats\n')

conds = {'Peripheral' 'MultiplePeripheral'};
clust = {'L_elecs' 'R_elecs'};

stats = struct();
for c = 1:length(conds)

    res_cond = cosmo_slice(res_all,strcmp(res_all.sa.cond,conds{c}));

    vfs = unique(res_cond.sa.vf); % central or LVF & RVF

    for v = 1:length(vfs)

        for p = 1:length(files)
            fprintf('Computing correlations for participant %d, %s, %s \n',p,conds{c},vfs{v})

            res_p = cosmo_slice(res_cond,res_cond.sa.subject==p);

            % correlate within a hemisphere: consistency
            for hemi = 1:2 % for each hemisphere

                % first half (odd/even)
                dat1 = cosmo_slice(res_p,strcmp(res_p.sa.vf,vfs{v})&strcmp(res_p.sa.clust,clust{hemi})&res_p.sa.half==1);

                % second half (odd/even)
                dat2 = cosmo_slice(res_p,strcmp(res_p.sa.vf,vfs{v})&strcmp(res_p.sa.clust,clust{hemi})&res_p.sa.half==2);

                % get data details and stats for the group
                x1 = dat1.samples;
                x2 = dat2.samples;

                % correlate halves per time point
                bytp = zeros(300,1);
                for t1 = 1:size(x1,2)
                    bytp(t1) = corr(x1(:,t1),x2(:,t1),'type','Spearman');
                end
                
                stats.(conds{c}).(vfs{v}).(clust{hemi}).mu_all(:,p) = bytp;
            end

            % correlate across hemispheres: shared information
            acr = zeros(300,2);
            for h = 1:2 % for each half of data -do twice

                % left hemi, one half of sequences
                dat1 = cosmo_slice(res_p,strcmp(res_p.sa.vf,vfs{v})&strcmp(res_p.sa.clust,clust{1})&res_p.sa.half==h);

                % right hemi, other half of sequences
                dat2 = cosmo_slice(res_p,strcmp(res_p.sa.vf,vfs{v})&strcmp(res_p.sa.clust,clust{2})&res_p.sa.half==(mod(h,2)+1));

                % get data details and stats for the group
                x1 = dat1.samples;
                x2 = dat2.samples;

                % correlate hemispheres per time point
                for t1 = 1:size(x1,2)
                    acr(t1,h) = corr(x1(:,t1),x2(:,t1),'type','Spearman');
                end

                stats.(conds{c}).(vfs{v}).across.raw(:,:,p) = acr;
                stats.(conds{c}).(vfs{v}).across.mu_all(:,p) = mean(acr,2);
            end
        end

        fprintf('Saving\n')
        save('results/stats_within.mat','stats','timevect');

    end
end

%% now combine contra-ipsi results

stats.Peripheral.contra.mu_all = (stats.Peripheral.LVF.R_elecs.mu_all+stats.Peripheral.RVF.L_elecs.mu_all)./2;
stats.Peripheral.ipsi.mu_all = (stats.Peripheral.LVF.L_elecs.mu_all+stats.Peripheral.RVF.R_elecs.mu_all)./2;
stats.Peripheral.across.mu_all = squeeze(mean(stats.Peripheral.LVF.across.raw,2)+mean(stats.Peripheral.RVF.across.raw,2))./2;

stats.MultiplePeripheral.contra.mu_all = (stats.MultiplePeripheral.LVF.R_elecs.mu_all+stats.MultiplePeripheral.RVF.L_elecs.mu_all)./2;
stats.MultiplePeripheral.ipsi.mu_all = (stats.MultiplePeripheral.LVF.L_elecs.mu_all+stats.MultiplePeripheral.RVF.R_elecs.mu_all)./2;
stats.MultiplePeripheral.across.mu_all = squeeze(mean(stats.MultiplePeripheral.LVF.across.raw,2)+mean(stats.MultiplePeripheral.RVF.across.raw,2))./2;

%% do stats for correlations above zero

con1 = [clust 'across'];
con2 = {'contra' 'ipsi' 'across'};

for c = 1:length(conds)
    % first vf x hemi
    for v = 1:length(vfs)
        fprintf('Computing stats for %s, %s\n',conds{c},vfs{v})

        for x = 1:length(con1) % for each hemisphere/across

            dat = stats.(conds{c}).(vfs{v}).(con1{x}).mu_all;
            clear s
            s.mu_all = dat;
            s.n = size(dat,2);
            s.mu = mean(dat,2);
            s.se = std(dat,[],2)./sqrt(s.n);

            % calculate bayesfactors
            s.bf = bayesfactor_R_wrapper(dat,...
                'returnindex',2,'verbose',false,...
                'args','mu=0,rscale="medium",nullInterval=c(-0.5,0.5)');
            stats.(conds{c}).(vfs{v}).(con1{x}) = s;
        end
    end

    %% now do contra/ipsi/across
    for ci = 1:length(con2)
        dat = stats.(conds{c}).(con2{ci}).mu_all;
        clear s
        s.mu_all = dat;
        s.n = size(dat,2);
        s.mu = mean(dat,2);
        s.se = std(dat,[],2)./sqrt(s.n);

        % calculate bayesfactors
        s.bf = bayesfactor_R_wrapper(dat,...
            'returnindex',2,'verbose',false,...
            'args','mu=0,rscale="medium",nullInterval=c(-0.5,0.5)');
        stats.(conds{c}).(con2{ci}) = s;
    end

    %% now do differences between contra-across and ipsi-across
    for ci = 1:(length(con2)-1)
        dat1 = stats.(conds{c}).(con2{ci}).mu_all; % contra/ipsi
        dat2 = stats.(conds{c}).across.mu_all; % across
        clear s
        dat = dat1-dat2;
        s.mu_all = dat;
        s.n = size(dat,2);
        s.mu = mean(dat,2);
        s.se = std(dat,[],2)./sqrt(s.n);

        % calculate bayesfactors
        s.bf = bayesfactor_R_wrapper(dat,...
            'returnindex',2,'verbose',false,...
            'args','mu=0,rscale="medium",nullInterval=c(-0.5,0.5)');
        stats.(conds{c}).withinVacross.(con2{ci}) = s;
    end
end

fprintf('Saving\n')
save('results/stats_within.mat','stats','timevect');
fprintf('Done\n')

