function stats_rsa_across()

% assess similarity of representations in the left and right hemispheres
% use decoding performed on two halves of the data to make 2x RDMs
% correlate RDM from hemisphere 1 (half 1) to hemisphere 2 (half 2) and vice versa
% (do on separate halves to avoid inflation of correlation due to
% correlated decoding performance from noise)
% use spearman correlation
% take mean of two correlations LHh1:RHh2 and LHh2:RHh1

if isempty(which('cosmo_wtf'))
    addpath('~/CoSMoMVPA/mvpa')
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

%% computing left-right rsa hemisphere correlations
fprintf('Computing stats\n')

conds = {'Central' 'Peripheral' 'MultiplePeripheral'};
clust = {'L_elecs' 'R_elecs'};

stats = struct();
stats.dims = 'LHtime_RHtime_participant';
stats.measure = 'Spearman correlations';
stats.conds = conds;
stats.electrodeclusters = clust;

for c = 1:length(conds)

    res_cond = cosmo_slice(res_all,strcmp(res_all.sa.cond,conds{c}));

    vfs = unique(res_cond.sa.vf); % central or LVF & RVF

    for v = 1:length(vfs)

        s = zeros(length(timevect),length(timevect),length(files));
        for p = 1:length(files)
            fprintf('Computing correlations for participant %d, %s, %s \n',p,conds{c},vfs{v})

            res_p = cosmo_slice(res_cond,res_cond.sa.subject==p);

            %% correlate LH half 1 to RH half 2 and LH half 2 to RH half 1
            r = nan(length(timevect),length(timevect),length(clust));
            for h = 1:2 % for each split-half cross-validation

                % left hemisphere
                dat1 = cosmo_slice(res_p,strcmp(res_p.sa.vf,vfs{v})&strcmp(res_p.sa.clust,clust{1})&res_p.sa.half==h);

                % right hemisphere
                dat2 = cosmo_slice(res_p,strcmp(res_p.sa.vf,vfs{v})&strcmp(res_p.sa.clust,clust{2})&res_p.sa.half==(mod(h,2)+1));

                % get data details and stats for the group
                x1 = dat1.samples;
                x2 = dat2.samples;

                % get correlations
                r(:,:,h) = corr(x1,x2,'type','spearman'); % rows x cols is LH x RH

            end

            % take mean of h1 and h2
            s(:,:,p) = mean(r,3);

            subplot(4,5,p)
            imagesc(mean(r,3));colorbar;

        end
        stats.data.(conds{c}).(vfs{v}) = s;
        fprintf('Saving\n')

    end
end
save('results/stats_across_hemi_correlations.mat','stats');

%% now combine VF results so contra-ipsi results and do stats

rsa.Central = stats.data.Central.central; %lh x rh

rsa.periphLVF = permute(stats.data.Peripheral.LVF,[2 1 3]); %make lvh rh x lh = contra x ipsi
rsa.periphRVF = stats.data.Peripheral.RVF; % rvf lh x rh = contra x ipsi
rsa.Peripheral =  mean(cat(4,rsa.periphLVF,rsa.periphRVF),4); % mean of two VFs

multLVF = permute(stats.data.MultiplePeripheral.LVF,[2 1 3]);
multRVF = stats.data.MultiplePeripheral.RVF;
rsa.MultiplePeripheral = mean(cat(4,multLVF,multRVF),4);

combinedstats=struct();
combinedstats.dims = 'LHtime_RHtime_participant';
combinedstats.measure = 'Spearman correlations';
combinedstats.conds = conds;
combinedstats.electrodeclusters = clust;
fprintf('Pre-saving\n')
save('results/stats_across_hemi_correlations.mat','combinedstats','stats','timevect');

for c = 1:length(conds)

    dat = rsa.(conds{c});
    clear s
    s.n = size(dat,3);
    s.tv = res_all.a.fdim.values{1};

    for t = 1:size(dat,1) % for each LH/contra hemi time point

        fprintf('Computing stats for %s, timepoint %d of %d\n',conds{c},t,size(dat,1))

        x = squeeze(dat(t,:,:))';

        s.mu(t,:) = mean(x);
        s.dat(t,:,:) = x';
        s.se(t,:) = std(x)./sqrt(s.n);

        % calculate bayesfactors
        s.bf(t,:) = bayesfactor_R_wrapper(x',...
            'returnindex',2,'verbose',false,...
            'args','mu=0,rscale="medium",nullInterval=c(-0.5,0.5)');

    end
    combinedstats.data.(conds{c}) = s;
end

timevect = s.tv;
fprintf('Saving\n')
save('results/stats_across_hemi_correlations.mat','combinedstats','stats','timevect');
fprintf('Done\n')

