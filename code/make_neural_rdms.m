function make_neural_rdms()

if isempty(which('cosmo_wtf'))
    addpath('~/Dropbox/MATLAB/CoSMoMVPA/mvpa')
end
load('imageorder.mat')

%% stack results
fprintf('Loading data\n')
ps = 1:20;
res_cell={};
cc = clock();mm='';
for f=1:length(ps)
    fn = sprintf('results/sub-%02i_decoding.mat',ps(f));
    load(fn,'res');
    res_cell{f} = res;
    mm = cosmo_show_progress(cc,f/length(ps),sprintf('%i/%i',f,length(ps)),mm);
end
res_all = cosmo_stack(res_cell);

%% collate left-right RDMs

conds = {'Central' 'Peripheral' 'MultiplePeripheral'};
clust = {'L_elecs' 'R_elecs'};

alldat = struct();
alldat.imagecombs = combnk(0:1:35,2);
alldat.time = res_all.a.fdim.values{1};

for c = 1:length(conds)

    res_cond = cosmo_slice(res_all,strcmp(res_all.sa.cond,conds{c}));

    vfs = unique(res_cond.sa.vf); % central or LVF & RVF

    for v = 1:length(vfs)

        for p = 1:length(ps)
            fprintf('Computing RDMs for participant %d, %s, %s \n',p,conds{c},vfs{v})

            res_p = cosmo_slice(res_cond,res_cond.sa.subject==p);

            % left hemisphere
            dat1 = cosmo_slice(res_p,strcmp(res_p.sa.vf,vfs{v})&strcmp(res_p.sa.clust,clust{1}));

            % right hemisphere
            dat2 = cosmo_slice(res_p,strcmp(res_p.sa.vf,vfs{v})&strcmp(res_p.sa.clust,clust{2}));

            % get data details and stats for the group
            alldat.(conds{c}).(vfs{v}).(clust{1})(:,:,p) = dat1.samples;
            alldat.(conds{c}).(vfs{v}).(clust{2})(:,:,p) = dat2.samples;

        end
        fprintf('Saving\n')
        save('results/rdms.mat','alldat');

    end
end
alldat.ims = [dat1.sa.cond1 dat2.sa.cond2];
alldat.imnames = [dat1.sa.im1 dat2.sa.im2];
alldat.measure = 'pairwise_decoding_accuracy';

fprintf('Saving\n')
save('results/rdms.mat','alldat');

%% periph by VF and hemi

addpath('~/Dropbox (Personal)/MATLAB/Colormaps/colormaps (5)/Colormaps/')
colormap(viridis(20))
rdm = struct();

perdat = alldat.Peripheral;
rdm.periph_LVFcontra = mean(perdat.LVF.R_elecs,3);
rdm.periph_LVFipsi = mean(perdat.LVF.L_elecs,3);
rdm.periph_RVFcontra = mean(perdat.RVF.L_elecs,3);
rdm.periph_RVFipsi = mean(perdat.RVF.R_elecs,3);

rdm.central_LH = mean(alldat.Central.central.L_elecs,3);
rdm.central_RH = mean(alldat.Central.central.R_elecs,3);

perdat = alldat.Peripheral;
rdm.periph_contra = (mean(perdat.LVF.R_elecs,3)+mean(perdat.RVF.L_elecs,3))/2;
rdm.periph_ipsi = (mean(perdat.LVF.L_elecs,3)+mean(perdat.RVF.R_elecs,3))/2;

mpdat = alldat.MultiplePeripheral;
rdm.multperiph_contra = (mean(mpdat.LVF.R_elecs,3)+mean(mpdat.RVF.L_elecs,3))/2;
rdm.multperiph_ipsi = (mean(mpdat.LVF.L_elecs,3)+mean(mpdat.RVF.R_elecs,3))/2;

fprintf('Saving\n')
save('results/rdms.mat','alldat','rdm');
fprintf('Done\n')


