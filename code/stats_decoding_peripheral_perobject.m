function stats_decoding_peripheral_perobject()

% ----- run stats split by face/object/word: peripheral conditions -------
% all condtions versus chance (one-tailed v 0.5)
% contra-ipsi or LVF/RVF for each hemisphere (two-tailed v 0)
% also RH-LH for contra/ipsi (two-tailed v 0)

if isempty(which('cosmo_wtf'))
    addpath('~/CoSMoMVPA/mvpa')
end

%% stack results and collate
fprintf('Loading data\n')

% get files
fileList = dir('results/sub-*_decoding.mat');
nofiles = dir('results/sub-*_half_decoding.mat'); % don't include split half decoding files
filenames = {fileList.name};
idx = ismember(filenames,{nofiles.name});
files = filenames(~idx);

doonsets=0;
combs = nchoosek(1:length(files),(length(files)-2)); %leave-2-participants-out for jackknifing

% conds
objconds = {'objects','faces','words'};
% condcomb = nchoosek([0:23 24:2:46],2);
stimnum.objects = [0:7 12:23];
stimnum.faces = 8:11;
stimnum.words = 24:2:46;

% collate
res_cell={};
cc = clock();mm='';
for f=1:length(files)
    fn = sprintf('results/%s',files{f});
    load(fn,'res');
    res.sa.objtype = cell(size(res.sa.clust));
    res.sa.objnum = zeros(size(res.sa.clust));
    for o = 1:length(objconds)
        idx = find(ismember(res.sa.cond1,stimnum.(objconds{o}))&...
            ismember(res.sa.cond2,stimnum.(objconds{o})));
        res.sa.objtype(idx) = objconds(o);
        res.sa.objnum(idx) = o;
    end
    res = cosmo_slice(res,res.sa.objnum>0);
    res_cell{f} = res;
    mm = cosmo_show_progress(cc,f/length(files),sprintf('%i/%i',f,length(files)),mm);
end
res_all = cosmo_stack(res_cell);

% average across all pairs
res_all = cosmo_average_samples(res_all,'split_by',{'clust','cond','objnum','subject','vf'});

%% set up stats structure
conds = {'Peripheral' 'MultiplePeripheral'};
clust = {'L_elecs' 'R_elecs'};
conip = {'contra' 'ipsi'};

stats = struct();
stats.format = 'stats.condition.objtype.electrodecluster.vf';
stats.timevect = res_all.a.fdim.values{1};

%% run stats
fprintf('Computing stats\n')

for c = 1:length(conds)

    res_cond = cosmo_slice(res_all,strcmp(res_all.sa.cond,conds{c}));

    vfs = unique(res_cond.sa.vf);
    for o = 1:length(objconds)
        for cl = 1:length(clust)
            for v = 1:length(vfs)

                %% get means and do stats
                fprintf('Computing stats for %s, %s, %s cluster, %s \n',conds{c},objconds{o},clust{cl},vfs{v})

                dat = cosmo_slice(res_cond,strcmp(res_cond.sa.vf,vfs{v})&strcmp(res_cond.sa.clust,clust{cl})&res_cond.sa.objnum==o);

                % get data details and stats for the group
                x = dat.samples;
                s = struct();
                s.condition = conds{c};
                s.eleccluster = clust{cl};
                s.objtype = unique(dat.sa.objtype);
                s.vf = vfs{v};
                s.n = size(x,1);
                s.mu = mean(x,1);
                s.mu_all = x;
                s.se = std(x)./sqrt(s.n);

                % calculate bayesfactors
                s.bf = bayesfactor_R_wrapper(x',...
                    'returnindex',2,'verbose',false,'args','mu=0.5,rscale="medium",nullInterval=c(-Inf,0.5)');

                s.tv = res_all.a.fdim.values{1};

                % group mean onset & peak
                [mo,ot] = max(movmean(s.bf>10,[0 9])==1); % onset 10 consecutive tp BF>10
                [~,peak] = max(s.mu); % peak
                ot(mo==0) = []; % remove "onsets" with max=0
                s.onset = stats.timevect(ot);
                s.peak = stats.timevect(peak);

                stats.(conds{c}).(objconds{o}).(clust{cl}).(vfs{v}) = s;
                fprintf('BFs done for %s, %s, %s cluster, %s\n',conds{c},objconds{o},clust{cl},vfs{v})

                %% calculate onsets and peaks by jackknifing
                if doonsets == 1
                    fprintf('Calculating onsets for %s, %s cluster, %s \n',conds{c},clust{cl},vfs{v})

                    % use jackknife method and calculate bayesfactors to get onsets
                    x = [];
                    for i=1:size(combs,1)
                        x = [x s.mu_all(combs(i,:),:)];
                    end
                    tic
                    b = bayesfactor_R_wrapper(x',...
                        'returnindex',2,'verbose',false,'args','mu=0.5,rscale="medium",nullInterval=c(-Inf,0.5)');
                    s.bf_jackknife = reshape(b,[],size(combs,1));
                    fprintf('jackknife BFs done for %s, %s cluster, %s \n',conds{c},clust{cl},vfs{v})
                    toc

                    % get onsets and peaks per jackknife
                    onset_jack = [];peak_jack=[];maxon=[];
                    for b=1:size(s.bf_jackknife,2)
                        [maxon(b),onset_jack(b)] = max(movmean(s.bf_jackknife(:,b)>10,[0 9])==1); % ten points in a row
                        [~,peak_jack(b)] = max(mean(s.mu_all(combs(b,:),:),1));
                    end
                    onset_jack(maxon==0) = []; % remove "onsets" with max=0

                    % convert onsets and peaks to time (rather than tp)
                    s.onsetjk = stats.timevect(onset_jack);
                    if isempty(onset_jack)
                        s.onsetci=[];
                    else
                        ci = prctile(onset_jack,[5 95]);
                        s.onsetci = stats.timevect([floor(ci(1)) ceil(ci(2))]);
                    end
                    s.peakjk = stats.timevect(peak_jack);
                    s.peakci = stats.timevect(round(prctile(peak_jack,[5 95])));

                    stats.(conds{c}).(objconds{o}).(clust{cl}).(vfs{v}) = s;
                    fprintf('Onsets done\n')
                end
            end

            %% calculate vf differences
            dat_cond = stats.(conds{c}).(objconds{o}).(clust{cl});

            % get data details and stats for the group difference
            datdiff = dat_cond.(vfs{1}).mu_all - dat_cond.(vfs{2}).mu_all; % LVF-RVF
            s = struct();
            s.condition = conds{c};
            s.eleccluster = clust{cl};
            s.comparison = 'LVF-RVF difference';
            s.n = size(datdiff,1);
            s.mu = mean(datdiff,1);
            s.mu_all = datdiff;
            s.se = std(datdiff)./sqrt(s.n);

            % calculate bayesfactors
            s.bf = bayesfactor_R_wrapper(datdiff',...
                'returnindex',2,'verbose',false,'args','mu=0,rscale="medium",nullInterval=c(-0.5,0.5)');

            % group mean onset & peak
            [mo,ot] = max(movmean(s.bf>10,[0 9])==1); % onset 10 consecutive tp BF>10
            [~,peak] = max(s.mu); % peak
            ot(mo==0) = []; % remove "onsets" with max=0
            s.onset = stats.timevect(ot);
            s.peak = stats.timevect(peak);

            if doonsets == 1
                fprintf('Calculating onsets for VF diff: %s, %s cluster\n',conds{c},clust{cl})

                % use jackknife method and calculate bayesfactors to get onsets
                x = [];
                for i=1:size(combs,1)
                    x = [x s.mu_all(combs(i,:),:)];
                end
                tic
                b = bayesfactor_R_wrapper(x',...
                    'returnindex',2,'verbose',false,'args','mu=0,rscale="medium",nullInterval=c(-0.5,0.5)');
                s.bf_jackknife = reshape(b,[],size(combs,1));
                fprintf('jackknife BFs done for %s, %s cluster, %s \n',conds{c},clust{cl},vfs{v})
                toc

                % get onsets and peaks per jackknife
                onset_jack = [];peak_jack=[];maxon=[];
                for b=1:size(s.bf_jackknife,2)
                    [maxon(b),onset_jack(b)] = max(movmean(s.bf_jackknife(:,b)>10,[0 9])==1); % ten points in a row
                    [~,peak_jack(b)] = max(mean(s.mu_all(combs(b,:),:),1));
                end
                onset_jack(maxon==0) = []; % remove "onsets" with max=0

                % convert onsets and peaks to time (rather than tp)
                s.onsetjk = stats.timevect(onset_jack);
                if isempty(onset_jack)
                    s.onsetci=[];
                else
                    ci = prctile(onset_jack,[5 95]);
                    s.onsetci = stats.timevect([floor(ci(1)) ceil(ci(2))]);
                end
                s.peakjk = stats.timevect(peak_jack);
                s.peakci = stats.timevect(round(prctile(peak_jack,[5 95])));

                fprintf('Onsets done\n')
            end
            stats.(conds{c}).(objconds{o}).(clust{cl}).LVFRVFdiff = s;
        end
        save('results/stats_decoding_peripheral_perobs.mat','stats');

        %% RH-LH hemispheric differences

        for ci = 1:length(conip) % contra then ipsi
            if ci == 1 % contra
                datdiff = stats.(conds{c}).(objconds{o}).R_elecs.LVF.mu_all - stats.(conds{c}).(objconds{o}).L_elecs.RVF.mu_all; % RH-LH contra
            else % ipsi
                datdiff = stats.(conds{c}).(objconds{o}).R_elecs.RVF.mu_all - stats.(conds{c}).(objconds{o}).L_elecs.LVF.mu_all; % RH-LH ipsi
            end

            s = struct();
            s.condition = conds{c};
            s.comparison = 'RH-LH difference';
            s.contraipsi = conip{ci};
            s.n = size(datdiff,1);
            s.mu = mean(datdiff,1);
            s.mu_all = datdiff;
            s.se = std(datdiff)./sqrt(s.n);

            % calculate bayesfactors
            s.bf = bayesfactor_R_wrapper(datdiff',...
                'returnindex',2,'verbose',false,'args','mu=0,rscale="medium",nullInterval=c(-0.5,0.5)');

            % group mean onset & peak
            [mo,ot] = max(movmean(s.bf>10,[0 9])==1); % onset 10 consecutive tp BF>10
            [~,peak] = max(s.mu); % peak
            ot(mo==0) = []; % remove "onsets" with max=0
            s.onset = stats.timevect(ot);
            s.peak = stats.timevect(peak);

            if doonsets == 1
                fprintf('Calculating onsets for hemi diff: %s\n',conds{c})

                % use jackknife method and calculate bayesfactors to get onsets
                x = [];
                for i=1:size(combs,1)
                    x = [x s.mu_all(combs(i,:),:)];
                end
                tic
                b = bayesfactor_R_wrapper(x',...
                    'returnindex',2,'verbose',false,'args','mu=0,rscale="medium",nullInterval=c(-0.5,0.5)');
                s.bf_jackknife = reshape(b,[],size(combs,1));
                fprintf('jackknife BFs done for %s, %s cluster, %s \n',conds{c},clust{cl},vfs{v})
                toc

                % get onsets and peaks per jackknife
                onset_jack = [];peak_jack=[];maxon=[];
                for b=1:size(s.bf_jackknife,2)
                    [maxon(b),onset_jack(b)] = max(movmean(s.bf_jackknife(:,b)>10,[0 9])==1); % ten points in a row
                    [~,peak_jack(b)] = max(mean(s.mu_all(combs(b,:),:),1));
                end
                onset_jack(maxon==0) = []; % remove "onsets" with max=0

                % convert onsets and peaks to time (rather than tp)
                s.onsetjk = stats.timevect(onset_jack);
                if isempty(onset_jack)
                    s.onsetci=[];
                else
                    conf = prctile(onset_jack,[5 95]);
                    s.onsetci = stats.timevect([floor(conf(1)) ceil(conf(2))]);
                end
                s.peakjk = stats.timevect(peak_jack);
                s.peakci = stats.timevect(round(prctile(peak_jack,[5 95])));

                fprintf('Onsets done\n')
            end

            stats.(conds{c}).(objconds{o}).RHLHdiff.(conip{ci}) = s;
            save('results/stats_decoding_peripheral_perobs.mat','stats');

        end

        %% finally do contra-ipsi for RH-LH
        datdiff = stats.(conds{c}).(objconds{o}).RHLHdiff.contra.mu_all - stats.(conds{c}).(objconds{o}).RHLHdiff.ipsi.mu_all; % contra-ipsi for RH-LH

        s = struct();
        s.condition = conds{c};
        s.comparison = 'RH-LH contra-ipsi difference';
        s.n = size(datdiff,1);
        s.mu = mean(datdiff,1);
        s.mu_all = datdiff;
        s.se = std(datdiff)./sqrt(s.n);

        % calculate bayesfactors
        s.bf = bayesfactor_R_wrapper(datdiff',...
            'returnindex',2,'verbose',false,'args','mu=0,rscale="medium",nullInterval=c(-0.5,0.5)');

        % group mean onset & peak
        [mo,ot] = max(movmean(s.bf>10,[0 9])==1); % onset 10 consecutive tp BF>10
        [~,peak] = max(s.mu); % peak
        ot(mo==0) = []; % remove "onsets" with max=0
        s.onset = stats.timevect(ot);
        s.peak = stats.timevect(peak);

        if doonsets == 1
            fprintf('Calculating onsets for contra-ipsi hemi diff: %s\n',conds{c})

            % use jackknife method and calculate bayesfactors to get onsets
            x = [];
            for i=1:size(combs,1)
                x = [x s.mu_all(combs(i,:),:)];
            end
            tic
            b = bayesfactor_R_wrapper(x',...
                'returnindex',2,'verbose',false,'args','mu=0,rscale="medium",nullInterval=c(-0.5,0.5)');
            s.bf_jackknife = reshape(b,[],size(combs,1));
            fprintf('jackknife BFs done for %s, %s cluster, %s \n',conds{c},clust{cl},vfs{v})
            toc

            % get onsets and peaks per jackknife
            onset_jack = [];peak_jack=[];maxon=[];
            for b=1:size(s.bf_jackknife,2)
                [maxon(b),onset_jack(b)] = max(movmean(s.bf_jackknife(:,b)>10,[0 9])==1); % ten points in a row
                [~,peak_jack(b)] = max(mean(s.mu_all(combs(b,:),:),1));
            end
            onset_jack(maxon==0) = []; % remove "onsets" with max=0

            % convert onsets and peaks to time (rather than tp)
            s.onsetjk = stats.timevect(onset_jack);
            if isempty(onset_jack)
                s.onsetci=[];
            else
                ci = prctile(onset_jack,[5 95]);
                s.onsetci = stats.timevect([floor(ci(1)) ceil(ci(2))]);
            end
            s.peakjk = stats.timevect(peak_jack);
            s.peakci = stats.timevect(round(prctile(peak_jack,[5 95])));

            fprintf('Onsets done\n')
        end

        stats.(conds{c}).(objconds{o}).RHLHdiff.conipdiff = s;
        save('results/stats_decoding_peripheral_perobs.mat','stats');

    end
end

fprintf('Saving\n')
save('results/stats_decoding_peripheral_perobs.mat','stats');
fprintf('Done\n')


