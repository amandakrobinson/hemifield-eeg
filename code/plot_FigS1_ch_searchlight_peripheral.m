
%% plot Figure S1 topographies

addpath('~/Dropbox/MATLAB/CoSMoMVPA/mvpa')
ftpath = '~/Dropbox/MATLAB/fieldtrip-20171104';
if ~isempty(which('ft_defaults'))
    rmpath(genpath(ftpath))
end
addpath(ftpath)
ft_defaults;
colormap viridis

%%
res_cell={};cc=clock();mm='';
for s=1:20
    fn = sprintf('results/sub-%02i_ch_searchlight_multiclass.mat',s);
    try
        load(fn,'res')
        res_cell{end+1} = res;
    catch
    end
    mm=cosmo_show_progress(cc,s/20,[],mm);
end
res = cosmo_stack(res_cell);

%% ch searchlight topomaps

timewins = (0:90:270)'+5*[-1 1];

layout=cosmo_meeg_find_layout(res);

h0mean = 1/36;
stc={'Peripheral: left visual field','Peripheral: right visual field',...
    'Multiple peripheral: left visual field','Multiple peripheral: right visual field'};

conds = {'LVF/RVF','LVF+RVF'};

plotnr=0;
f=figure(1);clf

f.Position=[f.Position(1:2) 1000 800];
% f.Position=[f.Position(1:2) 1000 1200];
f.Resize='off';
f.PaperPositionMode='auto';f.PaperOrientation='portrait';
n=0;
for ttt=1:length(timewins)
    for c=1:length(conds)
        xx = cosmo_slice(res,contains(res.sa.condlabel,conds{c}));
        vfs = unique(xx.sa.vfnum);

        for v = 1:length(vfs)
            n=n+1;


            plotnr=plotnr+1;
            subplot(length(timewins),4,plotnr)

            x = cosmo_slice(xx,xx.sa.vfnum==vfs(v));
            x.samples = x.samples-h0mean;
%             mrange = cosmo_average_samples(res,'split_by',{'cond'});
%             mrange = mrange.samples-h0mean;
%             mrange = [0 prctile(mrange(:)',99)];

            % map to FT struct for visualization
            ft = ft_timelockanalysis([],cosmo_map2meeg(x));

            co = viridis();
            % show figure with plots for each sensor
            cfg = [];
            cfg.interactive = 'yes';
            cfg.zlim = [0 0.01];%mrange;
            cfg.xlim = timewins(ttt,:);
            cfg.layout      = layout;
            cfg.showscale = 'yes';
            cfg.comment = 'no';
            cfg.colormap = co;
            cfg.style = 'straight';
            ft_topoplotER(cfg, ft);
            drawnow
%             a = gca;
%             text(mean(a.XLim),min(a.YLim-.1),sprintf('%ims',mean(timewins(ttt,:))),'HorizontalAlignment','center','FontSize',20)
            
%             if ttt==1
%                 tt = sprintf('%s',stc{n});
%                 t=title(sprintf('%s',tt),'Units','Normalized','Position',[-.065,1.02,1],'FontSize',20,'HorizontalAlignment','left')
%             end
        end
    end
end

%% annotate
annotation('textbox',[0.16 .9 .3 .1],...
    'String','Peripheral','FontSize',24,'LineStyle','none', ...
    'HorizontalAlignment','center')

annotation('textbox',[0.57 .9 .3 .1],...
    'String','Dual peripheral','FontSize',24,'LineStyle','none', ...
    'HorizontalAlignment','center')

annotation('textbox',[0.06 .865 .3 .1],...
    'String','Left visual field','FontSize',20,'LineStyle','none', ...
    'HorizontalAlignment','center')

annotation('textbox',[0.26 .865 .3 .1],...
    'String','Right visual field','FontSize',20,'LineStyle','none', ...
    'HorizontalAlignment','center')

annotation('textbox',[0.47 .865 .3 .1],...
    'String','Left visual field','FontSize',20,'LineStyle','none', ...
    'HorizontalAlignment','center')

annotation('textbox',[0.67 .865 .3 .1],...
    'String','Right visual field','FontSize',20,'LineStyle','none', ...
    'HorizontalAlignment','center')

for tw = 1:length(timewins)
    annotation('textbox',[0.02 .805-(tw-1)*.22 .1 .05],...
        'String',sprintf('%d ms',mean(timewins(tw,:))),'FontSize',20,'LineStyle','none', ...
        'HorizontalAlignment','center')
end

%% save plot
set(gcf,'PaperOrientation','portrait');
fn = 'figures/FigureS1_ch_searchlight_topo_peripheral';
tn = tempname;
print(gcf,'-dpng','-r500',tn)
im=imread([tn '.png']);
[i,j]=find(mean(im,3)<255);margin=2;
imwrite(im(min(i-margin):max(i+margin),min(j-margin):max(j+margin),:),[fn '.png'],'png');

