
function plot_FigS4_decoding_right_handers()

%% right handers only: plot peripheral decoding results for left and right hemispheres according to stimulus visual field
% FigS4 = single peripheral and dual peripheral

%% stuff for plotting
addpath('~/Dropbox/MATLAB/Colormaps/Colormaps (5)/Colormaps')
addpath('~/Dropbox/MATLAB/altmany-export_fig')
colormap viridis
co3=[];
cmap = viridis(3);
cmap(3,:) = [0 0 0];
for c = 1:3 % for each colour
    co3(:,:,c) = [linspace(1,cmap(c,1),5);...
        linspace(1,cmap(c,2),5);...
        linspace(1,cmap(c,3),5)]';
end
co3(1,:,:) = [];
co3 = flipud(co3); % shade x RGB x colour

%% load results
load('results/stats_decoding_peripheral.mat','stats');
timevect = stats.timevect;

%% get participant details
participants = readtable('data/participants.tsv',"FileType","text",'Delimiter', '\t');
righthanders = participants.Total>50;

%% set variables
conds = {'Peripheral' 'MultiplePeripheral'};
clustnames = {'Left hemisphere' 'Right hemisphere' 'Right-Left difference'};
contraipsi = [3 1; 1 3; 1 3]; % for changing colours from left/right to contra/ipsi

% bflabels = {'Contralateral','Ipsilateral'};

% plot details
xlims = [-100 800];
yplotlim = [0.48 0.58; 0.48 0.58; -.02 .08];
chance = [0.5 0.5 0];

%% plot peripheral, multipleperipheral paired decoding results
% initialise figure
fig=figure(1);clf;
fig.Resize = 'off';
fig.Position=[1 1 2000 900];

for cond = 1:length(conds) % peripheral, multi peripheral

    %% plot per condition
    dat_p = stats.(conds{cond});
    clust = fieldnames(dat_p); % LH, RH, LH-RH diff

    for cl = 1:length(clust)

        % get data to plot
        dat_cond = dat_p.(clust{cl});
        vfs = fieldnames(dat_cond);
        vfs = vfs(~ismember(vfs,{'LVFRVFdiff' 'conipdiff'})); % LVF/RVF or contra/ipsi

        sub=subplot(2,3,(cond-1)*3+cl);
        set(gca,'FontSize',25)
        set(gca,'linewidth',2)
        if cl<3
            sub.YTick = .5:.02:.58;
        else
            sub.YTick = -.02:.02:.1;
        end
        sub.TickDir = 'out';
        sub.XTick = 0:200:1000;

        hold on
        plot(timevect,timevect*0+chance(cl),'k','HandleVisibility','off')
        plot(yplotlim(cl,:)*0,yplotlim(cl,:),'Color',[.7 .7 .7],'HandleVisibility','off') % vertical line at 0ms

        for v = 1:length(vfs) % separate line for each VF

            col = co3(contraipsi(cl,v),:,cl);

            dat = dat_cond.(vfs{v}).mu_all(righthanders,:);
            mu = mean(dat,1);
            se = std(dat,[],1)/sqrt(size(dat,1));

            fill([timevect fliplr(timevect)],[mu-se fliplr(mu+se)],col,'FaceAlpha',0.2,'LineStyle','none','HandleVisibility','off')
            p(v) = plot(timevect,mu,'Color',col,'LineWidth',2);

            drawnow

        end

        ylim(yplotlim(cl,:));
        xlim(xlims);
        title(clustnames{cl})
        xlabel('Time (ms)')
        ylabel('Decoding accuracy')

        if cl==1
            lh=legend([p(2),p(1)],["Contralateral stimulus", "Ipsilateral stimulus"],'Box','off');
        else
            legend({'Contralateral stimulus' 'Ipsilateral stimulus'},'Box','off')
        end

    end
end

% add figure labels
annotation('textbox',[0.08 .89 .1 .1],...
    'String','A','FontSize',30,'LineStyle','none')

annotation('textbox',[0.36 .89 .1 .1],...
    'String','B','FontSize',30,'LineStyle','none')

annotation('textbox',[0.635 .89 .1 .1],...
    'String','C','FontSize',30,'LineStyle','none')

annotation('textbox',[0.08 .42 .1 .1],...
    'String','D','FontSize',30,'LineStyle','none')

annotation('textbox',[0.36 .42 .1 .1],...
    'String','E','FontSize',30,'LineStyle','none')

annotation('textbox',[0.635 .42 .1 .1],...
    'String','F','FontSize',30,'LineStyle','none')


%% save
fn = sprintf('figures/FigureS4_decoding_righthanders');
print(gcf,'-dpng','-r500',fn)
im=imread([fn '.png']);
[i,j]=find(mean(im,3)<255);margin=2;
imwrite(im(min(i-margin):max(i+margin),min(j-margin):max(j+margin),:),[fn '.png'],'png');