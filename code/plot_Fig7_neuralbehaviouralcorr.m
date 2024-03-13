function plot_Fig7_neuralbehaviouralcorr()
%% plot figure 7 of peripheral manuscript
% plot correlations between neural RDMs and behavioural task results
% collapsed across left and right hemispheres (contra v ipsi only)
load('results/stats_neuralbehavcorr.mat','stats','timevect')

%% plotting peripheral

lims = [-.04 .2];

fig=figure(1);clf;
fig.Position=[1 1 2000 1500];
fig.Resize = 'off';

cmap = tab10(6);
cmap(1:3,:)=[];
ci = {'contra' 'ipsi'};
pn = 13:14;

for plotnr = 1:length(ci)

    %% plot correlations
    sa = subplot(12,2,pn(plotnr):2:(pn(plotnr)+4));
    set(gca,'Fontsize',20)
    if plotnr<3
        sa.Position = [sa.Position(1),sa.Position(2)+.05,sa.Position(3:4)];
    else
         sa.Position = [sa.Position(1),sa.Position(2)-.04,sa.Position(3:4)];
    end
    hold on
    plot(timevect,timevect*0,'Color',[.7 .7 .7],'HandleVisibility','off')
    plot(lims*0,lims,'Color',[.7 .7 .7],'HandleVisibility','off') % vertical line at 0ms

    dat = stats.Peripheral.(ci{plotnr});
    for m = 1:size(dat.mu,2) % for each model
        col = cmap(m,:);
        mu = dat.mu(:,m);
        se = dat.se(:,m);

        fill([timevect fliplr(timevect)],[mu'-se' fliplr(mu'+se')],col,'FaceAlpha',0.2,'LineStyle','none','HandleVisibility','off')
        plot(timevect,mu,'Color',col,'Linewidth',2)
    end
    xlim([min(timevect) max(timevect)+1])
    legend({'Image task' 'Concept task'},'Box','off')
    % title('Left hemisphere')
    ylim(lims)
    ylabel('Neural-behaviour correlation')
    xlabel('Time (ms)')

    %% plot BFs
    for m = 1:3 % two models then differences
        
        % colour of line and figure position
        if m<3
            col = cmap(m,:);
            bf = dat.bf(:,m);
        else 
            col = [0 0 0];
            bf = dat.moddiff.bf;
        end
        if plotnr<3
            sp=subplot(12,2,pn(plotnr)+6+(m-1)*2);
            sp.Position = [sp.Position(1),sp.Position(2)+.02,sp.Position(3:4)];
        else
            sp=subplot('Position',[sa.Position(1),sa.Position(2)-.03-.07*m,sp.Position(3:4)]);
        end

        set(gca,'Fontsize',20)
        hold on

        % plot
        plot(timevect,1+0*timevect,'k-');
        cob = [.5 .5 .5;1 1 1;col];
        idx = [bf<1/10,1/10<bf & bf<10,bf>10]';
        for i=1:3
            x = timevect(idx(i,:));
            y = bf(idx(i,:));
            if ~isempty(x)
                stem(x,y,'Marker','o','Color',.6*[1 1 1],'BaseValue',1,'MarkerSize',4,'MarkerFaceColor',cob(i,:),'Clipping','off');
                plot(x,y,'o','Color',.6*[1 1 1],'MarkerSize',4,'MarkerFaceColor',cob(i,:),'Clipping','off');
            end
        end
        sp.YScale='log';
        ylim(10.^(6*[-1 1]))
        sp.YTick = 10.^([-5 0 5]);
        xlim([min(timevect) max(timevect)+1])
        ylabel('BF')
        sp.XTickLabel = '';
        if m == 3
            xlabel('Time (ms)')
        end
    end
end

annotation('textbox',[0.14 .505 .3 .1],...
    'String','Contralateral','FontSize',24,'LineStyle','none', ...
    'HorizontalAlignment','center')

annotation('textbox',[0.585 .505 .3 .1],...
    'String','Ipsilateral','FontSize',24,'LineStyle','none', ...
    'HorizontalAlignment','center')


%% plot neural and behaviour RDMs

load('results/rdms.mat','rdm')

t = find(timevect==170);
rdmdat = squareform(rdm.periph_LVFcontra(:,t));

a = axes('Units','pixels','Position',[350 1100 200 200],'Visible','off');
imagesc(rdmdat,[.5 .6]);
a.YDir='normal';
axis square
colormap viridis
axis off

annotation('textbox','Units','Pixels','Position', ...
    [350 1310 200 20],...
    'LineStyle','none','String','Neural dissimilarity',...
    'FontSize',20,'VerticalAlignment','middle','HorizontalAlignment','center');

annotation('textbox','Units','Pixels','Position', ...
    [350 1070 200 20],...
    'LineStyle','none','String','170 ms',...
    'FontSize',20,'VerticalAlignment','middle','HorizontalAlignment','center');

% image task
a = axes('Units','pixels','Position',[800 1250 200 200],'Visible','off');
imagesc(squareform(stats.models.models(:,1)));
a.YDir='normal';
axis square
colormap viridis
axis off

% category task
a = axes('Units','pixels','Position',[800 950 200 200],'Visible','off');
imagesc(squareform(stats.models.models(:,2)));
a.YDir='normal';
axis square
colormap viridis
axis off

annotation('textbox','Units','Pixels','Position', ...
    [800 1460 200 20],...
    'LineStyle','none','String','Image task',...
    'FontSize',20,'VerticalAlignment','middle','HorizontalAlignment','center');

annotation('textbox','Units','Pixels','Position', ...
    [800 1160 200 20],...
    'LineStyle','none','String','Concept task',...
    'FontSize',20,'VerticalAlignment','middle','HorizontalAlignment','center');

annotation('arrow',[.28 .39], [.8 .9],'LineWidth',2)
annotation('arrow',[.28 .39], [.79 .69],'LineWidth',2)


annotation('textbox','Units','Pixels','Position', ...
    [600 1310 100 20],...
    'LineStyle','none','String','Spearman correlation',...
    'FontSize',16,'VerticalAlignment','middle','HorizontalAlignment','center');

annotation('textbox','Units','Pixels','Position', ...
    [600 1060 100 20],...
    'LineStyle','none','String','Spearman correlation',...
    'FontSize',16,'VerticalAlignment','middle','HorizontalAlignment','center');

annotation('arrow',[.53 .65], [.79 .79],'LineWidth',2)
annotation('textbox','Units','Pixels','Position', ...
    [1120 1210 130 20],...
    'LineStyle','none','String','Repeat for all time points',...
    'FontSize',16,'VerticalAlignment','middle','HorizontalAlignment','center');


% now plot example correlation
sa = axes('Units','pixels','Position',[1350 1100 400 200]);
set(gca,'Fontsize',20)
hold on
plot(timevect,timevect*0,'Color',[.7 .7 .7],'HandleVisibility','off')
plot(lims*0,lims,'Color',[.7 .7 .7],'HandleVisibility','off') % vertical line at 0ms
xlim([-100 800])
mu = stats.Peripheral.L_elecs.RVF.mu(:,1);
plot(timevect,mu,'Color','k','Linewidth',2)
xlabel('Time (ms)')

annotation('textbox','Units','Pixels','Position', ...
    [1350 1310 400 20],...
    'LineStyle','none','String','Neural-behaviour correlation',...
    'FontSize',20,'VerticalAlignment','middle','HorizontalAlignment','center');

annotation('textbox','Units','Pixels','Position', ...
    [140 1470 20 20],...
    'LineStyle','none','String','A',...
    'FontSize',40,'VerticalAlignment','middle','HorizontalAlignment','center');
annotation('textbox','Units','Pixels','Position', ...
    [140 900 20 20],...
    'LineStyle','none','String','B',...
    'FontSize',40,'VerticalAlignment','middle','HorizontalAlignment','center');

%% save
set(gcf,'PaperOrientation','portrait');
fn = 'figures/Figure7_behaviouralRSA';
tn = tempname;
print(gcf,'-dpng','-r500',tn)
im=imread([tn '.png']);
[i,j]=find(mean(im,3)<255);margin=2;
imwrite(im(min(i-margin):max(i+margin),min(j-margin):max(j+margin),:),[fn '.png'],'png');

