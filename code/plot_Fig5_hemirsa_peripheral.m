function plot_Fig5_hemirsa_peripheral

%% plot Figure 5: shared structure of representations across the hemispheres

load('results/stats_across_hemi_correlations.mat','stats','combinedstats','timevect')

conds = {'Peripheral' 'MultiplePeripheral'};

%% first plot LH/RH neural RDMs
fig=figure(1);clf;
fig.Position=[1000 1 2000 1500];
fig.Resize = 'off';

cmap = brewermap(50,'PuOr');
cmap([1:5 46:50],:) = [];


load('results/rdms.mat','rdm')

t1 = find(timevect==170);
t2 = find(timevect==200);
rdm1dat = squareform(rdm.periph_RVFcontra(:,t1));
rdm2dat = squareform(rdm.periph_RVFipsi(:,t2));

% plot example LH RDM
a = axes('Units','pixels','Position',[350 1100 200 200],'Visible','off');
imagesc(rdm1dat,[.5 .58]);
a.YDir='normal';
axis square
set(gca,'Colormap',viridis)
axis off

annotation('textbox','Units','Pixels','Position', ...
    [350 1310 200 20],...
    'LineStyle','none','String','Left hemisphere',...
    'FontSize',20,'VerticalAlignment','middle','HorizontalAlignment','center');

annotation('textbox','Units','Pixels','Position', ...
    [350 1070 200 20],...
    'LineStyle','none','String','170 ms',...
    'FontSize',20,'VerticalAlignment','middle','HorizontalAlignment','center');

% plot example RH RDM
a = axes('Units','pixels','Position',[750 1100 200 200],'Visible','off');
imagesc(rdm2dat,[.5 .58]);
a.YDir='normal';
axis square
set(gca,'Colormap',viridis)
axis off

annotation('textbox','Units','Pixels','Position', ...
    [750 1310 200 20],...
    'LineStyle','none','String','Right hemisphere',...
    'FontSize',20,'VerticalAlignment','middle','HorizontalAlignment','center');

annotation('textbox','Units','Pixels','Position', ...
    [750 1070 200 20],...
    'LineStyle','none','String','200 ms',...
    'FontSize',20,'VerticalAlignment','middle','HorizontalAlignment','center');

annotation('arrow',[.29 .36], [.8 .8],'LineWidth',2)
annotation('arrow',[.36 .29], [.8 .8],'LineWidth',2)
annotation('textbox','Units','Pixels','Position', ...
    [550 1200 200 20],...
    'LineStyle','none','String','Correlate',...
    'FontSize',16,'VerticalAlignment','middle','HorizontalAlignment','center');

annotation('arrow',[.49 .56], [.8 .8],'LineWidth',2)
annotation('textbox','Units','Pixels','Position', ...
    [980 1210 150 20],...
    'LineStyle','none','String','Repeat for all pairs of time points',...
    'FontSize',16,'VerticalAlignment','middle','HorizontalAlignment','center');

% now plot timeXtime correlation matrix
ttc = mean(stats.data.Peripheral.RVF,3);
a = axes('Units','pixels','Position',[1200 1070 260 260],'Visible','off');
imagesc(timevect,timevect,ttc,[-.04 .04]);
a.YDir='normal';
axis square
set(gca,'Colormap',cmap)
cb = colorbar();
ylabel(cb,'Spearman correlation','FontSize',18,'Rotation',90)

set(gca,'FontSize',18)
ylabel('LH time (ms)')
xlabel('RH time (ms)')
hold on
xlim([-100 800])
ylim([-100 800])
% plot(timevect,timevect,'k')
plot([0 0],[-100 800],'--k')
plot([-100 800],[-100 800],'k')
plot([-100 800],[0 0],'--k')
plot([200],[170],'r','Marker','o','MarkerSize',15)

annotation('arrow',[.64 .73], [.77 .715],'LineWidth',2,'Color','r')
annotation('textbox','Units','Pixels','Position', ...
    [1460 1050 155 20],...
    'LineStyle','none','String','Correlation between LH at 170 ms and RH at 200 ms',...
    'FontSize',14,'VerticalAlignment','middle','HorizontalAlignment','center');

% A/B labels
annotation('textbox','Units','Pixels','Position', ...
    [300 1330 20 20],...
    'LineStyle','none','String','A',...
    'FontSize',40,'VerticalAlignment','middle','HorizontalAlignment','center');
annotation('textbox','Units','Pixels','Position', ...
    [300 980 20 20],...
    'LineStyle','none','String','B',...
    'FontSize',40,'VerticalAlignment','middle','HorizontalAlignment','center');
annotation('textbox','Units','Pixels','Position', ...
    [300 530 20 20],...
    'LineStyle','none','String','C',...
    'FontSize',40,'VerticalAlignment','middle','HorizontalAlignment','center');
annotation('textbox','Units','Pixels','Position', ...
    [300 250 20 20],...
    'LineStyle','none','String','D',...
    'FontSize',40,'VerticalAlignment','middle','HorizontalAlignment','center');

%% plot timegen correlations
lims = [-0.03 .03];

for c = 1:length(conds)

    dat = combinedstats.data.(conds{c}).mu;
    bf = combinedstats.data.(conds{c}).bf;
    dat(bf<3) = 0; % threshold for bf>3

    % plot time gen RSA correlations
    a = axes('Units','pixels','Position',[(c-1)*550+430 600 400 400],'Visible','off');
    imagesc(timevect,timevect,dat,lims)
    axis square
    hold on
    plot([0 0],[-100 800],'--k','LineWidth',2)
    plot([-100 800],[-100 800],'k','LineWidth',2)
    plot([-100 800],[0 0],'--k','LineWidth',2)

    ylabel('Contralateral hemisphere time (ms)')
    xlabel('Ipsilateral hemisphere time (ms)')
    ylim([-100 800])
    xlim([-100 800])
    colormap(gca,cmap)
    set(gca,'YDir','normal')
    cb=colorbar();
    cb.Ticks = [lims(1) 0 lims(2)];
    if c  == 2
        cb.Label.String = 'Spearman correlation';
    end
    set(gca,'FontSize',20)

end

annotation('textbox','Units','Pixels','Position', ...
    [400 980 400 20],...
    'LineStyle','none','String','Single Peripheral',...
    'FontSize',30,'VerticalAlignment','middle','HorizontalAlignment','center');

annotation('textbox','Units','Pixels','Position', ...
    [950 980 400 20],...
    'LineStyle','none','String','Dual Peripheral',...
    'FontSize',30,'VerticalAlignment','middle','HorizontalAlignment','center');

%% for each BF>3, what is the contra-ipsi offset?
% figure('Position',[297 206 1240 707])
cmap = tab10(3);

for c = 1:length(conds)
    
    mu = combinedstats.data.(conds{c}).mu;
    bf = combinedstats.data.(conds{c}).bf;

    % get times of reliable correlations
    [time_contra,time_ipsi] = find(mu>0&bf>3);

    % restrict to times 0-500ms
    idx = find(timevect(time_contra)<0|timevect(time_contra)>500|timevect(time_ipsi)<0|timevect(time_ipsi)>500);
    time_contra(idx)=[];
    time_ipsi(idx) = [];

    % plot times of contra/ipsi correlations
    axes('Units','pixels','Position',[(c-1)*550+430 340 380 200],'Visible','off');
    histogram(timevect(time_contra),-100:5:800,'FaceColor',cmap(1,:),'EdgeColor','none')
    hold on
    histogram(timevect(time_ipsi),-100:5:800,'FaceColor',cmap(2,:),'EdgeColor','none')
    xlim([0 500])
    ylim([0 800])
    set(gca,'FontSize',18)
    legend({'Contralateral' 'Ipsilateral'},'Box','off')
    xlabel('Time of reliable interhemispheric correlation (ms)')
    ylabel('Number of time points')

    % plot ipsi-contra time difference
    axes('Units','pixels','Position',[(c-1)*550+430 60 380 200],'Visible','off');
    histogram(timevect(time_ipsi)-timevect(time_contra),-1000:5:1000,'FaceColor','k','EdgeColor','none')
    xlim([-250 250])
    ylim([0 550])
    set(gca,'FontSize',18)
    legend({'Ipsi-contra difference'},'Box','off','Location','northwest')
    hold on
    plot([0 0],[0 550],'k','HandleVisibility','off');
    xlabel('Ipsilateral-Contralateral time difference (ms)')
    ylabel('Number of time points')

end

%% save

fn = 'figures/Figure5_hemisphere_correlations_peripheral';
print(gcf,'-dpng','-r500',fn)
im=imread([fn '.png']);
[i,j]=find(mean(im,3)<255);margin=2;
imwrite(im(min(i-margin):max(i+margin),min(j-margin):max(j+margin),:),[fn '.png'],'png');



