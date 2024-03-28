
function plot_withinVacross()

%% plot common/shared information across hemispheres compared to reliability within each hemisphere
% correlations in information (RSA) for same hemisphere (split-half) versus
% correlations in information between the left and right hemispheres

%% load results
load('results/stats_within.mat','stats','timevect');

%% set variables
conds = 'Peripheral';
vfs = {'LVF' 'RVF'};

% plot details
xlims = [-100 800];
yplotlim = [-.04 0.2];

%% plot within vs across hemisphere correlations
% initialise figure
fig=figure(1);clf;
fig.Resize = 'off';
fig.Position=[1 1 1500 600];

comps = {'Within LH' 'Within RH' 'Across LH-RH'};
cmap = viridis(3);

% plot per condition
for v = 1:length(vfs)
    % plot correlations
    sub=subplot(2,length(vfs),v);
    set(gca,'FontSize',20)
    set(gca,'linewidth',2)
    sub.YTick = 0:.1:1;
    sub.TickDir = 'out';
    sub.XTick = 0:200:1000;

    hold on
    plot(timevect,timevect*0,'k','HandleVisibility','off')
    plot(yplotlim*0,yplotlim,'Color',[.7 .7 .7],'HandleVisibility','off') % vertical line at 0ms

    % get within/across dat
    con = {'L_elecs','R_elecs','across'};
    dat_vf = stats.(conds).(vfs{v});
    for b = 1:length(con)

        dat = dat_vf.(con{b});
        col = cmap(b,:);

        mu = dat.mu';
        se = dat.se';

        fill([timevect fliplr(timevect)],[mu-se fliplr(mu+se)],col,'FaceAlpha',0.2,'LineStyle','none','HandleVisibility','off')
        plot(timevect,mu,'Color',col,'LineWidth',2)

    end

    ylim(yplotlim);
    xlim(xlims);
    xlabel('Time (ms)')
    ylabel('Correlation')

    legend(comps,'Box','off')

    %% plot bayes factors

    for b = 1:3

        bf = dat_vf.(con{b}).bf;

        col = cmap(b,:);

        sp=subplot(6,length(vfs),length(vfs)*3+v+(b-1)*length(vfs));
        set(sp,'FontSize',20);
        set(sp,'YScale','log')
        hold on

        % plot BFs
        plot(timevect,1+0*timevect,'k-');
        cob = [.5 .5 .5;1 1 1;col];
        idx = [bf<1/10,1/10<bf & bf<10,bf>10]';
        for i=1:3
            x = timevect(idx(i,:));
            y = bf(idx(i,:));
            if ~isempty(x)
                stem(x,y,'Marker','o','Color',.6*[1 1 1],'BaseValue',1,'MarkerSize',5,'MarkerFaceColor',cob(i,:),'Clipping','off');
                plot(x,y,'o','Color',.6*[1 1 1],'MarkerSize',5,'MarkerFaceColor',cob(i,:),'Clipping','off');
            end
        end
        sp.YScale='log';
        ylim(10.^(6.5*[-1 1]))
        sp.YTick = 10.^([-5 0 5]);
        xlim([min(timevect) max(timevect)+1])
        ylabel('BF')
        if ismember(b,[1 2])
            set(sp,'XTickLabel','');
        end

        % add labels
        % bflab = sprintf('%s within > across',bflabels{b});
        text(300,1e4,comps{b},'VerticalAlignment','middle','HorizontalAlignment','left','FontSize',20,'Color',col);

    end
    xlabel('Time (ms)')
end
%% save
fn = 'figures/FigureSX_withinVacross_perVF';
print(gcf,'-dpng','-r500',fn)
im=imread([fn '.png']);
[i,j]=find(mean(im,3)<255);margin=2;
imwrite(im(min(i-margin):max(i+margin),min(j-margin):max(j+margin),:),[fn '.png'],'png');


%% plot within/across for contra/ipsi
% initialise figure
fig=figure(2);clf;
fig.Resize = 'off';
fig.Position=[1 1 800 900];

comps = {'Within contralateral' 'Within ipsilateral' 'Across contra-ipsi'};
bfcomps = {'Within versus across: contra' 'Within versus across: ipsi'};

% plot correlations
sub=subplot(2,1,1);
set(gca,'FontSize',20)
set(gca,'linewidth',2)
sub.YTick = 0:.1:1;
sub.TickDir = 'out';
sub.XTick = 0:200:1000;

hold on
plot(timevect,timevect*0,'k','HandleVisibility','off')
plot(yplotlim*0,yplotlim,'Color',[.7 .7 .7],'HandleVisibility','off') % vertical line at 0ms

% get dat
cont = {'contra' 'ipsi' 'across'};
cmap = tab10(3);

for x = 1:3
    dat = stats.(conds).(cont{x});

    col = cmap(x,:);

    mu = dat.mu';
    se = dat.se';

    fill([timevect fliplr(timevect)],[mu-se fliplr(mu+se)],col,'FaceAlpha',0.2,'LineStyle','none','HandleVisibility','off')
    plot(timevect,mu,'Color',col,'LineWidth',2)

end

ylim(yplotlim);
xlim(xlims);
xlabel('Time (ms)')
ylabel('Correlation')

legend(comps,'Box','off')

% plot bayes factors

for b = 1:2

    bf = stats.(conds).withinVacross.(cont{b}).bf;

    col = cmap(b,:);

    sp=subplot(8,1,4+b);
    set(sp,'FontSize',20);
    set(sp,'YScale','log')
    hold on

    % plot BFs
    plot(timevect,1+0*timevect,'k-');
    cob = [.5 .5 .5;1 1 1;col];
    idx = [bf<1/10,1/10<bf & bf<10,bf>10]';
    for i=1:3
        x = timevect(idx(i,:));
        y = bf(idx(i,:));
        if ~isempty(x)
            stem(x,y,'Marker','o','Color',.6*[1 1 1],'BaseValue',1,'MarkerSize',5,'MarkerFaceColor',cob(i,:),'Clipping','off');
            plot(x,y,'o','Color',.6*[1 1 1],'MarkerSize',5,'MarkerFaceColor',cob(i,:),'Clipping','off');
        end
    end
    sp.YScale='log';
    ylim(10.^(6.5*[-1 1]))
    sp.YTick = 10.^([-5 0 5]);
    xlim([min(timevect) max(timevect)+1])
    ylabel('BF')
    if ismember(b,[1 2])
        set(sp,'XTickLabel','');
    end

    % add labels
    % bflab = sprintf('%s within > across',bflabels{b});
    text(300,1e4,bfcomps{b},'VerticalAlignment','middle','HorizontalAlignment','left','FontSize',20,'Color',col);

    if b == 3
        xlabel('Time (ms)')
    end
end

%% save
fn = 'figures/Figure6_withinVacross';
print(gcf,'-dpng','-r500',fn)
im=imread([fn '.png']);
[i,j]=find(mean(im,3)<255);margin=2;
imwrite(im(min(i-margin):max(i+margin),min(j-margin):max(j+margin),:),[fn '.png'],'png');

end

