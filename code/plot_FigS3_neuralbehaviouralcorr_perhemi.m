function plot_FigS3_neuralbehaviouralcorr_perhemi()

%% plot figure S3: neural-behavioural correlations per hemisphere 
% plot correlations between neural RDMs and behavioural task results

load('results/stats_neuralbehavcorr.mat','stats','timevect')

%% plotting peripheral

lims = [-.04 .2];

fig=figure(1);clf;
fig.Position=[1 1 2000 1500];
fig.Resize = 'off';
    
cmap = tab10(4);
cmap([2 3],:)=[];

clusts = {'L_elecs','R_elecs','L_elecs','R_elecs'}; % order for plotting
vfs = {'RVF' 'LVF' 'LVF' 'RVF'}; % order for plotting
pn = [1 2 13 14]; % subplot number

for plotnr = 1:length(clusts)

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

    dat = stats.Peripheral.(clusts{plotnr}).(vfs{plotnr});
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
%             sp=subplot(12,2,pn(plotnr)+6+(m-1)*2,'Position',[sa.Position(1),sa.Position(2)-.3-.02*m,sp.Position(3:4)]);
            sp=subplot('Position',[sa.Position(1),sa.Position(2)-.03-.07*m,sp.Position(3:4)]);

%              sp.Position = [sp.Position(1),sp.Position(2)-.1-.02*m,sp.Position(3:4)];
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

annotation('textbox',[0.14 .905 .3 .1],...
    'String','Left hemisphere','FontSize',24,'LineStyle','none', ...
    'HorizontalAlignment','center')

annotation('textbox',[0.585 .905 .3 .1],...
    'String','Right hemisphere','FontSize',24,'LineStyle','none', ...
    'HorizontalAlignment','center')

annotation('textbox',[0.02 .8 .05 .1],...
    'String','Contralateral stimuli','FontSize',24,'LineStyle','none', ...
    'HorizontalAlignment','center')

annotation('textbox',[0.02 .3 .05 .1],...
    'String','Ipsilateral stimuli','FontSize',24,'LineStyle','none', ...
    'HorizontalAlignment','center')


saveas(gcf,'figures/FigureS3_neuralbehaviour_correlation_perhemi.png')