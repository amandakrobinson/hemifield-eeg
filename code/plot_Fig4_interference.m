function plot_Fig4_interference()

%% plot Figure 4: decoding of single peripheral image versus dual peripheral images

%% stuff for plotting
addpath('~/Dropbox/MATLAB/Colormaps/Colormaps (5)/Colormaps')
colormap viridis
cmap = tab10(3);

%% load results
load('results/stats_decoding_peripheral.mat','stats');
timevect = stats.timevect;

%% set up

conds = {'Peripheral' 'MultiplePeripheral'};

xlims = [-100 800];
yplotlim = [0.48 0.58];

%% plot peripheral, dual peripheral paired decoding results together

fig=figure(2);clf;
fig.Resize = 'off';
fig.Position=[1 1 800 900];

linetype = {'-','--'};
labels = {'Contralateral (single)' 'Ipsilateral (single)','Contralateral (dual)' 'Ipsilateral (dual)'};

% average across clusters
% get data
for c = 1:length(conds)
    
    sub=subplot(8,1,2:5);
    set(gca,'FontSize',25)
    set(gca,'linewidth',2)
    
    sub.YTick = .5:.02:.58;
    sub.TickDir = 'out';
    sub.XTick = 0:200:1000;
    
    hold on
    plot(timevect,timevect*0+0.5,'k','HandleVisibility','off')
    plot(yplotlim*0,yplotlim,'Color',[.7 .7 .7],'HandleVisibility','off') % vertical line at 0ms

    % draw two traces - contra and ipsi
    for ic = 1:2
        if ic == 1
            % contra
            dat1 = stats.(conds{c}).L_elecs.RVF.mu_all;
            dat2 = stats.(conds{c}).R_elecs.LVF.mu_all;
        else
            % contra
            dat1 = stats.(conds{c}).L_elecs.LVF.mu_all;
            dat2 = stats.(conds{c}).R_elecs.RVF.mu_all;
        end
        
        dat = mean(cat(3,dat1,dat2),3);
        
        alldat(:,:,c,ic) = dat;
        
        mu = mean(dat);
        n = size(dat,1);
        se = std(dat)./sqrt(n);
        
        % plot decoding
        col = cmap(ic,:);
        fill([timevect fliplr(timevect)],[mu-se fliplr(mu+se)],col,'FaceAlpha',0.2,'LineStyle','none','HandleVisibility','off')
        plot(timevect,mu,'Color',col,'LineWidth',2,'LineStyle',linetype{c})

    end
end
ylim(yplotlim);
xlim(xlims);
xlabel('Time (ms)')
ylabel('Decoding accuracy')

legend(labels,'Box','off')

labels = {'Single-dual difference: contralateral','Single-dual difference: ipsilateral'};
for ic = 1:2
    
    dat = alldat(:,:,1,ic)-alldat(:,:,2,ic); % ipsi/contra dat for both single and multiple conditions
    
    % calculate bayesfactors
    bf = bayesfactor_R_wrapper(dat',...
        'returnindex',2,'verbose',false,...
        'args','mu=0,rscale="medium",nullInterval=c(-0.5,0.5)');    
    
    % group mean onset & peak
    [~,ot] = max(movmean(bf>10,[0 9])==1); % onset 10 consecutive tp BF>10
    [~,peak] = max(mean(dat,1)); % peak
    onset = stats.timevect(ot);
    peak = stats.timevect(peak);

    % plot bayes
    figure(2)
    sp=subplot(8,1,5+ic);
    sp.Position(2) = sp.Position(2)-.01-.005*ic;
    set(sp,'FontSize',25);
    hold on
    col = cmap(ic,:);
    
    % plot bayes ----------------
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
    ylim(10.^(6.5*[-1 1]))
    sp.YTick = 10.^([-5 0 5]);
    sp.TickDir = 'out';

    xlim([min(timevect) max(timevect)+1])
    ylabel('BF')
    
    set(gca,'XMinorTick','on')
    if ic==1 % x axis label
        sp.XTickLabel = '';
    else
        xlabel('Time (ms)')
    end
    
    % add labels
    text(300,1e5,labels{ic},'VerticalAlignment','middle','HorizontalAlignment','left','FontSize',20,'Color',col);
   
    % now that we have onsets and peaks, add them to decoding plot
%     subplot(2,1,1)
%     % plot onset and peak CIs
%     yval = 0.5-.012-(ic-1)*.004;
% 
%     plot(onset,yval,'.','Color',col,'MarkerSize',12,'HandleVisibility','off')
%     plot(peak,yval,'-s','Color',col,'MarkerSize',12,'HandleVisibility','off')
%     %             line(dat.onsetci+[-.5 .5],[yval yval],'Color',co(t,:),'LineWidth',2)
%     %             line(dat.peakci+[-.5 .5],[yval yval],'Color',co(t,:),'LineWidth',2)
    drawnow
    
end

%% now draw images of different conditions at the top
im_single = imread('figures/screen_00419.png');
im_dual = imread('figures/screen_00007.png');

qs = round(prctile(1:1400,[15 85]));

axes('Units','pixels','Position',[200 700 150 200],'Visible','off');
imshow(im_single(qs(1):qs(2),:,:))

axes('Units','pixels','Position',[450 700 150 200],'Visible','off');
imshow(im_dual(qs(1):qs(2),:,:))

% add figure labels
annotation('textbox',[0.19 .88 .3 .1],...
    'String','Single peripheral','FontSize',20,'LineStyle','none', ...
    'HorizontalAlignment','center')

annotation('textbox',[0.505 .88 .3 .1],...
    'String','Dual peripheral','FontSize',20,'LineStyle','none', ...
    'HorizontalAlignment','center')

%% save
fn = 'figures/Figure4_decoding_interference';
print(gcf,'-dpng','-r500',fn)
im=imread([fn '.png']);
[i,j]=find(mean(im,3)<255);margin=2;
imwrite(im(min(i-margin):max(i+margin),min(j-margin):max(j+margin),:),[fn '.png'],'png');


