
% plot mean decoding accuracy for peripheral stimuli

function plot_Fig2and3_decoding_peripheral()

%% plot peripheral decoding results for left and right hemispheres according to stimulus visual field
% figure 2 = single peripheral
% figure 3 = dual peripheral

%% stuff for plotting
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

%% set variables
conds = {'Peripheral' 'MultiplePeripheral'};
clustnames = {'Left hemisphere' 'Right hemisphere' 'Right-Left difference'};
% clust = {'L_elecs' 'R_elecs'};
contraipsi = [3 1; 1 3; 1 3]; % for changing colours from left/right to contra/ipsi

bflabels = {'Contralateral','Ipsilateral'};

% plot details
xlims = [-100 800];
yplotlim = [0.48 0.58; 0.48 0.58; -.02 .08];
chance = [0.5 0.5 0];

%% plot peripheral, multipleperipheral paired decoding results

for cond = 1:length(conds) % peripheral, multi peripheral
    
    %% initialise figure
    fig=figure(cond);clf;
    fig.Resize = 'off';
    fig.Position=[1 1 2000 900];
    
    %% plot per condition
    dat_p = stats.(conds{cond});
    clust = fieldnames(dat_p); % LH, RH, LH-RH diff
    
    for cl = 1:length(clust)
        
        % get data to plot
        dat_cond = dat_p.(clust{cl});
        vfs = fieldnames(dat_cond);
        vfs = vfs(~ismember(vfs,{'LVFRVFdiff' 'conipdiff'})); % LVF/RVF or contra/ipsi
        
        sub=subplot(2,3,cl);
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
            
            dat = dat_cond.(vfs{v});
            mu = dat.mu;
            se = dat.se;
            
            fill([timevect fliplr(timevect)],[mu-se fliplr(mu+se)],col,'FaceAlpha',0.2,'LineStyle','none','HandleVisibility','off')
            p(v) = plot(timevect,mu,'Color',col,'LineWidth',2);
            
            % % plot onset and peak CIs
            % yval = chance(cl)-.01-contraipsi(cl,v)*.002;
            % if ~isempty(dat.onset)
            %     plot(dat.onset,yval,'.','Color',col,'MarkerSize',12,'HandleVisibility','off')
            % end
            % plot(dat.peak,yval,'-s','Color',col,'MarkerSize',12,'HandleVisibility','off')
            % 
            % if contraipsi(cl,v) == 1 % top onset point
            %     text(dat.onset-90,yval-.002,'Onset')
            % else
            %     text(dat.peak+30,yval+.002,'Peak')
            % end
%             if ~isempty(dat.onsetci)
%                 line(dat.onsetci+[-.5 .5],[yval yval],'Color',co(t,:),'LineWidth',2)
%             end
%             line(dat.peakci+[-.5 .5],[yval yval],'Color',co(t,:),'LineWidth',2)
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
        
        % plot bayes factors
        for v = 1:length(vfs) % new plot for each
            
            % colour of line and figure position (contra is upper)
            col = co3(contraipsi(cl,v),:,cl);
            if cl == 1 % left cluster
                ypos = mod(v,2)+1;
            else
                ypos = v;
            end
            
            bf = dat_cond.(vfs{v}).bf;
            
            sp=subplot(6,3,(3+ypos-1)*3+cl); % contra is upper plot
            set(sp,'FontSize',25);
            hold on
            
            % new bits ----------------
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
            xlim([min(timevect) max(timevect)+1])
            ylabel('BF')
            sp.XTickLabel = '';
            
            % add labels
            if cl == 1
                bflab = bflabels{mod(v,2)+1};
            else
                bflab = bflabels{v};
            end
            
            text(500,1.9e4,bflab,'VerticalAlignment','middle','HorizontalAlignment','left','FontSize',25,'Color',col);
            
        end
        
        %% now do difference bayes factors
        if cl<3
            bf = dat_cond.LVFRVFdiff.bf;
        else
            bf = dat_cond.conipdiff.bf;
        end
        
        sp=subplot(6,3,15+cl);
        set(sp,'FontSize',25);
        set(gca,'YScale','log')
        hold on
        
        plot(timevect,1+0*timevect,'k-');
        cob = [.5 .5 .5;1 1 1; 0 0 0];
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
        xlim([min(timevect) max(timevect)+1])
        ylabel('BF')
        
        xlabel('Time (ms)')
        
        % add labels
        text(500,1e5,'Difference','VerticalAlignment','middle','HorizontalAlignment','left','FontSize',25,'Color','k');
        
    end
    
    % add figure labels
    annotation('textbox',[0.08 .89 .1 .1],...
        'String','A','FontSize',30,'LineStyle','none')
    
    annotation('textbox',[0.36 .89 .1 .1],...
        'String','B','FontSize',30,'LineStyle','none')

    annotation('textbox',[0.635 .89 .1 .1],...
        'String','C','FontSize',30,'LineStyle','none')

    
    fn = sprintf('figures/Figure%d_decoding_%s',cond+1,conds{cond});
    print(gcf,'-dpng','-r500',fn)
    im=imread([fn '.png']);
    [i,j]=find(mean(im,3)<255);margin=2;
    imwrite(im(min(i-margin):max(i+margin),min(j-margin):max(j+margin),:),[fn '.png'],'png');
end

end

