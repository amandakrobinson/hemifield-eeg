%% tabulate onsets and peaks
load('./results/stats_decoding_peripheral.mat')

conds = {'Peripheral' 'MultiplePeripheral'};
vf = {'LVF' 'RVF'};
hemi = {'L_elecs' 'R_elecs'};
heminames = {'Left' 'Right'};

%% get onsets and peaks for each cond/hemi/VF combination

allinfo=struct();
x=0;
for c = 1:length(conds)
    for h = 1:length(hemi)
        for v = 1:length(vf)

            s = stats.(conds{c}).(hemi{h}).(vf{v});

            x=x+1;
            allinfo.Condition{x,1} = conds{c};
            allinfo.Hemisphere{x,1} = heminames{h};
            allinfo.ImageSide{x,1} = vf{v};

            if h == v
                allinfo.HemiXvf{x,1} = 'ipsi';
            else
                allinfo.HemiXvf{x,1} = 'contra';
            end
            allinfo.Onset(x,1) = s.onset;
            allinfo.Onset_lowCI(x,1) = s.onsetci(1);
            allinfo.Onset_highCI(x,1) = s.onsetci(2);
            allinfo.Peak(x,1) = s.peak;
            allinfo.Peak_lowCI(x,1) = s.peakci(1);
            allinfo.Peak_highCI(x,1) = s.peakci(2);
            allinfo.PeakAcc(x,1) = max(s.mu);
        end
    end
end

allinfo = struct2table(allinfo);

save('results/onsets_peaks.mat','allinfo');

%% plot
cmap = viridis(3);

f=figure(1);clf
set(f,'Position',[1657 150 1248 643])

for c = 1:length(conds)
    sub=subplot(2,1,c);
    set(gca,'FontSize',20)
    set(gca,'linewidth',2)
    hold on
    ylim([0 5])
    xlim([0 300])

    for h = 1:length(hemi)
        for ci = 1:2
            
            if ci == 1 % contra
                v = mod(h,2)+1;
            else
                v = h;
            end

            s = stats.(conds{c}).(hemi{h}).(vf{v});
            ypos= 5-((h-1)*2+ci);
            col = cmap(ci,:);
            plot(s.onsetci,[ypos ypos],'LineWidth',3,'Color',col)
            plot(s.onset,ypos,'o','MarkerSize',20,'Color',col)

            plot(s.peakci,[ypos ypos],'LineWidth',3,'Color',col)
            plot(s.peak,ypos,'x','MarkerSize',20,'Color',col)

        end
    end
        
    sub.YTick = 1:4;
    sub.YTickLabel = {'Right hemisphere: ipsi'  'Right hemisphere: contra'  'Left hemisphere: ipsi' 'Left hemisphere: contra'};
    title(['Decoding time: ' conds{c}])
    
    xlabel('Time (ms)');

end
