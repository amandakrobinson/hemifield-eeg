
%% plots for figure 7: behaviour
% A - triplet task design with questions for two tasks 
% B - concept and image/word models
% C - behavioural RDMs for two tasks
% D - embedding for two tasks

% get stimuli and behavioural data
tasks = {'image_similarity' 'concept_similarity'};

load('results/imageorder.mat','imnames') % column 2 is stimuli names
stims = imnames(2,:);
stimorder = {'fish','bird','face','boats','tree','tools'};


f=figure('Resize','off');clf;
f.Position = [f.Position(1:2) 1000 1500];

%% make category model RDMs for B

% higher level category model - fish/bird/face/boats/tree/tools
m_cat = ones(36,36);
for c = 1:6
    init = (c-1)*4+1;
    widx = 25+(c-1)*2;
    idx = [init:(init+3) widx widx+1];
    m_cat(idx,idx) = 0;
end

% image vs word model
m_imw = ones(36,36);
m_imw(1:24,1:24) = 0;
m_imw(25:36,25:36) = 0;


%% A: plot triplet task design

% plot A) task design
clf
axes('Units','pixels','Position',[105 965 250 250],'Visible','off');
im = imread('figures/behavioural_design.png');
imshow(im)

annotation('textbox','Units','Pixels','Position', ...
    [5 935 450 20],...
    'LineStyle','none','String','"Choose the one that looks different"',...
    'FontSize',20,'VerticalAlignment','middle','HorizontalAlignment','center');
annotation('textbox','Units','Pixels','Position', ...
    [5 910 450 20],...
    'LineStyle','none','String','(Image task)',...
    'FontSize',20,'VerticalAlignment','middle','HorizontalAlignment','center');
annotation('textbox','Units','Pixels','Position', ...
    [5 865 450 20],...
    'LineStyle','none','String','"Choose the one that is conceptually different"',...
    'FontSize',20,'VerticalAlignment','middle','HorizontalAlignment','center');
annotation('textbox','Units','Pixels','Position', ...
    [5 840 450 20],...
    'LineStyle','none','String','(Concept task)',...
    'FontSize',20,'VerticalAlignment','middle','HorizontalAlignment','center');


%% plot B: category models
% im/word model
a = axes('Units','pixels','Position',[520 950 200 200],'Visible','off');
imagesc(m_imw);
a.YDir='normal';
axis square
colormap viridis
axis off

% plot category model
a = axes('Units','pixels','Position',[770 950 200 200],'Visible','off');
imagesc(m_cat);
a.YDir='normal';
axis square
colormap viridis
axis off

annotation('textbox','Units','Pixels','Position', ...
    [520 1165 200 20],...
    'LineStyle','none','String','Image/word model',...
    'FontSize',20,'VerticalAlignment','middle','HorizontalAlignment','center');

annotation('textbox','Units','Pixels','Position', ...
    [770 1165 200 20],...
    'LineStyle','none','String','Concept model',...
    'FontSize',20,'VerticalAlignment','middle','HorizontalAlignment','center');


%% plot RDMs and embeddings

medcat = [2.5:4:24 25.5:2:36];

for t = 1:length(tasks)
    load(sprintf('results/stats_%s.mat',tasks{t}),'RDMmean')

    % RDM
    a = axes('Units','pixels','Position',[110+(t-1)*450 450 350 350],'Visible','off');
    a.FontSize = 18;
    imagesc(RDMmean,[0 1]);axis square
    a.YDir='normal';
    colorbar
    colormap viridis
%     axis off
    set(a,'Fontsize',18)
    set(a,'XTick',medcat)
    set(a,'XTickLabels',[stimorder stimorder])

    set(a,'YTick',medcat)

    if t == 1 % do y-labels
        set(a,'YTickLabels',[stimorder stimorder])    
    else
        set(a,'YTickLabels',[])
    end

    annotation('textbox','Units','Pixels','Position', ...
        [310+(t-1)*450 410 50 20],...
        'LineStyle','none','String','Words',...
        'FontSize',20,'VerticalAlignment','middle','HorizontalAlignment','center');
    annotation('textbox','Units','Pixels','Position', ...
        [170+(t-1)*450 410 50 20],...
        'LineStyle','none','String','Images',...
        'FontSize',20,'VerticalAlignment','middle','HorizontalAlignment','center');

    % embedding
    fprintf('create embedding\n')
    Y = RDMmean;
    Y(eye(size(Y))==1)=0;
    X1 = mdscale(Y,2,'Start','random','Criterion','metricstress','Replicates',10);

    a = axes('Units','pixels','Position',[110+(t-1)*450 30 330 330]);
    a.FontSize = 18;
    a.LineWidth = 2;
    hold on
    X=X1;
    aw=range(X(:))*.04;
    axis equal
    axis square
    imrange = 30:220;
    plot(X(:,1),X(:,2),'.')
    for i=1:length(stims)
        IM = imread(sprintf('stimuli/stimuli/%s.png',stims{i}));
        IM = flipud(IM);
        if size(IM,3)==1
            IM = repmat(IM,1,1,3);
        end
        
        image('XData',X(i,1)+aw*[-1 1],'YData',X(i,2)+aw*[-1 1],'CData',IM(imrange,imrange,1:3))
    end
    xlabel('PC2')
    ylabel('PC1')
    a.XTick = [];
    a.YTick = [];

end

annotation('textbox','Units','Pixels','Position', ...
    [150 780 200 20],...
    'LineStyle','none','String','Image task',...
    'FontSize',20,'VerticalAlignment','middle','HorizontalAlignment','center');

annotation('textbox','Units','Pixels','Position', ...
    [600 780 200 20],...
    'LineStyle','none','String','Concept task',...
    'FontSize',20,'VerticalAlignment','middle','HorizontalAlignment','center');

%% annotations
text(-580,640,'Words','Rotation',90,'FontSize',20,'Units','pixels')
text(-580,490,'Images','Rotation',90,'FontSize',20,'Units','pixels')

annotation('textbox','Units','Pixels','Position', ...
    [10 1200 20 20],...
    'LineStyle','none','String','A',...
    'FontSize',40,'VerticalAlignment','middle','HorizontalAlignment','center');
annotation('textbox','Units','Pixels','Position', ...
    [480 1200 20 20],...
    'LineStyle','none','String','B',...
    'FontSize',40,'VerticalAlignment','middle','HorizontalAlignment','center');
annotation('textbox','Units','Pixels','Position', ...
    [10 780 20 20],...
    'LineStyle','none','String','C',...
    'FontSize',40,'VerticalAlignment','middle','HorizontalAlignment','center');
annotation('textbox','Units','Pixels','Position', ...
    [10 380 20 20],...
    'LineStyle','none','String','D',...
    'FontSize',40,'VerticalAlignment','middle','HorizontalAlignment','center');

set(gcf,'PaperOrientation','portrait');
fn = 'figures/Figure7_behaviour';
tn = tempname;
print(gcf,'-dpng','-r500',tn)
im=imread([tn '.png']);
[i,j]=find(mean(im,3)<255);margin=2;
imwrite(im(min(i-margin):max(i+margin),min(j-margin):max(j+margin),:),[fn '.png'],'png');
