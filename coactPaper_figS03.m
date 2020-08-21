function coactPaper_figS03()
labelWeight=true;
labelCase=true;
labelSize=9;
letGapX=6;
letGapY=3;

close all
fh=initFig('height',10);


x=9;y=3;
panel_01(x,y);
panelLetter2(x-letGapX-2.5,y-letGapY+1,alphabet(1,labelCase),'fontSize',labelSize,'isBold',labelWeight)
drawnow();

x=9;y=30;
panel_02(x,y);
panelLetter2(x-letGapX-2.5,y-letGapY+1,alphabet(2,labelCase),'fontSize',labelSize,'isBold',labelWeight)
drawnow();


print(fh,'figS03.pdf','-dpdf','-painters','-r300')

end
%%

function panel_01(x,y)
load('~/data/Fear/triple/innis190601/innis190601.basicMetaData.mat')

load([basicMetaData.AnalysesName '-synCCG.mat']);
load([basicMetaData.AnalysesName '-synConn.mat']);
%%
nm=[65 56
    56 34];

yRange=[0,1200;
        0,500];
width=18;
height=15;
hGap=10;
wGap=9;
for idx=1:2
    subplotInMM(x+(width+wGap)*(idx-1),y,width,height)
    bar(synCCG.t,    synCCG.cnt(:,nm(idx,1),nm(idx,2)),1,'FaceColor',0*[1,1,1],'linestyle','none')
    hold on
    fill([synCCG.t,fliplr(synCCG.t)],[synCCG.confInt(:,nm(idx,1),nm(idx,2),1)',flipud(synCCG.confInt(:,nm(idx,1),nm(idx,2),2))'],...
        0.5*[1,1,1],'FaceAlpha',0.5,'linestyle','none')
    plot(synCCG.t([1,end]),    synCCG.globalBand(nm(idx,1),nm(idx,2),2)+[0,0],'c-','LineWidth',0.5)
    plot(synCCG.t([1,end]),    synCCG.globalBand(nm(idx,1),nm(idx,2),1)+[0,0],'r-','LineWidth',0.5)
    xlim(20*[-1,1])
    ylim(yRange(idx,:))
    ax=fixAxis;
    rectangle('position',[synConn.param.judgeWindow(1),ax(3),diff(synConn.param.judgeWindow),diff(ax(3:4))],...
                'linestyle','none','facecolor',[1,0.8,0])    
    h=get(gca,'Children');
    isRec=ismember(h,findobj(h,'type','rectangle'));
    set(gca,'Children',[h(~isRec);h(isRec)])
    
    
    ylabel('Counts')
    xlabel('\DeltaTime (ms)')
    box off
end

end
%%

function panel_02(x,y)
heightScatter=12;
heightWave=9;

wGap=6;
hGapInter=1;
hGapIntra=11;
nCol=5;
width=(110-wGap*(nCol-1))/nCol;

synConn=poolVar('synConn.mat');
unitInfo=poolVar('okUnit.cellinfo.mat');
ratList=fieldnames(synConn)
reg={};

for ratIdx=1:length(ratList)
    rat=ratList{ratIdx};
    reg=[reg,unitInfo.(rat).region];
end
[reg,regList]=relabel_region(reg);
reg(ismember(reg,regList(9:end)))=regList(end);
regList(9:end-1)=[];
regList{end}='Other regions';
reg(strcmp(reg,'other'))=regList(end);
%%
rise=[];
decay=[];
halfwidth=[];
troughAmp=[];
peakTroughAmp=[];
isPositive=[];
fr=[];
nzFr=[];
for ratIdx=1:length(ratList)
    rat=ratList{ratIdx};
    
    rise=[rise,unitInfo.(rat).waveform.rise'];
    decay=[decay,unitInfo.(rat).waveform.decay'];
    halfwidth=[halfwidth,unitInfo.(rat).waveform.halfwidth'];
    troughAmp=[troughAmp,unitInfo.(rat).waveform.troughAmp'];
    peakTroughAmp=[peakTroughAmp,unitInfo.(rat).waveform.peakTroughAmp'];
    isPositive=[isPositive,unitInfo.(rat).waveform.positiveSpike'];
    
    fr=[fr,unitInfo.(rat).FR.mean];
    nzFr=[nzFr,unitInfo.(rat).FR.nonZeroMean];
end
cType=[];
for ratIdx=1:length(ratList)
    rat=ratList{ratIdx};
    
    nEx=sum(synConn.(rat).inspected.map==1,2);
    nInh=sum(synConn.(rat).inspected.map==-1,2);
    
    cType=[cType,((nEx>0 & nInh==0)-(nEx==0 & nInh>0))'+2];
end
% %%
waveForm=[];
for ratIdx=1:length(ratList)
    rat=ratList{ratIdx};
    
    for n=1:size(unitInfo.(rat).waveform.wave,2)
        [~,ch]=max(range(unitInfo.(rat).waveform.wave(n).mean,2));
        waveForm(end+1,:)=unitInfo.(rat).waveform.wave(n).mean(ch,:);
    end
end
%%
col=setCoactColor;
colCellType=[col.cellType.inh;col.cellType.nc;col.cellType.ex];

for regIdx=1:9
    subplotInMM(x+(width+wGap)*mod((regIdx-1),nCol),y+(heightWave+heightScatter+hGapInter+hGapIntra)*(ceil(regIdx/nCol)-1),width,heightWave)
    hold on
    t=(-16:24)/20;
    for n=[2,3,1]
        target=find(strcmp(reg,regList{regIdx})& ~isPositive & cType==n)
        if isempty(target)
            continue
        end
        wav=waveForm(target,:)';
        plot(t,wav./(max(-wav,[],1)),'-','color',colCellType(n,:))
    end
    ylim([-1,1])
    xlim([-0.8,1.2])
    title(regList{regIdx},'fontsize',5,'fontweight','normal')
    plot(0.6+[0,0.2],-0.6+[0,0],'k-')
    text(0.6+0.1,-0.6,'0.2 ms', 'VerticalAlignment','top','HorizontalAlignment','center','fontsize',5)
    axis off
    
    ax=fixAxis;
    if regIdx==9
        legTxt{1}=sprintf('\\color[rgb]{%f %f %f}Classified as excitatory',col.cellType.ex);
        legTxt{2}='based on CCG';
        legTxt{3}=sprintf('\\color[rgb]{%f %f %f}Classified as inhibitory',col.cellType.inh);
        legTxt{4}='based on CCG';
        legTxt{5}=sprintf('\\color[rgb]{%f %f %f}Not classified',col.cellType.nc);
        legTxt{6}='based on CCG';
        text2(1.1,0.9,legTxt,ax,'verticalALign','top','fontsize',5)
    end
    
    subplotInMM(x+(width+wGap)*mod((regIdx-1),nCol),y+(heightWave+heightScatter+hGapInter+hGapIntra)*(ceil(regIdx/nCol)-1)+heightWave+hGapInter,width,heightScatter)
    threshold=[0.5,0.6];
    
            X=decay;
            xTxt='Spike width (ms)';
            xbin=0:0.05:1.5;
            xRange=[0,1.0];

            Y=log10(fr);
            yTxt='Mean FR (Hz)';
            yLog=true;
            yTick=-2:2;
            yRange=[-2,2];
            ybin=-3:0.1:3;
            
    target=strcmp(reg,regList{regIdx})& ~isPositive;
hold on
        scatter(X(target),Y(target),4,colCellType(cType(target),:),'fill','MarkerFaceAlpha',0.5)
        xlabel(xTxt)
        xlim(xRange)
        if mod(regIdx,nCol)==1
            ylabel('FR (Hz)')
        end
        ylim(yRange)
        set(gca,'ytick',yTick,'yTickLabel',10.^yTick)        
        xlim(xRange)
        box off


    ax=fixAxis;
    plot(threshold(1)+[0,0],ax(3:4),'-','color',0.5*[1,1,1])
    plot(threshold(2)+[0,0],ax(3:4),'-','color',0.5*[1,1,1])
end
end



