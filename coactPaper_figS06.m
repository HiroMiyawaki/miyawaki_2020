function coactPaper_figS06()
labelWeight=true;
labelCase=true;
labelSize=9;
letGapX=5;
letGapY=4;
%
close all
fh=initFig('height',6);

x=8;y=7;
panel_01(x,y);
panelLetter2(x-letGapX,y-letGapY,alphabet(1,labelCase),'fontSize',labelSize,'isBold',labelWeight)
drawnow();

x=8+48;y=7;
panel_02(x,y);
panelLetter2(x-letGapX+2,y-letGapY,alphabet(2,labelCase),'fontSize',labelSize,'isBold',labelWeight)
drawnow();

x=8+52;y=7+25;
panel_03(x,y);
panelLetter2(x-letGapX-2,y-letGapY-1,alphabet(3,labelCase),'fontSize',labelSize,'isBold',labelWeight)
drawnow();

print(fh,'figS06.pdf','-dpdf','-painters','-r300')
%
end

%%
function panel_01(xOri,yOri)

width=15;
totalHeigh=38;

tempIdx=2;
beh='rem';
tempName='conditioning';
x=xOri;
y=yOri+(totalHeigh+20)*(panelN-1);
gapY=1;
gapX=3.5;

smSigma=20; %in ms
cLim=0.01*[-1,1];
nShowBin=21;
ccg=poolVar('icaReacZNCCG.mat');
ccgSig=poolVar('icaReacCCG_sig.mat');

col=setCoactColor;

ratList=fieldnames(ccg);

reg={};
peakVal=[];
sig=[];
ccgVal=[];
for ratIdx=1:length(ratList)
    rat=ratList{ratIdx};
    reg=[reg;ccg.(rat)(tempIdx).region(ccg.(rat)(tempIdx).pairID)];
    peakVal=[peakVal;ccgSig.(rat)(tempIdx).(beh).peakValue(:,[2,3])];
    sig=[sig;ccgSig.(rat)(tempIdx).(beh).significance(:,[2,3])];
    ccgVal=cat(1,ccgVal,ccg.(rat)(tempIdx).(beh).real.ccg(:,:,2:3));
end

tBinSize=ccg.(ratList{1})(tempIdx).tBinSize*1e3;

nSm=ceil(smSigma/tBinSize);
xSm=(-nSm*4:nSm*4)*tBinSize;
smCore=normpdf(xSm,0,smSigma);
smCore=smCore/sum(smCore);

cBin=(size(ccgVal,2)+1)/2;

tBin=(-nShowBin:nShowBin)*tBinSize;

for n=1:2
    ccgVal(:,:,n)=Filter0(smCore,ccgVal(:,:,n));
end

ccgVal=ccgVal(:,cBin+(-nShowBin:nShowBin),:);
nPair=0;
for n=1:3
    switch n
        case 1
            target={'BLA','PrL L5'};
        case 2
            target={'vCA1','PrL L5'};
        case 3
            target={'vCA1','BLA'};
        otherwise
    end
    nPair=nPair+sum(strcmp(reg(:,1),target{1})&strcmp(reg(:,2),target{2}));
end

eachHight=(totalHeigh-gapY*2)/nPair;
%%
totalY=0;
for n=1:3
    switch n
        case 1
            target={'BLA','PrL L5'};
        case 2
            target={'vCA1','PrL L5'};
        case 3
            target={'vCA1','BLA'};
        otherwise
            continue
    end
    
    idx=find(strcmp(reg(:,1),target{1})&strcmp(reg(:,2),target{2}));
    subSig=sig(idx,:);
    subPeak=peakVal(idx,:);
    
    [~,order]=sort(mean(ccgVal(idx,nShowBin+1+(-3:3),2),2),'descend');
    idx=idx(order);
    subSig=subSig(order,:);
    
    height=length(idx)*eachHight;% cm
    for m=0:1
        subplotInMM(x+(width+gapX)*m,y+totalY,width,height,true)
        imagesc(tBin,1:length(idx),ccgVal(idx,:,1+m))
        box off
        set(gca,'ytick',[])
        if n~=3
            set(gca,'xtick',[])
        else
            set(gca,'xtick',-400:400:400)
        end
        xlim(400*[-1,1])
        set(gca,'clim',cLim)
        colormap(gca,col.coact.map)
        
        if n==1
            if m==0
                title({['Pre-' tempName], upper(beh)},'fontweight','normal','fontsize',5)
            else
                title({['Post-' tempName], upper(beh)},'fontweight','normal','fontsize',5)
            end
        end
        if n==3
            xlabel('\Deltatime (ms)','fontsize',5)
        end
        if m==0
            ylabel(join(target, ' - '),'fontsize',5)
        end
    end
    subplotInMM(x+width+0.5,y+totalY,2.5,height)
    imagesc([subSig(:,1),-2*ones(size(subSig(:,1))),subSig(:,2)])
    set(gca,'clim',[-2,1])
    colormap(gca,[1,1,1;flipud(col.pVal)])
    box off
    axis off
    
    totalY=totalY+height+gapY;
end

subplotInMM(x+width*2+gapX+0.5,y,1,totalY-gapY)
imagescXY([],cLim,linspace(cLim(1),cLim(2),512));
set(gca,'clim',cLim)
colormap(gca,col.coact.map)
box off
set(gca,'XTick',[])
set(gca,'YAxisLocation','right')
set(gca,'YTick',[cLim(1),cLim(1)/2,0,cLim(2)/2,cLim(2)])
ax=fixAxis;
text2(7,0.5,'Correlation',ax,'horizontalALign','center','Rotation',-90)
end
%%
function panel_02(xOri,yOri)
coact=poolVar('icaReacZNCCG_sig.mat');

width=9;
height=11;
withinGap=0;
acrossGap=4;

tempIdx=2;
beh='rem';

x=xOri;
y=yOri+(38+20)*(panelN-1);
sig=[];
reg={};
ratList=fieldnames(coact);
for ratIdx=1:length(ratList)
    rat=ratList{ratIdx};
    
    sig=[sig;coact.(rat)(tempIdx).(beh).significance(:,2:3)];
    
    tempReg=relabel_region(coact.(rat)(tempIdx).region,'minCellNum',0);
    reg=[reg;tempReg(coact.(rat)(tempIdx).pairID)];
    
end

[pairList,~,pairIdx]=uniqueCellRows(reg);

targetPair={'BLA','PrL L5'
    'vCA1','PrL L5'
    'vCA1','BLA'};
targetPairIdx=[]
for n=1:3
    targetPairIdx(n)=find(strcmp(pairList(:,1),targetPair{n,1})&strcmp(pairList(:,2),targetPair{n,2}));
end
temp=setCoactColor();
col=flipud(temp.pVal);
pFrac=[];
frac=[];

for n=1:length(targetPairIdx)
    subSig=sig(pairIdx==targetPairIdx(n),:);
    
    observed=[histcounts(subSig(:,1),-1.5:1.5);
        histcounts(subSig(:,2),-1.5:1.5)];
    frac(:,:,n)=observed;
    
    
    observed(:,sum(observed,1)==0)=[];
    
    if size(observed,2)<2
        pFrac(n)=1;
        continue
    end
    
    expect=sum(observed,2)*sum(observed,1)/sum(observed(:));
    
    chi2=sum((observed(:)-expect(:)).^2./expect(:));
    df=prod(size(observed)-1);
    pFrac(n)=chi2cdf(chi2,df,'upper');
end
prePost={'Pre','Post'};
for n=1:3
    subplotInMM(x+(2*width+withinGap+acrossGap)*(n-1),y,2*width+withinGap,height)
    
    sigFontSize=7;
    sigYshift=1.25;
    if pFrac(n)<0.001
        sigTxt='***';
    elseif pFrac(n)<0.01
        sigTxt='**';
    elseif pFrac(n)<0.05
        sigTxt='*';
    else
        sigTxt='';
        sigFontSize=5;
        sigYshift=0;
    end
    if ~isempty(sigTxt)
        plot(width/2+[0,0,1,1]*(width+withinGap),height-2+[0,1,1,0]*0.5,'k-','LineWidth',0.5)
        text(width+withinGap/2,height-1.5-sigYshift,sigTxt,'HorizontalAlignment','center','VerticalAlignment','bottom','FontSize',sigFontSize)
    end
    
    title(join(pairList(targetPairIdx(n),:),' - '),'fontweight','normal','fontsize',5)
    xlim([0,width*2+withinGap])
    ylim([0,height])
    axis off
    for m=1:2
        subplotInMM(x+(width+withinGap)*(m-1)+(2*width+withinGap+acrossGap)*(n-1),y+2,width,height-2,[],true)
        h=pie(frac(m,:,n),{'','',''});
        for hIdx=1:length(h)
            if strcmpi(h(hIdx).Type,'patch')
                h(hIdx).LineStyle='none';
            end
        end
        tempCol=col;
        tempCol(frac(m,:,n)<1,:)=[];
        colormap(gca,tempCol)
        ax=fixAxis;
        text2(0.5,0,prePost{m},ax,'verticalAlign','top','horizontalAlign','center')
    end
end

end

%%
function panel_03(xOri,yOri)
width=20;
height=12;

coact=poolVar('icaReacZNCCG_sig.mat');

tempIdx=2;
beh='rem';


x=xOri;
y=yOri+(38+20)*(panelN-1);

peak=[];
reg={};
ratList=fieldnames(coact);
for ratIdx=1:length(ratList)
    rat=ratList{ratIdx};
    
    peak=[peak;coact.(rat)(tempIdx).(beh).peakValue(:,2:3)];
    
    tempReg=relabel_region(coact.(rat)(tempIdx).region,'minCellNum',0);
    reg=[reg;tempReg(coact.(rat)(tempIdx).pairID)];
    
end

[pairList,~,pairIdx]=uniqueCellRows(reg);


targetPair={'BLA','PrL L5'
    'vCA1','PrL L5'
    'vCA1','BLA'};
targetPairIdx=[];
for n=1:3
    targetPairIdx(n)=find(strcmp(pairList(:,1),targetPair{n,1})&strcmp(pairList(:,2),targetPair{n,2}));
end
col=setCoactColor();

peakDiffMean=[];
peakDiffSTE=[];
pDiff=[];
for n=1:length(targetPairIdx)
    subPeak=peak(pairIdx==targetPairIdx(n),:);
    
    peakDiffMean(:,n)=mean(diff(subPeak,1,2));
    peakDiffSTE(:,n)=ste(diff(subPeak,1,2),[],1);
    if any(~isnan(diff(subPeak,1,2)))
        pDiff(n)=signrank(diff(subPeak,1,2));
    else
        pDiff(n)=1;
    end
    
end

subplotInMM(x,y,width,height)
hold on
for n=1:3
    
    temp=strrep(strrep(targetPair(n,:),' ',''),'/','');
    pairName=[temp{:}];
    bar(n,peakDiffMean(n),0.8,'LineStyle','none','FaceColor',col.pair.(pairName))
    posErr=peakDiffMean(n)+(2*(peakDiffMean(n)>0)-1)*peakDiffSTE(n);
    plot(n+[0,0],[peakDiffMean(n),posErr],'-','color',col.pair.(pairName))
    
    
    sigFontSize=7;
    sigYshift=0.2e-3;
    if pDiff(n)<0.001
        sigTxt='***';
    elseif pDiff(n)<0.01
        sigTxt='**';
    elseif pDiff(n)<0.05
        sigTxt='*';
    else
        sigTxt='';
        sigFontSize=5;
        sigYshift=1e-3;
    end
    if ~isempty(sigTxt)
        text(n,posErr+sigYshift,sigTxt,'FontSize',sigFontSize,'HorizontalAlignment','center')
    end
end
set(gca,'xtick',1:3,'XTickLabel',join(targetPair, ' - '),'XTickLabelRotation',30)
xlim([0,4])
ylim([-0.004,0.008])
set(gca,'YTick',-0.004:0.004:0.008)
ylabel('\Deltapeak correlation')
end
%%





