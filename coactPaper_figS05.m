function coactPaper_figS05()
labelWeight=true;
labelCase=true;
labelSize=9;
letGapX=5;
letGapY=4;
%
close all
fh=initFig('height',6);
%
x=15;y=5;
panel_01(x,y);
drawnow();

print(fh,'figS05.pdf','-dpdf','-painters','-r300')
%
end

%%
function panel_01(x,y)

width=15;
totalHeigh=38*2+20;


tempIdx=2;
tempName='conditioning';
beh='nrem';
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

%%
[regList,~,regIdx]=uniqueCellRows(reg);
nCol=0;
totalY=0;
pairList=[];
for n=1:size(regList,1)
    if strcmp(regList{n,1},regList{n,2})
        continue
    end
    if all(strcmp(regList(n,:),{'BLA','PrL L5'}))
        continue
    end
    if all(strcmp(regList(n,:),{'vCA1','PrL L5'}))
        continue
    end
    if all(strcmp(regList(n,:),{'vCA1','BLA'}))
        continue
    end
    pairList(end+1)=n;
end




eachHight=(totalHeigh-gapY*(length(pairList)-1))/sum(ismember(regIdx,pairList));


for nn=1:length(pairList)
    n=pairList(nn);
    idx=find(regIdx==n);
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
        if nn~=length(pairList) && nn~=9
            set(gca,'xtick',[])
        else
            set(gca,'xtick',-400:400:400)
        end
        xlim(400*[-1,1])
        set(gca,'clim',cLim)
        colormap(gca,col.coact.map)
        
        if totalY==0
            if m==0
                title({['Pre-' tempName], upper(beh)},'fontweight','normal','fontsize',5)
            else
                title({['Post-' tempName], upper(beh)},'fontweight','normal','fontsize',5)
            end
        end
        if nn==length(pairList) || nn==9
            xlabel('\Deltatime (ms)','fontsize',5)
        end
        if m==0
            ax=fixAxis;
            text2(-0.05,0.5,join(regList(n,:) , ' - '),ax,'fontsize',5,'horizontalALign','right')
        end
    end

        subplotInMM(x+width+0.5,y+totalY,2.5,height)
        imagesc([subSig(:,1),-2*ones(size(subSig(:,1))),subSig(:,2)])
        set(gca,'clim',[-2,1])
        colormap(gca,[1,1,1;flipud(col.pVal)])
        box off
        axis off
        
    totalY=totalY+height+gapY;
    
    if totalY>48 &&nCol==0
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
        totalY=0;
        x=x+width*2+gapX+0.5+25;
        nCol=1;
    end
    
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






