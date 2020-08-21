function coactPaper_figS11()
labelWeight=true;
labelCase=true;
labelSize=9;
letGapX=3;
letGapY=6;

close all
fh=initFig('height',6);

% 
x=4;y=7;
panel_01(x,y);
panelLetter2(x-letGapX,y-letGapY+2,alphabet(1,labelCase),'fontSize',labelSize,'isBold',labelWeight)
drawnow();

x=4+40;y=7;
panel_02(x,y);
panelLetter2(x-letGapX,y-letGapY+2,alphabet(2,labelCase),'fontSize',labelSize,'isBold',labelWeight)
drawnow();

x=4+40*2;y=7;
panel_03(x,y);
panelLetter2(x-letGapX,y-letGapY+2,alphabet(3,labelCase),'fontSize',labelSize,'isBold',labelWeight)
drawnow();

print(fh,'figS11.pdf','-dpdf','-painters','-r300')

end

%%
function panel_01(x,y)

width=15;
totalHeigh=38;

gapY=1;
gapX=3.5;
smSigma=20; %in ms
cLim=0.01*[-1,1];
nShowBin=21;
ccg=poolVar('icaReacZNCCG.mat');
ccgSig=poolVar('icaReacCCG_sig.mat');

ccgAll=ccg;

col=setCoactColor;

ratList=fieldnames(ccg);

reg={};
sig=[];
ccgVal=[];
for ratIdx=1:length(ratList)
    rat=ratList{ratIdx};
    reg=[reg;ccg.(rat)(2).region(ccg.(rat)(2).pairID)];
    sig=[sig;ccgSig.(rat)(2).nrem.significance(:,[2,3])];
    ccgVal=cat(1,ccgVal,ccg.(rat)(2).nrem.real.ccg(:,:,2:3));
end

tBinSize=ccg.(ratList{1})(2).tBinSize*1e3;

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



ccgValAll=[];
for ratIdx=1:length(ratList)
    rat=ratList{ratIdx};
    ccgValAll=cat(1,ccgValAll,ccgAll.(rat)(2).nrem.real.ccg(:,:,2:3));
end

for n=1:2
    ccgValAll(:,:,n)=Filter0(smCore,ccgValAll(:,:,n));
end

ccgValAll=ccgValAll(:,cBin+(-nShowBin:nShowBin),:);
peakVal=mean(ccgValAll(:,nShowBin+1+(-3:3),2),2);


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
    subPeak=peakVal(idx);
    
    [~,order]=sort(subPeak,'descend');
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
                title({'Pre-conditioning'},'fontweight','normal','fontsize',5)
                textInMM(x+(width+gapX/2),y-4,'Entire NREM','fontsize',5,'horizontalAlign','center')
            else
                title({'Post-conditioning'},'fontweight','normal','fontsize',5)
            end
        end
        if n==3
            xlabel('\Deltatime (ms)','fontsize',5)
        end
        if m==0
            ylabel(join(target, ' - '),'fontsize',5)
        end
        if m==1 & n==2
            ps=get(gcf,'PaperSize')*10;
            xMM=x+(width+gapX)*m+width/2+[-1,0]-1.5;
            yMM=y+totalY+[2,1]+1;
            annotation('arrow',xMM/ps(1),1-yMM/ps(2),'color','w','HeadWidth',4,'HeadLength',4,...
            'LineWidth',1,'LineStyle','none')
        end
        if m==1 & n==1
            ps=get(gcf,'PaperSize')*10;
            xMM=x+(width+gapX)*m+width/2+[-2,0]-1.5;
            yMM=y+totalY+[4,2]+0.5;
            annotation('arrow',xMM/ps(1),1-yMM/ps(2),'color','w','HeadWidth',4,'HeadLength',4,...
            'LineWidth',1,'LineStyle','-')
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


subplotInMM(x,y+totalY+6,width*2+gapX,1)
imagescXY(cLim,[],linspace(cLim(1),cLim(2),512)');
set(gca,'clim',cLim)
colormap(gca,col.coact.map)
box off
set(gca,'YTick',[])
set(gca,'XTick',[cLim(1),cLim(1)/2,0,cLim(2)/2,cLim(2)])
xlabel('Correlation')
end

%%
function panel_02(x,y)

width=15;
totalHeigh=38;

gapY=1;
gapX=3.5;
smSigma=20; %in ms
cLim=0.01*[-1,1];
nShowBin=21;

ccg=poolVar('icaReacZNCCG_exSWR.mat');
ccgSig=poolVar('icaReacZNCCG_exSWR_sig.mat');
ccgAll=poolVar('icaReacZNCCG.mat');

col=setCoactColor;

ratList=fieldnames(ccg);

reg={};
sig=[];
ccgVal=[];
for ratIdx=1:length(ratList)
    rat=ratList{ratIdx};
    reg=[reg;ccg.(rat)(2).region(ccg.(rat)(2).pairID)];
    sig=[sig;ccgSig.(rat)(2).exRipple.significance(:,[2,3])];
    ccgVal=cat(1,ccgVal,ccg.(rat)(2).exRipple.real.ccg(:,:,2:3));
end

tBinSize=ccg.(ratList{1})(2).tBinSize*1e3;

nSm=ceil(smSigma/tBinSize);
xSm=(-nSm*4:nSm*4)*tBinSize;
smCore=normpdf(xSm,0,smSigma);
smCore=smCore/sum(smCore);

cBin=(size(ccgVal,2)+1)/2;

tBin=(-nShowBin:nShowBin)*tBinSize;

for n=1:2
    ccgVal(:,:,n)=Filter0(smCore,ccgVal(:,:,n));
end

ccgValAll=[];
for ratIdx=1:length(ratList)
    rat=ratList{ratIdx};
    ccgValAll=cat(1,ccgValAll,ccgAll.(rat)(2).nrem.real.ccg(:,:,2:3));
end

for n=1:2
    ccgValAll(:,:,n)=Filter0(smCore,ccgValAll(:,:,n));
end

ccgValAll=ccgValAll(:,cBin+(-nShowBin:nShowBin),:);
peakVal=mean(ccgValAll(:,nShowBin+1+(-3:3),2),2);

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
    subPeak=peakVal(idx);
    
    [~,order]=sort(subPeak,'descend');
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
                title({'Pre-conditioning'},'fontweight','normal','fontsize',5)
                textInMM(x+(width+gapX/2),y-4,'NREM excluding SWR','fontsize',5,'horizontalAlign','center')
            else
                title({'Post-conditioning'},'fontweight','normal','fontsize',5)
            end
        end
        if n==3
            xlabel('\Deltatime (ms)','fontsize',5)
        end
        if m==0
            ylabel(join(target, ' - '),'fontsize',5)
        end
        if m==1 & n==2
            ps=get(gcf,'PaperSize')*10;
            xMM=x+(width+gapX)*m+width/2+[-1,0]-1.5;
            yMM=y+totalY+[2,1]+1;
            annotation('arrow',xMM/ps(1),1-yMM/ps(2),'color','w','HeadWidth',4,'HeadLength',4,...
            'LineWidth',1,'LineStyle','none')
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


subplotInMM(x,y+totalY+6,width*2+gapX,1)
imagescXY(cLim,[],linspace(cLim(1),cLim(2),512)');
set(gca,'clim',cLim)
colormap(gca,col.coact.map)
box off
set(gca,'YTick',[])
set(gca,'XTick',[cLim(1),cLim(1)/2,0,cLim(2)/2,cLim(2)])
ax=fixAxis;
xlabel('Correlation')
end
%%
function panel_03(x,y)

width=15;
totalHeigh=38;

gapY=1;
gapX=3.5;
smSigma=20; %in ms
cLim=0.01*[-1,1];
nShowBin=21;

ccg=poolVar('icaReacZNCCG_exHFObaseCond.mat');
ccgSig=poolVar('icaReacZNCCG_exHFObaseCond_sig.mat');
ccgAll=poolVar('icaReacZNCCG.mat');

col=setCoactColor;

ratList=fieldnames(ccg);

reg={};
sig=[];
ccgVal=[];
for ratIdx=1:length(ratList)
    rat=ratList{ratIdx};
    reg=[reg;ccg.(rat)(2).region(ccg.(rat)(2).pairID)];
    sig=[sig;ccgSig.(rat)(2).exHFO.significance(:,[2,3])];
    ccgVal=cat(1,ccgVal,ccg.(rat)(2).exHFO.real.ccg(:,:,2:3));
end

tBinSize=ccg.(ratList{1})(2).tBinSize*1e3;

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


ccgValAll=[];
for ratIdx=1:length(ratList)
    rat=ratList{ratIdx};
    ccgValAll=cat(1,ccgValAll,ccgAll.(rat)(2).nrem.real.ccg(:,:,2:3));
end

for n=1:2
    ccgValAll(:,:,n)=Filter0(smCore,ccgValAll(:,:,n));
end

ccgValAll=ccgValAll(:,cBin+(-nShowBin:nShowBin),:);
peakVal=mean(ccgValAll(:,nShowBin+1+(-3:3),2),2);

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
    subPeak=peakVal(idx);
    
    [~,order]=sort(subPeak,'descend');
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
                title({'Pre-conditioning'},'fontweight','normal','fontsize',5)
                textInMM(x+(width+gapX/2),y-4,'NREM excluding HFO','fontsize',5,'horizontalAlign','center')
            else
                title({'Post-conditioning'},'fontweight','normal','fontsize',5)
            end
        end
        if n==3
            xlabel('\Deltatime (ms)','fontsize',5)
        end
        if m==0
            ylabel(join(target, ' - '),'fontsize',5)
        end
        if m==1 & n==1
            ps=get(gcf,'PaperSize')*10;
            xMM=x+(width+gapX)*m+width/2+[-2,0]-1.5;
            yMM=y+totalY+[4,2]+0.5;
            annotation('arrow',xMM/ps(1),1-yMM/ps(2),'color','w','HeadWidth',4,'HeadLength',4,...
            'LineWidth',1,'LineStyle','-')
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

subplotInMM(x,y+totalY+6,width*2+gapX,1)
imagescXY(cLim,[],linspace(cLim(1),cLim(2),512)');
set(gca,'clim',cLim)
colormap(gca,col.coact.map)
box off
set(gca,'YTick',[])
set(gca,'XTick',[cLim(1),cLim(1)/2,0,cLim(2)/2,cLim(2)])
xlabel('Correlation')
end
%%
