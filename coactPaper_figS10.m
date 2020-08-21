function coactPaper_figS10()
labelWeight=true;
labelCase=true;
labelSize=9;
letGapX=5;
letGapY=6;

close all
fh=initFig('height',17);

x=8;y=3;
panel_01(x,y);
panelLetter2(x-letGapX-2,y-letGapY+4,alphabet(1,labelCase),'fontSize',labelSize,'isBold',labelWeight)
drawnow();

x=69;y=5;
panel_02(x,y);
panelLetter2(x-letGapX-1,y-letGapY+2,alphabet(2,labelCase),'fontSize',labelSize,'isBold',labelWeight)
drawnow();

y=5+57;
panel_03(x,y);
panelLetter2(x-letGapX-1,y-letGapY+2,alphabet(3,labelCase),'fontSize',labelSize,'isBold',labelWeight)
drawnow();

y=5+57*2;
panel_04(x,y);
panelLetter2(x-letGapX-1,y-letGapY+2,alphabet(4,labelCase),'fontSize',labelSize,'isBold',labelWeight)
drawnow();

print(fh,'figS10.pdf','-dpdf','-painters','-r300')

end

%%
function panel_01(x,y)
wGapInter=13;
wGapIntra=3;
wForScale=3;
width=((170-wGapInter*2-wForScale*3)/3-wGapIntra*2)/3;

height=13;
hGap=5;
%%
hfoTrigHist= poolVar('hfoTrigHist.mat');
swrTrigHist= poolVar('swrTrigHist.mat');
spdlTrigHist= poolVar('spindleTrigHist.mat');
CellInfo= poolVar('okUnit.cellinfo.mat');


ratList=fieldnames(hfoTrigHist);
%%
reg={};
peth.hfo=[];
peth.swr=[];
peth.spdl=[];

cellType=[];
for ratIdx=1:length(ratList)
    ratName=ratList{ratIdx};
        
    reg=[reg,CellInfo.(ratName).region];
    
    nTrig=0;
    sumSpk=zeros(size(hfoTrigHist.(ratName).smZ.nrem,1),size(hfoTrigHist.(ratName).smZ.nrem,2));
    for n=1:size(hfoTrigHist.(ratName).smZ.nrem,3)    
        sumSpk=sumSpk+hfoTrigHist.(ratName).smZ.nrem(:,:,n)*hfoTrigHist.(ratName).triger.nrem.n(n);
        nTrig=nTrig+hfoTrigHist.(ratName).triger.nrem.n(n);
    end    
    peth.hfo=cat(1,peth.hfo,sumSpk/nTrig);
    
    if isfield(swrTrigHist,ratName)
        nTrig=0;
        sumSpk=zeros(size(swrTrigHist.(ratName).smZ.nrem,1),size(swrTrigHist.(ratName).smZ.nrem,2));
        for n=1:size(swrTrigHist.(ratName).smZ.nrem,3)    
            sumSpk=sumSpk+swrTrigHist.(ratName).smZ.nrem(:,:,n)*swrTrigHist.(ratName).triger.nrem.n(n);
            nTrig=nTrig+swrTrigHist.(ratName).triger.nrem.n(n);
        end    
        peth.swr=cat(1,peth.swr,sumSpk/nTrig);    
    else
        sumSpk=nan(size(hfoTrigHist.(ratName).smZ.nrem,1),size(hfoTrigHist.(ratName).smZ.nrem,2));
        peth.swr=cat(1,peth.swr,sumSpk/nTrig);            
    end

    nTrig=0;
    sumSpk=zeros(size(spdlTrigHist.(ratName).pfc.smZ,1),size(spdlTrigHist.(ratName).pfc.smZ,2));
    for n=1:size(spdlTrigHist.(ratName).pfc.smZ,3)    
        sumSpk=sumSpk+spdlTrigHist.(ratName).pfc.smZ(:,:,n)*spdlTrigHist.(ratName).pfc.triger.n(n);
        nTrig=nTrig+spdlTrigHist.(ratName).pfc.triger.n(n);
    end
    
    peth.spdl=cat(1,peth.spdl,sumSpk/nTrig);    
    
    cellType=[cellType,CellInfo.(ratName).cellType.type];
    
end
[reg,regList]=relabel_region(reg);
reg(ismember(reg,regList(9:end)))=regList(end);
regList(9:end-1)=[];
regList{end}='Other regions';
reg(strcmp(reg,'other'))=regList(end);
%%
trigList={'hfo','swr','spdl'};
evtName.hfo='HFO';
evtName.swr='SWR';
evtName.spdl='spindle';

tWin.spdl=spdlTrigHist.(ratList{1}).param.nHalfwin*[-1,1]*spdlTrigHist.(ratList{1}).param.tBinSize;
tWin.swr=swrTrigHist.(ratList{1}).param.nHalfwin*[-1,1]*swrTrigHist.(ratList{1}).param.tBinSize*1e3;
tWin.hfo=hfoTrigHist.(ratList{1}).param.nHalfwin*[-1,1]*hfoTrigHist.(ratList{1}).param.tBinSize*1e3;

col=setCoactColor();
for n=1:9
    target=find(strcmp(reg,regList{n})&abs(cellType)==1);
    
    peak=mean(peth.hfo(target,501+(-1:1)),2);
    cType=cellType(target)';
    [~,order]=sortrows([cType,peak],'descend');
    
    posY=n-1;
    posX=0;
    
    temp=[peth.hfo(target,:),peth.swr(target,:) peth.spdl(target,:)];
    cLim=prctile(temp(:),[1,99]);
    cLim=[floor(cLim(1)*10)/10,ceil(cLim(2)*10)/10];
    
    nEx=sum(cType==1);

    
    for m=1:3
        subset=peth.(trigList{m})(target,:);
        subset=subset(order,:);
        
        subplotInMM(x+(3*width+2*wGapIntra+wGapInter+wForScale)*posX+(width+wGapIntra)*(m-1),...
                    y+(height+hGap)*posY,width,height)
        imagesc(tWin.(trigList{m}),[],subset)
        box off
        if m==1
            ylabel(regList{n})
        else
            set(gca,'YTickLabel',[])
        end
            if m==3
                xlim(2*[-1,1])
            else
                xlim(420*[-1,1])
            end
        if mod(n,9)==0
            if m==3
                xlabel({'Time from' [evtName.(trigList{m}) ' peak (s)']})
            else
                xlabel({'Time from' [evtName.(trigList{m}) ' peak (ms)']})
            end
        end
        hold on
        plot(tWin.(trigList{m}),nEx+[0,0],'c-','linewidth',0.5)
        set(gca,'cLim',cLim)
        colormap(gca,col.fr)
    end
    
    subplotInMM(x+(3*width+2*wGapIntra+wGapInter+wForScale)*posX+width*3+wGapIntra*2+1,...
                y+(height+hGap)*posY,1,height);
    imagescXY([0,1],cLim,linspace(cLim(1),cLim(2),size(col.fr,1)))
    colormap(gca,col.fr)
    set(gca,'CLim',cLim,'XTick',[],'YAxisLocation','right')
    box off
    xlim([0,1])
    ax=fixAxis;
    if n<4
        text(5.5,mean(cLim),'Firing rate (z)','horizontalAlign','center','rotation',-90,'fontsize',5)        
    elseif n>3 && n<7
        text(5.5,mean(cLim),'Firing rate (z)','horizontalAlign','center','rotation',-90,'fontsize',5)
    else
        text(5.5,mean(cLim),'Firing rate (z)','horizontalAlign','center','rotation',-90,'fontsize',5)
    end
    
end

end
%%
function panel_02(x,y);
xGap=2;
width=(45-xGap)/2-1;
yGap=1;
totalHeight=40;
eachPETH=poolVar('hfoTrigIcaReac.mat');

ratList=fieldnames(eachPETH);

col=setCoactColor;
%%
smSigma=20/1000; %in s
cLim=4*[-1,1];


sesIdx=2;
    tempName=eachPETH.(ratList{1})(sesIdx).tempName;
    tBinSize=eachPETH.(ratList{1})(sesIdx).tBinSize;
    tWin=eachPETH.(ratList{1})(sesIdx).param.tWindow;
    nWin=ceil(tWin/tBinSize);
    tBin=(-nWin:nWin)*tBinSize*1e3;
    peth=[];
    reg={};
    for ratIdx=1:length(ratList)
        rat=ratList{ratIdx};
        peth=cat(1,peth,eachPETH.(rat)(sesIdx).nrem.mean);
        reg=[reg,eachPETH.(rat)(sesIdx).region];
    end
    %%
    [reg,regList]=relabel_region(reg);
    regID=zeros(size(reg));
    for regIdx=1:length(regList)
        regID(strcmp(reg,regList{regIdx}))=regIdx;
    end
    
    nPair=sum(ismember(reg,{'vCA1','BLA','PrL L5'}));
    eachHeight=(totalHeight-yGap*2)/nPair;

    
    xRange=[-420,420];
    
    periodText={'Pre','Post'};
    xTxt={'Time form' 'HFO peak (ms)'};
    
    nSm=ceil(smSigma/tBinSize);
    xSm=(-nSm*4:nSm*4)*tBinSize;
    smCore=normpdf(xSm,0,smSigma);
    smCore=smCore/sum(smCore);
   

        totalY=0;
        for n=1:3
            switch n
                case 1
                    target='vCA1';
                case 2
                    target='BLA';
                case 3
                    target='PrL L5';
                otherwise
                    continue
            end
            idx=find(strcmp(reg,target));
        
            temp=zscore([peth(idx,:,1)';peth(idx,:,2)']);
        subsetZ{1}=Filter0(smCore,temp(1:nWin*2+1,:));
        subsetZ{2}=Filter0(smCore,temp(nWin*2+1+1:end,:));
        
        peak=mean(subsetZ{2}(nWin+1+(-1:1),:),1);
        [~,order]=sort(peak);
        
            height=length(idx)*eachHeight;
        for prePost=1:2
            
            subplotInMM(x+(width+xGap)*(prePost-1),y+totalY,width,height)

            imagescXY(tBin,[],subsetZ{prePost}(:,order));
            box off
            xlim(xRange)
            set(gca,'YTick',[])
            set(gca,'clim',cLim)
            ax=fixAxis;
            if n==3
                xlabel(xTxt)
            else
                set(gca,'XTickLabel',[])
            end
            if n==1
                text2(0.5,1,{sprintf('%s-conditioning',periodText{prePost}),'NREM'},ax,...
                    'fontsize',5,'horizontalAlign','center','verticalAlign','bottom')
            end
            if prePost==1
                ylabel(sprintf('%s ensemble',target),'fontsize',5,'fontweight','normal')
            end
            colormap(gca,col.react.map)            
        end        
        totalY=totalY+height+yGap;
    end
    subplotInMM(x+width*2+xGap+1,y,1,totalY-yGap)
    imagescXY([0,1],cLim,linspace(cLim(1),cLim(2),256));
    box off
    set(gca,'XTick',[],'YTick',[cLim(1),cLim(1)/2,0,cLim(2)/2,cLim(2)],'YAxisLocation','right')
    colormap(gca,col.react.map)   
    ax=fixAxis;
    text2(5,0.5,'Ensemble activation strength (z)',ax,'fontsize',5,'horizontalALign','center','rotation',-90)
end

function panel_03(x,y);
xGap=2;
width=(45-xGap)/2-1;
yGap=1.5;
totalHeight=40;
eachPETH=poolVar('swrTrigIcaReac.mat');

ratList=fieldnames(eachPETH);

col=setCoactColor;
%%
smSigma=20/1000; %in s
cLim=4*[-1,1];


sesIdx=2;
    tempName=eachPETH.(ratList{1})(sesIdx).tempName;
    tBinSize=eachPETH.(ratList{1})(sesIdx).tBinSize;
    tWin=eachPETH.(ratList{1})(sesIdx).param.tWindow;
    nWin=ceil(tWin/tBinSize);
    tBin=(-nWin:nWin)*tBinSize*1e3;
    peth=[];
    reg={};
    for ratIdx=1:length(ratList)
        rat=ratList{ratIdx};
        peth=cat(1,peth,eachPETH.(rat)(sesIdx).nrem.mean);
        reg=[reg,eachPETH.(rat)(sesIdx).region];
    end
    %%
    [reg,regList]=relabel_region(reg);
    regID=zeros(size(reg));
    for regIdx=1:length(regList)
        regID(strcmp(reg,regList{regIdx}))=regIdx;
    end
    
    nPair=sum(ismember(reg,{'vCA1','BLA','PrL L5'}));
    eachHeight=(totalHeight-yGap*2)/nPair;

    
    xRange=[-420,420];
    
    periodText={'Pre','Post'};
    xTxt={'Time form' 'SWR peak (ms)'};
    
    nSm=ceil(smSigma/tBinSize);
    xSm=(-nSm*4:nSm*4)*tBinSize;
    smCore=normpdf(xSm,0,smSigma);
    smCore=smCore/sum(smCore);
   

        totalY=0;
        for n=1:3
            switch n
                case 1
                    target='vCA1';
                case 2
                    target='BLA';
                case 3
                    target='PrL L5';
                otherwise
                    continue
            end
            idx=find(strcmp(reg,target));
        
            temp=zscore([peth(idx,:,1)';peth(idx,:,2)']);
        subsetZ{1}=Filter0(smCore,temp(1:nWin*2+1,:));
        subsetZ{2}=Filter0(smCore,temp(nWin*2+1+1:end,:));
        
        peak=mean(subsetZ{2}(nWin+1+(-1:1),:),1);
        [~,order]=sort(peak);
        
            height=length(idx)*eachHeight;
        for prePost=1:2
            
            subplotInMM(x+(width+xGap)*(prePost-1),y+totalY,width,height)

            imagescXY(tBin,[],subsetZ{prePost}(:,order));
            box off
            xlim(xRange)
            set(gca,'YTick',[])
            set(gca,'clim',cLim)
            ax=fixAxis;
            if n==3
                xlabel(xTxt)
            else
                set(gca,'XTickLabel',[])
            end
            if n==1
                text2(0.5,1,{sprintf('%s-conditioning',periodText{prePost}),'NREM'},ax,...
                    'fontsize',5,'horizontalAlign','center','verticalAlign','bottom')
            end
            if prePost==1
                ylabel(sprintf('%s ensemble',target),'fontsize',5,'fontweight','normal')
            end
            colormap(gca,col.react.map)            
        end        
        totalY=totalY+height+yGap;
    end
    subplotInMM(x+width*2+xGap+1,y,1,totalY-yGap)
    imagescXY([0,1],cLim,linspace(cLim(1),cLim(2),256));
    box off
    set(gca,'XTick',[],'YTick',[cLim(1),cLim(1)/2,0,cLim(2)/2,cLim(2)],'YAxisLocation','right')
    colormap(gca,col.react.map)   
    ax=fixAxis;
    text2(5,0.5,'Ensemble activation strength (z)',ax,'fontsize',5,'horizontalALign','center','rotation',-90)
end

%%
function panel_04(x,y);
xGap=2;
width=(45-xGap)/2-1;
yGap=1.5;
eachHeight=0.32;
totalHeight=40;
eachPETH=poolVar('spdlTrigIcaReac.mat');

ratList=fieldnames(eachPETH);

col=setCoactColor;
%%
smSigma=100/1000; %in s
cLim=4*[-1,1];


sesIdx=2;
    tempName=eachPETH.(ratList{1})(sesIdx).tempName;
    tBinSize=eachPETH.(ratList{1})(sesIdx).tBinSize;
    tWin=eachPETH.(ratList{1})(sesIdx).param.tWindow;
    nWin=ceil(tWin/tBinSize);
    tBin=(-nWin:nWin)*tBinSize;
    peth=[];
    reg={};
    for ratIdx=1:length(ratList)
        rat=ratList{ratIdx};
        peth=cat(1,peth,eachPETH.(rat)(sesIdx).pfc.mean);
        reg=[reg,eachPETH.(rat)(sesIdx).region];
    end
    %%
    [reg,regList]=relabel_region(reg);

        
    nPair=sum(ismember(reg,{'vCA1','BLA','PrL L5'}));
    eachHeight=(totalHeight-yGap*2)/nPair;
    
    regID=zeros(size(reg));
    for regIdx=1:length(regList)
        regID(strcmp(reg,regList{regIdx}))=regIdx;
    end
    xRange=[-2,2];
    
    periodText={'Pre','Post'};
    xTxt={'Time form' 'Spindle peak (s)'};
    
    nSm=ceil(smSigma/tBinSize);
    xSm=(-nSm*4:nSm*4)*tBinSize;
    smCore=normpdf(xSm,0,smSigma);
    smCore=smCore/sum(smCore);
   

        totalY=0;
        for n=1:3
            switch n
                case 1
                    target='vCA1';
                case 2
                    target='BLA';
                case 3
                    target='PrL L5';
                otherwise
                    continue
            end
            idx=find(strcmp(reg,target));
        
            temp=zscore([peth(idx,:,1)';peth(idx,:,2)']);
        subsetZ{1}=Filter0(smCore,temp(1:nWin*2+1,:));
        subsetZ{2}=Filter0(smCore,temp(nWin*2+1+1:end,:));
        
        peak=mean(subsetZ{2}(nWin+1+(-1:1),:),1);
        [~,order]=sort(peak);
        
            height=length(idx)*eachHeight;
        for prePost=1:2
            
            subplotInMM(x+(width+xGap)*(prePost-1),y+totalY,width,height)

            imagescXY(tBin,[],subsetZ{prePost}(:,order));
            box off
            xlim(xRange)
            set(gca,'YTick',[])
            set(gca,'clim',cLim)
            ax=fixAxis;
            if n==3
                xlabel(xTxt)
            else
                set(gca,'XTickLabel',[])
            end
            if n==1
                text2(0.5,1,{sprintf('%s-conditioning',periodText{prePost}),'NREM'},ax,...
                    'fontsize',5,'horizontalAlign','center','verticalAlign','bottom')
            end
            if prePost==1
                ylabel(sprintf('%s ensemble',target),'fontsize',5,'fontweight','normal')
            end
            colormap(gca,col.react.map)            
        end        
        totalY=totalY+height+yGap;
    end
    subplotInMM(x+width*2+xGap+1,y,1,totalY-yGap)
    imagescXY([0,1],cLim,linspace(cLim(1),cLim(2),256));
    box off
    set(gca,'XTick',[],'YTick',[cLim(1),cLim(1)/2,0,cLim(2)/2,cLim(2)],'YAxisLocation','right')
    colormap(gca,col.react.map)   
    ax=fixAxis;
    text2(5,0.5,'Ensemble activation strength (z)',ax,'fontsize',5,'horizontalALign','center','rotation',-90)
end

