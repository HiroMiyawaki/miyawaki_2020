function coactPaper_figS08()
labelWeight=true;
labelCase=true;
labelSize=9;
letGapX=4;
letGapY=6;

close all
fh=initFig('height',4);

x=19;y=5;
panel_01(x,y);
drawnow();

print(fh,'figS08.pdf','-dpdf','-painters','-r300')
%
end

%%
function panel_01(x,y)
xGap=3;
width=(100-xGap)/2;
height=32;
doUpdate=false;

col=setCoactColor;
behList=fieldnames(col.state);
for n=1:3
    beh=behList{n};
    legSlp{n}=sprintf('\\color[rgb]{%f %f %f}%s',col.state.(beh),upper(beh))
end

if ~doUpdate &&...
        exist('~/Dropbox/FearAnalyses/png/example-coac-pre.png','file') && ...
        exist('~/Dropbox/FearAnalyses/png/example-coac-pre-info.mat','file')&& ...
        exist('~/Dropbox/FearAnalyses/png/example-coac-post.png','file') && ...
        exist('~/Dropbox/FearAnalyses/png/example-coac-post-info.mat','file')
else
    %%
    basename='~/data/Fear/triple/hoegaarden181115/hoegaarden181115';
    
    load([basename '.basicMetaData.mat'])
    load([basicMetaData.Basename '.sleepstate.states.mat'])
    load([basicMetaData.Basename '.sessions.events.mat'])
    load([basicMetaData.AnalysesName '-icaReac.mat'])
    load([basicMetaData.AnalysesName '-icaCoactTimeCondHT.mat'])
    load([basicMetaData.AnalysesName '-icaReacZNCCG_sig.mat'])
    
    %%
    slp=relabel_ma2sleep(SleepState.MECE.timestamps);
    slp(:,3)=(slp(:,3)+1)/2;
    slp(:,1:2)=slp(:,1:2)/60;
    slpCol=[1,0.5,0.5;
        0.5,0.5,1;
        0.7,0.7,1];
    %%
    tempSes=icaCoactTimeCond.param.templateIdx;
    tBin=(1:size(icaReac(tempSes).strength,2))*20e-3/60;
    
    withSlpState=true;

    exID=[6,8,14,25,28,33];
    regName=icaReac(tempSes).region(exID);
    regList=unique(regName);
    num=zeros(size(regList))
    for regID=1:length(regName)
        idx=find(strcmp(regName{regID},regList));
        num(idx)=num(idx)+1;
        enName{regID}=[regName{regID} '_{En' num2str(num(idx)) '}'];
    end
    
    dur=12;
    tRange=[287;471]+[0,dur];
    
    reac=zscore(icaReac(tempSes).strength(exID,:),[],2);
    
    yGapUit=ceil(diff(prctile(reac(:),[0.01,99.99]))/10)*10;    
    
    coact=[];
    cnt=0;
    for n=1:length(exID)-1
        for m=n+1:length(exID)
            if strcmp(regName{n},regName{m})
                continue
            end
            pId=find(any((icaCoactTimeCond.reacID==exID(n)),2)&any((icaCoactTimeCond.reacID==exID(m)),2));
            gap=icaCoactTimeCond.tGap(pId);
            if gap<0
                yVal=[reac(m,1-gap:end),zeros(1,-gap)];
            else
                yVal=[zeros(1,gap),reac(m,1:end-gap)];
            end
            cnt=cnt+1;
            coact(cnt,:)=reac(n,:).*yVal;
            pName{cnt}=[enName{n} ' - ' enName{m}];
        end
    end
    
    tTxt={'Pre-conditioning homecage session','Post-conditioning homecage session'};
    
    yGapUit=ceil(diff(prctile(coact(:),[0.01,99.99]))/10)*10;
    yGapStep=0:-1:-size(coact,1)+1
    
    [ytick,order]=sort(yGapStep);
    ytickLabel=pName(order);
    
    for prePost=1:2
        fhTemp=figure();
        set(fhTemp, 'paperUnit','centimeters','Units','centimeters')
        set(fhTemp,'position',[0,20,width/10,height/10])
        set(fhTemp, 'Units','centimeters')
        set(fhTemp,'PaperSize',[width/10,height/10])
        set(fhTemp,'paperPosition',[0,0,width/10,height/10])
        
        set(fhTemp,'defaultAxesFontName','Helvetica')
        set(fhTemp,'defaultTextFontName','Helvetica')
        
        set(fhTemp,'defaultAxesXColor',[0,0,0]); % factory is [0.15,0.15,0.15]
        set(fhTemp,'defaultAxesYColor',[0,0,0]);
        set(fhTemp,'defaultAxesZColor',[0,0,0]);
        
        set(fhTemp,'defaultAxesFontSize',5);
        set(fhTemp,'defaultTextFontSize',5);
        set(fhTemp,'defaultAxesLineWidth', 0.5);
        subplot('Position',[0,0,1,1]);
        
        subSlp=slp(slp(:,2)>tRange(prePost,1)&slp(:,1)<tRange(prePost,2),:);
        if subSlp(1,1)<tRange(prePost,1);subSlp(1,1)=tRange(prePost,1);end
        if subSlp(end,2)>tRange(prePost,2);subSlp(end,2)=tRange(prePost,2);end
        subSlp(:,1:2)=subSlp(:,1:2)-tRange(prePost,1);
        if withSlpState
            for sIdx=1:size(subSlp,1)
                rectangle('Position',[subSlp(sIdx,1),yGapUit*(min(yGapStep)-2),diff(subSlp(sIdx,1:2)),yGapUit*(-min(yGapStep)+4)],'LineStyle','none','FaceColor',slpCol(subSlp(sIdx,3),:))
            end
        end
        toShow=(tBin>=tRange(prePost,1)&tBin<=tRange(prePost,2));
        hold on
        plot(tBin(toShow)-tRange(prePost,1), coact(:,toShow)+yGapUit*yGapStep','k-','linewidth',0.25)
        ylim(yGapUit*[min(yGapStep)-2,2])
        
        axis off
        %         ax=fixAxis;
        xRange=get(gca,'XLim');
        yRange=get(gca,'YLim');
        tText=tTxt{prePost};
        if prePost==1
            fName='pre'
        else
            fName='post'
        end
        
        print(fhTemp,['~/Dropbox/FearAnalyses/png/example-coac-' fName '.png'],'-dpng','-r600')
        save(['~/Dropbox/FearAnalyses/png/example-coac-' fName '-info.mat'],...
            'xRange','yRange','tText','dur','ytickLabel','ytick','yGapUit')
        close(fhTemp)
    end
end
im{1}=imread('~/Dropbox/FearAnalyses/png/example-coac-pre.png');
info{1}=load('~/Dropbox/FearAnalyses/png/example-coac-pre-info.mat');
im{2}=imread('~/Dropbox/FearAnalyses/png/example-coac-post.png');
info{2}=load('~/Dropbox/FearAnalyses/png/example-coac-post-info.mat');

ytick=info{1}.ytick;
yGapUit=info{1}.yGapUit;
scalePos=1.1;
for prePost=1:2
    subplotInMM(x+(prePost-1)*(width+xGap),y,width,height)
    
    image(info{prePost}.xRange,info{prePost}.yRange,flipud(im{prePost}))
    set(gca,'YDir','normal')
    hold on
    plot(info{prePost}.dur*0.8+[0,1],yGapUit*(min(ytick)-scalePos)+[0,0],'k-','LineWidth',1)
    text(info{prePost}.dur*0.8+0.5,yGapUit*(min(ytick)-scalePos),'1 min','horizontalAlign','center','verticalAlign','top')
    plot(info{prePost}.dur*0.8+1.05+[0,0],yGapUit*(min(ytick)-scalePos)+[0,50],'k-','LineWidth',1)
    text(info{prePost}.dur*0.8+1.05,yGapUit*(min(ytick)-scalePos)+25,' 50 z^2','horizontalAlign','left','verticalAlign','middle')
    
    if prePost==1
        for n=1:length(info{prePost}.ytickLabel)
            text(info{prePost}.xRange(1)-diff(info{prePost}.xRange)*0.01,yGapUit*ytick(n),info{prePost}.ytickLabel{n},'horizontalAlign','right','fontsize',5)
        end
    end
    axis off
    title(info{prePost}.tText)
    ax=fixAxis;
        text2(0.98,-0.01,join(legSlp, ' '),ax,'verticalAlign','top','horizontalAlign','right','fontsize',5)
    
end

end
%%
