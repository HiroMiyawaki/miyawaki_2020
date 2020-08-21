function coactPaper_science_figS13()
labelWeight=true;
labelCase=true;
labelSize=9;
letGapX=5;
letGapY=6;

close all
fh=initFig('height',12);

x=13;y=5;
panel_01(x,y);
panelLetter2(x-letGapX-4,y-letGapY+2,alphabet(1,labelCase),'fontSize',labelSize,'isBold',labelWeight)
drawnow();

x=13;y=5+33+2;
panel_02(x,y);
panelLetter2(x-letGapX-4,y-letGapY+1,alphabet(2,labelCase),'fontSize',labelSize,'isBold',labelWeight)

x=13+70;y=5;
panel_03(x,y);
panelLetter2(x-letGapX-4,y-letGapY+2,alphabet(3,labelCase),'fontSize',labelSize,'isBold',labelWeight)
drawnow();

print(fh,'figS13.pdf','-dpdf','-painters','-r300')
%
end
%%

function panel_01(x,y)

width=18-8/3;
height=8;
xGapIntra=4.5;
yGapIntra=4;

evtTrigReact=poolVar('shockTrigIcaReact.mat');

ratList=fieldnames(evtTrigReact);
%%

t=evtTrigReact.(ratList{1}).time;
evtRate=[];
evtPeak=[];
evtStr=[];
reg={};
withPartner=[];
for ratIdx=1:length(ratList)
    rat=ratList{ratIdx};
    evtRate=cat(1,evtRate,evtTrigReact.(rat).rate);
    evtPeak=cat(1,evtPeak,evtTrigReact.(rat).peak);
    evtStr=cat(1,evtStr,evtTrigReact.(rat).strength);
    
    reg=cat(1,reg,evtTrigReact.(rat).region');
    withPartner=cat(1,withPartner,~cellfun(@isempty,evtTrigReact.(rat).partner)');
end
%%

[regList,~,regListIdx]=unique(reg);
%%
tBinSize=mean(diff(t));
smSigma=40/1000;
smBin=(0:ceil(smSigma*4/tBinSize))*tBinSize;
smBin=[-fliplr(smBin),smBin(2:end)]
smCore=normpdf(smBin,0,smSigma);
smCore=smCore/sum(smCore);

targetReg={'vCA1','BLA','PrL L5'};

for yType=1:2
    avg={};
    err={};
    if yType==1
        val=evtRate;
        yTxt={'Occurrence' 'rate (1/s)'};
        yLim=[0,6];
    else
        val=evtPeak;
        yTxt='Peak (z)';
        yLim=[0,20];
    end
    
    for regIdx=1:size(regList,1);
        targetBool{1}=find(regListIdx==regIdx&withPartner);
        targetBool{2}=find(regListIdx==regIdx&~withPartner);
        for pType=1:2
            target=targetBool{pType};
            peth=val(target,:);
            peth(isnan(peth))=0;
            
            
            if smSigma>0
                for n=1:size(peth,1)
                    peth(n,:)=Filter0(smCore,peth(n,:));
                end
            end
            avg{regIdx,pType}=nanmean(peth,1);
            err{regIdx,pType}=nanste(peth,[],1);
        end
        
        p=ones(1,size(val,2));
        for tIdx=1:size(val,2)
            if sum(~isnan(val(targetBool{1},tIdx)))>0 && sum(~isnan(val(targetBool{2},tIdx)))>0
                p(tIdx)=ranksum(val(targetBool{1},tIdx),val(targetBool{2},tIdx));
            end
        end
        pSig=p<0.01;
        
        sigOnset{regIdx}=find(diff(pSig)==1)+1;
        if pSig(1); sigOnset{regIdx}=[1,sigOnset{regIdx}];end
        
        sigOffset{regIdx}=find(diff(pSig)==-1);
        if pSig(end); sigOffset{regIdx}=[sigOffset{regIdx},length(pSig)];end
    end
    
    colDef=setCoactColor();
    
    col=[colDef.region.vCA1
        0.3*[1,1,1]
        colDef.region.BLA
        0.3*[1,1,1]
        colDef.region.PrLL5
        0.3*[1,1,1]
        ];
    
    for pairIdx=1:size(targetReg,2)
        regIdx=find(strcmp(regList,targetReg{pairIdx}));
        
        subplotInMM(x+(width+xGapIntra)*(pairIdx-1),...
            y+(height+yGapIntra)*(yType-1),width,height)
        
        hold on
        fill([t,fliplr(t)],[avg{regIdx,2}+err{regIdx,2},...
            fliplr(avg{regIdx,2}-err{regIdx,2})],...
            col(pairIdx*2,:),'linestyle','none','FaceAlpha',0.5)
        plot(t,avg{regIdx,2},'-','color',col(pairIdx*2,:))
        
        fill([t,fliplr(t)],[avg{regIdx,1}+err{regIdx,1},...
            fliplr(avg{regIdx,1}-err{regIdx,1})],...
            col(pairIdx*2-1,:),'linestyle','none','FaceAlpha',0.5)
        plot(t,avg{regIdx,1},'-','color',col(pairIdx*2-1,:))
        xlim([-1,2])
        ylim(yLim)
        ax=fixAxis;
        if pairIdx==1
            ylabel(yTxt)
        end
        if yType==2
            xlabel({'Time from' 'shock onset (s)'})
        end
        
        if ~isempty(sigOnset{regIdx})
            plot(t([sigOnset{regIdx};sigOffset{regIdx}])+tBinSize*[-1;1]/2,yLim(2)*0.95+[0,0],'k-','linewidth',0.5)
        end
        
        if yType==1
            title(regList{regIdx},'fontsize',5,'fontweight','normal')
        end
    end
end
end
%%
function panel_02(x,y)
%%
width=29-4;
xGapInter=5;

height=6;
hipHight=6;
yGapIntra=4;
yGapInter=6;

coact=poolVar('icaCoactTimeCondHT.mat');
react=poolVar('icaReacTimeCond.mat');
ses=poolVar('sessions.events.mat','base');
beh=poolVar('sleepState.states.mat','base');
ratList=fieldnames(coact)

tBinSize=3;
tWin=[-141,102];
tBorderZero=[fliplr(-(0:tBinSize:-min(tWin))),tBinSize:tBinSize:max(tWin)];
tBin=(tBorderZero(2:end)+tBorderZero(1:end-1))/2;

%%

col=setCoactColor();
for exIdx=1:2
    switch exIdx
        case 1

            rat=ratList{13};
            pairID=144;
        case 2
             rat=ratList{9};
            pairID=76;
    end
    idx=find(coact.(rat).pairID==pairID);
    tEvt{1}=coact.(rat).timestamp{idx};
    pEvt{1}=coact.(rat).peakHeight{idx};
    
    tEvt(2:3)=react.(rat).timestamps(coact.(rat).reacID(idx,:));
    
    pEvt(2:3)=react.(rat).peakHeight(coact.(rat).reacID(idx,:));
    
    tBorder=tBorderZero*60+ses.(rat).timestamps(2,2);
    
    reg=coact.(rat).region(idx,:);
    
    coactName=join(reg,'');
    coactName=strrep(strrep(coactName{1},' ',''),'/','');
    
    tempCol=[col.pair.(coactName);
        col.region.(strrep(strrep(reg{1},' ',''),'/',''));
        col.region.(strrep(strrep(reg{2},' ',''),'/',''))];
    for evtType=1:3
        cnt=histcounts(tEvt{evtType},[-inf,tBorder,inf]);
        evtBorder=cumsum(cnt);
        peakAvg=zeros(size(cnt,2)-2,1);
        peakSum=zeros(size(cnt,2)-2,1);
        
        for tIdx=1:length(evtBorder)-2
            peakAvg(tIdx)=mean(pEvt{evtType}(evtBorder(tIdx)+1:evtBorder(tIdx+1)));
        end
        rate=cnt(2:end-1)/tBinSize;
        
        for yType=1:2
            if yType==1
                val=rate;
            else
                val=peakAvg;
            end
            subplotInMM(x+(width+xGapInter)*(exIdx-1),...
                y+(height+yGapIntra)*(yType-1)+(2*height+yGapIntra+yGapInter)*(evtType-1),...
                width,height)
            bar(tBin,val,1,'linestyle','none','facecolor',tempCol(evtType,:))
            box off
            xlim(tBorderZero([1,end]))
            ax=fixAxis;
            hold on
            rectangle('Position',[-diff(ses.(rat).timestamps(2,:)/60),ax(3),diff(ses.(rat).timestamps(2,:)/60),diff(ax(3:4))],...
                'linestyle','none','facecolor',0.7*[1,1,1]);
            h=get(gca,'Children');
            isRec=ismember(h,findobj(h,'type','rectangle'));
            set(gca,'Children',[h(~isRec);h(isRec)])
            if yType==1
                if exIdx==1
                    ylabel({'Occurrence' 'rate' '(1/min)'})
                end
                if evtType==1
                    title(join(reg, ' - '),'fontsize',5,'fontweight','normal')
                else
                    title(reg{evtType-1},'fontsize',5,'fontweight','normal')
                end
            else
                if exIdx==1
                if evtType==1
                    ylabel({'Peak' '(z^2)'})
                else
                    ylabel({'Peak' '(z)'})
                end
                end
            end
        end
    end
    subplotInMM(x+(width+xGapInter)*(exIdx-1),...
        y+(2*height+yGapIntra)*3+yGapInter*2+yGapIntra,...
        width,hipHight)
    hipno=relabel_ma2sleep(beh.(rat).MECE.timestamps);
    hipno(:,1:2)=(hipno(:,1:2)-ses.(rat).timestamps(2,2))/60;
    hipno(hipno(:,2)<tBorderZero(1),:)=[];
    hipno(hipno(:,1)>tBorderZero(end),:)=[];
    if hipno(1,1)<tBorderZero(1);hipno(1,1)=tBorderZero(1);end
    if hipno(end,2)>tBorderZero(end);hipno(end,2)=tBorderZero(end);end
    for hipIdx=1:size(hipno,1)
        rectangle('Position',[hipno(hipIdx,1),(6-hipno(hipIdx,3))/2,diff(hipno(hipIdx,1:2)),1],...
                    'linestyle','none','facecolor','k')
    end
    xlim(tBorderZero([1,end]))
    ylim([0,4])
    if exIdx==1
    set(gca,'YTick',1:3,'YTickLabel',{'REM','NREM','WAKE'})
    else
        set(gca,'YTick',1:3,'YTickLabel',{''})
    end
    
    
    xlabel({'Time from' 'cond. session end (min)'})
end


end
%%
function panel_03(x,y)
width=18-8/3;
xGapIntra=4;

height=10;
yGapIntra=5;
yGapInter=14.5;
%%
basic=poolVar('basicMetaData.mat','base');

ica=poolVar('icaReacTimeCond.mat');

ses=poolVar('sessions.events.mat','base');
beh=poolVar('sleepState.states.mat','base');
ratList=fieldnames(basic);
%%


tRange{1}=[-102,0];
tRange{2}=[0,102];
tBinSize=3;

reg={};
sig=[];
for prePost=1:2
    tBorder{prePost}=-fliplr(0:tBinSize:-min(tRange{prePost}));
    tBorder{prePost}=[tBorder{prePost},tBinSize:tBinSize:max(tRange{prePost})];
    tBin{prePost}=(tBorder{prePost}(1:end-1)+tBorder{prePost}(2:end))/2;
    tBorder{prePost}=[-inf,tBorder{prePost},inf];
    
    cnt{prePost}=[];
    peak{prePost}=[];
    strength{prePost}=[];
    n=0;
    
    for rIdx=1:length(ratList);
        rat=ratList{rIdx};
        
        target=1:length(ica.(rat).region);
        
        slp=relabel_ma2sleep(beh.(rat).MECE.timestamps);
        slp(:,1:2)=slp(:,1:2);
 
        nrem=slp(slp(:,3)==3,1:2);
        if prePost==1
            t0=ses.(rat).timestamps(2,1);
        else
            t0=ses.(rat).timestamps(2,2);
        end
        
        algNrem=(nrem-t0)/60;
        
        tNrem=nan(1,length(tBin{prePost}));
        for tIdx=2:length(tBorder{prePost})-2
            subset=algNrem(algNrem(:,2)>tBorder{prePost}(tIdx)&algNrem(:,1)<tBorder{prePost}(tIdx+1),:);
            if isempty(subset); tNrem(tIdx-1)=0; continue; end
            if subset(1)<tBorder{prePost}(tIdx); subset(1)=tBorder{prePost}(tIdx); end
            if subset(end)>tBorder{prePost}(tIdx+1); subset(end)=tBorder{prePost}(tIdx+1); end
            tNrem(tIdx-1)=sum(diff(subset,1,2));
        end
        
        for idx=1:length(target)
            tEvt=ica.(rat).timestamps{target(idx)};
            pEvt=ica.(rat).peakHeight{target(idx)};
                        
            n=n+1;
            inNREM=any(tEvt>nrem(:,1) & tEvt<nrem(:,2));
            
            tEvt=(tEvt(inNREM)-t0)/60;
            pEvt=pEvt(inNREM);
            
            [tEvt,order]=sort(tEvt); % make sure it's sorted
            pEvt=pEvt(order);
            
            temp=histcounts(tEvt,tBorder{prePost});
            cnt{prePost}(n,:)=temp(2:end-1)./tNrem;
            
            evtIdx=cumsum(temp);
            
            for m=1:length(evtIdx)-2
                peak{prePost}(n,m)=nanmean(pEvt(evtIdx(m)+1:evtIdx(m+1)));
                strength{prePost}(n,m)=sum(pEvt(evtIdx(m)+1:evtIdx(m+1)))/tNrem(m);
            end
            if prePost==1
                reg(n,:)=ica.(rat).region(target(idx));
                sig(n)=~isempty(ica.(rat).partner.pos{target(idx)});
            end
        end
    end
end
[regPairList,~,pairID]=uniqueCellRows(reg);
%%

targetReg={'vCA1','BLA','PrL L5'};

colTempate=setCoactColor();
col=[
    0.5*[1,1,1];
    colTempate.region.vCA1;
    0.5*[1,1,1];
    colTempate.region.BLA;
    0.5*[1,1,1];
    colTempate.region.PrLL5;
    ];
%
yRange.cnt=[0,60;
    0,60
    0,60];
yRange.peak=[0,30;
    0,30;
    0,30];
yRange.str=[0,50;
    0,50;
    0,50];

for idx=1:length(targetReg)
    targetRegID=find(strcmpi(regPairList,targetReg{idx}));
    
    subSig=sig(pairID==targetRegID);
    for prePost=1:2
        for k=1:2
            subplotInMM(x+(width+xGapIntra)*(prePost-1),...
                y+(k-1)*(height+yGapIntra)+(2*height+yGapIntra+yGapInter)*(idx-1),width,height,true)
            switch k
                case 1
                    val=cnt{prePost}(pairID==targetRegID,:);
                    yTxt={'Occurrence' 'rate (1/min)'};
                    Ylim=yRange.cnt(idx,:);
                case 2
                    val=peak{prePost}(pairID==targetRegID,:);
                    yTxt={'Peak (z)'};
                    Ylim=yRange.peak(idx,:);
                case 3
                    val=strength{prePost}(pairID==targetRegID,:);
                    yTxt={'Strength' '(z/min)'};
                    Ylim=yRange.str(idx,:);
            end
            hold on
            np=[];
            targetBool{1}=(subSig~=1);
            targetBool{2}=(subSig==1);
            p=ones(1,size(val,2));
            for tIdx=1:size(val,2)
                if sum(~isnan(val(targetBool{1},tIdx)))>0 && sum(~isnan(val(targetBool{2},tIdx)))>0
                    p(tIdx)=ranksum(val(targetBool{1},tIdx),val(targetBool{2},tIdx));
                end
            end
            pSig=p<0.01;
            
            sigOnset=find(diff(pSig)==1)+1;
            if pSig(1); sigOnset=[1,sigOnset];end
            
            sigOffset=find(diff(pSig)==-1);
            if pSig(end); sigOffset=[sigOffset,length(pSig)];end
            
            
            for sigType=1:2
                target=targetBool{sigType};
                np(sigType)=sum(target);
                avg=nanmean(val(target,:),1);
                err=real(nanste(val(target,:),[],1));
                err(isnan(err))=0;
                avg(isnan(avg))=0;
                if sum(target)>1
                    patch([tBin{prePost},fliplr(tBin{prePost})],[avg+err,fliplr(avg-err)],col(2*(idx-1)+sigType,:),'linestyle','none','facealpha',0.5)
                end
                if sum(target)>0
                    plot(tBin{prePost},avg,'-','color',col(2*(idx-1)+sigType,:))
                end
                if sigType==2
                    ylim(Ylim)
                    xlim(tBorder{prePost}([2,end-1]))
                    if prePost==1
                        ylabel(yTxt)
                    else
                    end
                    if k==1
                        if prePost==1
                            tTxt=targetReg{idx};
                            textInMM(x+width+xGapIntra/2,...
                                y+(k-1)*(height+yGapIntra)-0.5+(2*height+yGapIntra+yGapInter)*(idx-1),tTxt,...
                                'fontsize',5,'fontweight','normal','verticalAlign','bottom','horizontalAlign','center')
                            
                        end
                    elseif k==2 
                        if prePost==1
                            xlabel({'Time to cond.' 'session start (min)'})
                        else
                            xlabel({'Time from cond.' 'session end (min)'})
                        end
                        
                    end
                    if ~isempty(sigOnset)
                        plot(tBin{prePost}([sigOnset;sigOffset])+tBinSize*[-1;1]/2,Ylim(2)*0.95+[0,0],'k-','linewidth',0.5)
                    end
                end
            end
        end
    end
end
end

