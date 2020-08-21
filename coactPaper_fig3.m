function coactPaper_fig3()
labelWeight=true;
labelCase=true;
labelSize=9;
letGapX=5;
letGapY=6;

close all
fh=initFig('height',10);

x=15;y=6;
panel_01_02(x,y);
panelLetter2(x-letGapX-2,y-letGapY+2,alphabet(1,labelCase),'fontSize',labelSize,'isBold',labelWeight)
panelLetter2(x-letGapX-2+52,y-letGapY+2,alphabet(2,labelCase),'fontSize',labelSize,'isBold',labelWeight)
drawnow();

x=15;y=33;
panel_03_04(x,y);
panelLetter2(x-letGapX-2,y-letGapY+1,alphabet(3,labelCase),'fontSize',labelSize,'isBold',labelWeight)
panelLetter2(x-letGapX-2+52,y-letGapY+1,alphabet(4,labelCase),'fontSize',labelSize,'isBold',labelWeight)
drawnow();

x=15;y=73;
panel_05(x,y);
panelLetter2(x-letGapX-2,y-letGapY,alphabet(5,labelCase),'fontSize',labelSize,'isBold',labelWeight)
drawnow();



print(fh,'fig03.pdf','-dpdf','-painters','-r300')
end

function panel_01_02(x,y)

width=19;
height=12;
xGapIntra=4;
xGapInter=12;

evtTrigReact=poolVar('shockTrigIcaCoact.mat');

ratList=fieldnames(evtTrigReact);
%%

t=evtTrigReact.(ratList{1}).time;
evtRate=[];
evtPeak=[];
evtStr=[];
reg={};
sig=[];
for ratIdx=1:length(ratList)
    rat=ratList{ratIdx};
    evtRate=cat(1,evtRate,evtTrigReact.(rat).rate);
    evtPeak=cat(1,evtPeak,evtTrigReact.(rat).peak);
    evtStr=cat(1,evtStr,evtTrigReact.(rat).strength);
    
    reg=cat(1,reg,evtTrigReact.(rat).region);
    sig=cat(1,sig,evtTrigReact.(rat).sigLevel);
end
%%

[regPairList,~,regPairIdx]=uniqueCellRows(reg);
%%
tBinSize=mean(diff(t));
smSigma=40/1000;
smBin=(0:ceil(smSigma*4/tBinSize))*tBinSize;
smBin=[-fliplr(smBin),smBin(2:end)]
smCore=normpdf(smBin,0,smSigma);
smCore=smCore/sum(smCore);

targetPair={'BLA','PrL L5';
    'vCA1','PrL L5'};

for yType=1:2
    avg={};
    err={};
    if yType==1
        val=evtRate;
        yTxt={'Occurrence' 'rate (1/s)'};
        yLim=[0,5];
    else
        val=evtPeak;
        yTxt='Peak (z^2)';
        yLim=[0,300];
    end
    
    for regIdx=1:size(regPairList,1);
        targetBool{1}=find(regPairIdx==regIdx&sig==1);
        targetBool{2}=find(regPairIdx==regIdx&sig==~1);
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
    
    col=[colDef.pair.BLAPrLL5
        0.3*[1,1,1]
        colDef.pair.vCA1PrLL5
        0.3*[1,1,1]];
    
    for pairIdx=1:size(targetPair,1)
        regIdx=find(strcmp(regPairList(:,1),targetPair{pairIdx,1})&strcmp(regPairList(:,2),targetPair{pairIdx,2}));
        
        subplotInMM(x+(width+xGapIntra)*(pairIdx-1)+(2*width+xGapIntra+xGapInter)*(yType-1),y,width,height)
        
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
        xlabel({'Time from' 'shock onset (s)'})
        
        if ~isempty(sigOnset{regIdx})
            plot(t([sigOnset{regIdx};sigOffset{regIdx}])+tBinSize*[-1;1]/2,yLim(2)*0.95+[0,0],'k-','linewidth',0.5)
        end
        
        pairName=join(regPairList(regIdx,:),' - ');
        title(pairName,'fontsize',5,'fontweight','normal')
    end
end
end

function panel_03_04(x,y)
width=18;
xGapIntra=4;
xGapInter=12;

height=10;
yGap=5;

%%
basic=poolVar('basicMetaData.mat','base');

ica=poolVar('icaCoactTimeCondHT.mat');

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
    % th=10
    
    for rIdx=1:length(ratList);
        rat=ratList{rIdx};
        
        target=find(cellfun(@(x,y) ~strcmp(x,y), ...
            ica.(rat).region(:,1),...
            ica.(rat).region(:,2)));
        
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
            tEvt=ica.(rat).timestamp{target(idx)};
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
                reg(n,:)=ica.(rat).region(target(idx),:);
                sig(n)=ica.(rat).sigLevel(target(idx));
            end
        end
    end
end
[regPairList,~,pairID]=uniqueCellRows(reg);
%%
targetPair={'BLA','PrL L5'
    'vCA1','PrL L5'};

colTempate=setCoactColor();
col=[
    0.5*[1,1,1];
    colTempate.pair.BLAPrLL5;
    0.5*[1,1,1];
    colTempate.pair.vCA1PrLL5;
    ];
%
yRange.cnt=[0,5;
    0,5];
yRange.peak=[0,400;
    0,400]
yRange.str=[0,700;
    0,500];

for idx=1:size(targetPair,1)
    targetPairID=find(strcmpi(regPairList(:,1),targetPair{idx,1})&strcmpi(regPairList(:,2),targetPair{idx,2}));
    
    subSig=sig(pairID==targetPairID);
    for prePost=1:2
        for k=1:2%3
            subplotInMM(x+(width+xGapIntra)*(prePost-1)+(2*width+xGapInter+xGapIntra)*(idx-1),...
                y+(k-1)*(height+yGap),width,height,true)
            switch k
                case 1
                    val=cnt{prePost}(pairID==targetPairID,:);
                    val2=cnt{3-prePost}(pairID==targetPairID,:);
                    yTxt={'Occurrence' 'rate (1/min)'};
                    Ylim=yRange.cnt(idx,:);
                case 2
                    val=peak{prePost}(pairID==targetPairID,:);
                    val2=peak{3-prePost}(pairID==targetPairID,:);
                    yTxt={'Peak (z^2)'};
                    Ylim=yRange.peak(idx,:);
                case 3
                    val=strength{prePost}(pairID==targetPairID,:);
                    val2=strength{3-prePost}(pairID==targetPairID,:);
                    yTxt={'Strength' '(z^2/min)'};
                    Ylim=yRange.str(idx,:);
            end
            if prePost==1
                
                temp=[nanmean(val(subSig==1,:));nanmean(val2(subSig==1,:))];
                temp2=join(yTxt);
                fprintf('\\Delta%s in %s - %s: %f +/- %f, p=%f\n',[temp2{:}],targetPair{idx,:},...
                    nanmean(diff(temp,1,1)),nanste(diff(temp,1,1)),signrank(temp(1,:),temp(2,:)))
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
                            tTxt=join(targetPair(idx,:), ' - ');
                            textInMM(x+width+xGapIntra/2+(2*width+xGapInter+xGapIntra)*(idx-1),...
                                y+(k-1)*(height+yGap)-0.5,tTxt,...
                                'fontsize',5,'fontweight','normal','verticalAlign','bottom','horizontalAlign','center')
                            
                        end
                    elseif k==2%3
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

function panel_05(x,y)
width=19*2+4;
xGapInter=12;
height=12;
sigCue=poolVar('icaReacZNCCGchamberCue_sig.mat');
sigChamber=poolVar('icaReacZNCCGchamber_sig.mat');
sigHomecage=poolVar('icaReacZNCCG_sig.mat');

tempSes=2;

ratList=fieldnames(sigHomecage);
preHC=[1,2,3,3,3,3,4];
col=[1,0,1;
    0.4*[1,1,1]];

sigJudgeIdx=preHC(tempSes)+1;
chID=[6:8,9,13];
chName={'Baseline','Conditioning','Context','Cue ses. before first tone','Cue ses. after first tone'};

sig=[];
reg={};

for ratIdx=1:length(ratList)
    rat=ratList{ratIdx};
    sig=[sig;sigHomecage.(rat)(tempSes).nrem.significance,sigChamber.(rat)(tempSes).significance,sigCue.(rat)(tempSes).significance];
    reg=[reg;sigHomecage.(rat)(tempSes).region(sigHomecage.(rat)(tempSes).pairID)];
end
tempName= sigHomecage.(ratList{1})(tempSes).template;

[regPairList,~,regID]=uniqueCellRows(reg);


targetPair={'BLA','PrL L5';
    'vCA1','PrL L5'};

colTempate=setCoactColor();
col=[
    colTempate.pair.BLAPrLL5;
    0.5*[1,1,1];
    colTempate.pair.vCA1PrLL5;
    0.5*[1,1,1];
    ];
yRange=[0,40;
    0,60];
for pairIdx=1:size(targetPair,1)
    subplotInMM(x+(width+xGapInter)*(pairIdx-1),y,width,height)
    
    pID=find(strcmp(regPairList(:,1),targetPair{pairIdx,1})&strcmp(regPairList(:,2),targetPair{pairIdx,2}));
    
    subSig=sig(regID==pID,:);
    grp=subSig(:,sigJudgeIdx);
    
    totN=[sum(grp==1);sum(grp~=1)];
    
    sigN=[sum(subSig(grp==1,chID)==1,1);
        sum(subSig(grp~=1,chID)==1,1)];
    
    for n=1:size(sigN,2)
        [~,p(n)]=fishertest([sigN(:,n),totN-sigN(:,n)]);
    end
    frac=[mean(subSig(grp==1,chID)==1,1)*100;
        mean(subSig(grp~=1,chID)==1,1)*100];
    
    bar(1:size(frac,2),frac','LineStyle','none')
    xlim([0,size(frac,2)+1])
    ylim(yRange(pairIdx,:))
    ax=fixAxis;
    if max(frac(:))>0.9*ax(4); ylim([ax(3),max(frac(:))+ax(4)*0.1]); end
    for n=1:size(sigN,2)
        if p(n)<0.001
            sTxt='***';
        elseif p(n)<0.01
            sTxt='**';
        elseif p(n)<0.05
            sTxt='*';
        else
            continue
        end
        text(n,max(frac(:,n))+diff(ax(3:4))*0.05,sTxt,'HorizontalAlignment','center','fontsize',7)
    end
    box off
    set(gca,'XTick',1:length(chName),'XTickLabel',chName,'XTickLabelRotation',20)
    colormap(gca,col(2*(pairIdx-1)+(1:2),:))
    title(join(regPairList(pID,:),' - '),'fontsize',5,'fontweight','normal');
    ylabel({'Fraction of pairs with' 'significant peak (%)'})
end
end


