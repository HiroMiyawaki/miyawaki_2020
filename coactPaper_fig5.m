function coactPaper_fig5()
labelWeight=true;
labelCase=true;
labelSize=9;
letGapX=7;
letGapY=5;
%
close all
fh=initFig('nColumn',2,'height',6.5);

x=6;y=5;
panel_01(x,y)
panelLetter2(x-letGapX+2,y-letGapY+1,alphabet(1,labelCase),'fontSize',labelSize,'isBold',labelWeight)
drawnow();

x=11;y=44;
panel_02_03(x,y);
panelLetter2(x-letGapX-3,y-letGapY,alphabet(2,labelCase),'fontSize',labelSize,'isBold',labelWeight)
panelLetter2(x-letGapX+52,y-letGapY,alphabet(3,labelCase),'fontSize',labelSize,'isBold',labelWeight)
drawnow();

x=102;y=44;
panel_04(x,y);
panelLetter2(x-letGapX-4,y-letGapY,alphabet(4,labelCase),'fontSize',labelSize,'isBold',labelWeight)
drawnow();


print(fh,'fig05.pdf','-dpdf','-painters','-r300')

end
%%
function panel_01(x,y)
nremFRwidth=23;
nremFRheigth=15;
nremFRGapX=17;

coact=poolVar('coactCompCell.mat');
meanFR=poolVar('meanFR.mat');
info=poolVar('okUnit.cellinfo.mat');


ratList=fieldnames(coact);

tempSes=2;
sigHC=3;
%% get partners
partner={};

for ratIdx=1:length(ratList)
    rat=ratList{ratIdx};
    
    
    region=relabel_region(coact.(rat).region,'minCellNum',0);
    
    temp=cell(size(info.(rat).channel));
    
    [cIdx,rIdx]=find(coact.(rat).ica(tempSes).homecage(sigHC).nrem);
    cList=unique(cIdx)';
    for cc=cList
        temp{cc}={region{ rIdx(cIdx==cc)}};
    end
    
    partner=[partner,temp];
end

for cIdx=1:length(partner)
    if isempty(partner{cIdx})
        continue
    end
    partner{cIdx}(~ismember(partner{cIdx},{'BLA','vCA1','PrL L5'}))=[];
    partner{cIdx}=unique(partner{cIdx});
end



%% get cellinfo
cellType=[];
reg=[];
for ratIdx=1:length(ratList)
    rat=ratList{ratIdx};
    
    cellType=[cellType;info.(rat).cellType.type'];
    reg=[reg,info.(rat).region];
end

[reg,regList]=relabel_region(reg,'minCellNum',0);
reg=reg';
%% getFR in HC
nremFR=[];

for ratIdx=1:length(ratList)
    rat=ratList{ratIdx};
    
    slp='nrem';
    nremFR=[nremFR;meanFR.(rat).Hz.(slp)([4:5,7:8],:)'];
    
end

%%
colTemp=setCoactColor();


targetReg={'PrL L5','BLA','vCA1'};
cellTypeList={'excitatory cells','inhibitory cells'};





for n=1:2


    yTickPos.nrem.PrLL5{n}=[];
    yTickPos.nrem.vCA1{n}=[];
    yTickPos.nrem.BLA{n}=[];
end

yTickPos.nrem.PrLL5{1}=[1,2,5];
yTickPos.nrem.vCA1{1}=[0.5,1,2,5,10,20];

yTickPos.nrem.vCA1{2}=[2,5,10,20];
%%
yExpandSum=[0,0];
for targetIdx=1:length(targetReg)
    yTop=y;
    target=find(strcmp(reg,targetReg{targetIdx}));
    partnerList=targetReg;
    partnerList(strcmp(partnerList,targetReg{targetIdx}))=[]
    
    eiNrem.mean=zeros(length(partnerList)+1,4,2);
    eiNrem.ste=zeros(length(partnerList)+1,4,2);
    eiNrem.raw=cell(length(partnerList)+1,4,2);
    
    partnerName={};
    partnerNameCore={};
    pairCol=[];
    for n=1:length(partnerList)+1
        if n>length(partnerList)
            id=target(cellfun(@isempty,partner(target)));
            partnerName{n}='Others';
            partnerNameCore{n}='Others';
            pairCol(n,:)=colTemp.region.others;
            pairLeg{2*n-1}=sprintf('\\color[rgb]{%f %f %f}%s',pairCol(n,:),'Others')
        else
            id=target(cellfun(@(x) any(strcmp(x,partnerList{n})), partner(target)));
            partnerName{n}=['Coupled with ' partnerList{n}];
            partnerNameCore{n}=partnerList{n};
            pairCol(n,:)=colTemp.pair.(strrep(strrep([targetReg{targetIdx} partnerList{n}],' ',''),'/',''));
            pairLeg{2*n-1}=sprintf('\\color[rgb]{%f %f %f}%s',pairCol(n,:),'Coupled');
            pairLeg{2*n}=sprintf('\\color[rgb]{%f %f %f}with %s', pairCol(n,:),partnerNameCore{n});
        end
        
        eiId{1}=id(cellType(id)==1);
        eiId{2}=id(cellType(id)==-1);
        
        cnt=histcounts(cellType(id),-1.5:1.5);
        eiRatio(n,:)=cnt/sum(cnt)*100;

        fprintf('in %s, %s cells : n=%d\n',targetReg{targetIdx},partnerNameCore{n},sum(cnt))
        
        for eiIdx=1:2           

            val=log10(nremFR(eiId{eiIdx},:));
            val(isinf(val))=nan;
            eiNrem.mean(n,:,eiIdx)=nanmean(val,1);
            eiNrem.ste(n,:,eiIdx)=nanste(val,[],1);
            for nremIdx=1:4
                eiNrem.raw{n,nremIdx,eiIdx}=nremFR(eiId{eiIdx},nremIdx)';
            end

        end
    end
 
    
    eiNrem.p.all=ones(4,2);
    eiNrem.p.each=ones(4,2,3,3);
    eiNrem.p.couple.all=ones(3,2);
    eiNrem.p.couple.each=ones(3,2,4*3/2,3);
    for nremIdx=1:4
        for eiIdx=1:2
            temp=cellfun(@(x) ones(size(x)),eiNrem.raw(:,nremIdx,eiIdx),'UniformOutput',false);
            grp=[];
            for regIdx=1:size(temp,1)
                grp=[grp,regIdx*temp{regIdx}]
            end
            val=cat(2,eiNrem.raw{:,nremIdx,eiIdx});
            [p,~,s]=kruskalwallis(val,grp,'off');
            eiNrem.p.all(nremIdx,eiIdx)=p;
            idx=0;
            nGrp=3;
            for n=1:nGrp-1
                for m=n+1:nGrp
                    idx=idx+1;
                    p=ranksum(eiNrem.raw{n,nremIdx,eiIdx},eiNrem.raw{m,nremIdx,eiIdx});
                    p=p*nGrp*(nGrp-1)/2;
                    eiNrem.p.each(nremIdx,eiIdx,idx,:)=[n,m,p];
                end
            end
        end
    end
    for eiIdx=1:2
        for n=1:3
            temp=cat(1,eiNrem.raw{n,:,eiIdx});
            temp(:,any(isnan(temp),1))=[];
            eiNrem.p.couple.all(n,eiIdx)=friedman(temp',1,'off');
            
            idx=0;
            for nremIdx1=1:3
                for nremIdx2=nremIdx1+1:4;
                    idx=idx+1;
                    p=signrank(temp(nremIdx1,:),temp(nremIdx2,:))*6
                    eiNrem.p.couple.each(n,eiIdx,idx,:)=[nremIdx1,nremIdx2,p];
                end
            end
        end
    end

    % FR within NREM
   
    for eiIdx=1%:2
        if targetIdx==1
            yExpand=1.5;
        else
            yExpand=1;
        end
        
        if yExpand==1
            yShift=nremFRheigth*0.5;
        else
            yShift=0;
        end
        xShift=(nremFRwidth+nremFRGapX)*(targetIdx-1)-(targetIdx==2)*1;    
        
        subplotInMM(x+xShift+(nremFRwidth+nremFRGapX)*(eiIdx-1),yTop+yShift,nremFRwidth,nremFRheigth*yExpand)
        yExpandSum(eiIdx)=yExpandSum(eiIdx)+yExpand;
        hold on
        posPool=[];
        poolAvg=[];
        for n=1:length(partnerName)
            avg=10.^eiNrem.mean(n,:,eiIdx);
            pos=10.^(eiNrem.mean(n,:,eiIdx)+eiNrem.ste(n,:,eiIdx))-avg;
            neg=avg-10.^(eiNrem.mean(n,:,eiIdx)-eiNrem.ste(n,:,eiIdx));
            errorbar(([1:2,3.5:4.5])+(n-2)*0.1,avg,neg,pos,...
                'color',pairCol(n,:),'CapSize',0,'Marker','.','MarkerSize',4,'linewidth',0.5,...
                'linestyle','none')
            posPool(n,:)=pos+avg;
            poolAvg(n,:)=avg;
        end
        
        for nremIdx=1:4
            if eiNrem.p.all(nremIdx,eiIdx)<0.001
                sigTxt='###';
            elseif eiNrem.p.all(nremIdx,eiIdx)<0.01
                sigTxt='##';
            elseif eiNrem.p.all(nremIdx,eiIdx)<0.05
                sigTxt='#';
            else
                sigTxt='';
            end
            
            if ~isempty(sigTxt)
                text(nremIdx+(nremIdx>2)*0.5,max(posPool(:,nremIdx))*1.1,sigTxt,'horizontalAlign','center','fontsize',5,'verticalAlign','middle')

            
                if any(eiNrem.p.each(nremIdx,eiIdx,:,3)<0.05)
                 [sigPos,sigTxt]=findSig(squeeze(eiNrem.p.each(nremIdx,eiIdx,:,:)),poolAvg(:,nremIdx));   
                 
                    for sigIdx=1:size(sigPos,1)
                        if sigPos(sigIdx,1)<0
                            sigX=[-0.22,-0.28,-0.28,-0.222];
                            sigTxtX=-0.31;
                            vAlign='middle';
                        else
                            sigX=[0.22,0.28,0.28,0.22];
                            sigTxtX=0.28;
                            vAlign='top';
                        end
                        sigY=sort(sigPos(sigIdx,3:4)).*[1.00,1.00];
                        
                        plot(nremIdx+(nremIdx>2)*0.5+sigX,sigY([1,1,2,2]),'k-','linewidth',0.5)
                        text(nremIdx+(nremIdx>2)*0.5+sigTxtX,geomean(sigPos(sigIdx,3:4)),sigTxt{sigIdx},'fontsize',6,...
                            'Rotation',90,'VerticalAlignment',vAlign,'HorizontalAlignment','center')
                    end
                    
                end
            end
                
        end
        sigPosY=max(posPool(:));
        sigPosX=[1:2,3.5:4.5];
        curPos=0;
        sigTxtPool=cell(1,3);
        for n=length(partnerName):-1:1
            sigTxt=getSigTxt(eiNrem.p.couple.all(n,eiIdx));
            if ~isempty(sigTxt)
                sigTxtPool{n}=repmat('\dag',1,length(sigTxt));
                
                [sigPos,sigTxt]=findSigPos(squeeze(eiNrem.p.couple.each(n,eiIdx,:,:)));
                if ~isempty(sigPos)
                    for idx=1:size(sigPos,1)
                        plot(sigPosX(sigPos(idx,[1,1,2,2])), sigPosY*1.2^(sigPos(idx,3)+curPos)*[1,1.05,1.05,1],'-',...
                            'LineWidth',0.5,'color',pairCol(n,:))
                        text(mean(sigPosX(sigPos(idx,[1,2]))),sigPosY*1.2^(sigPos(idx,3)+curPos)*1.025,sigTxt{idx},...
                            'HorizontalAlignment','center','VerticalAlignment','baseline', 'Color',pairCol(n,:))
                    end                
                    curPos=curPos+max(sigPos(:,3));
                end
            end
            
        end
        
        set(gca,'XTick',[1:2,3.5:4.5],'XTickLabel',{'First','Last'},'XTickLabelRotation',0)
        set(gca,'YScale','log')
        title([targetReg{targetIdx},' ' cellTypeList{eiIdx}],'fontsize',5,'fontweight','normal')
        axis tight
        xlim([0.5,5.2])
        ax=fixAxis;
                
        ylim(exp(log(ax(3:4))+diff(log(ax(3:4)))*[-1,1]/10))
        ax=fixAxis;             

        tempTick=yTickPos.nrem.(strrep(targetReg{targetIdx},' ','')){eiIdx};
        if ~isempty(tempTick)
            set(gca,'YTick',tempTick)
        end        
        
        if eiIdx==1
            ylabel('Firing rate (Hz)')
        end
        for n=1:5
            text2(1,1/yExpand-0.15*(n-1)/yExpand,pairLeg{n},ax,'verticalALign','top')
        end
        for n=1:3
            if ~isempty(sigTxtPool{n})
                text2(5/5.2,1/yExpand-0.15*(2*n-2)/yExpand,sigTxtPool{n},ax,'verticalALign','top',...
                    'color',pairCol(n,:),'Interpreter','latex','horizontalALign','center')
            end
        end
        text2(1/4.7,-0.25/yExpand,'Pre-cond.',ax,'horizontalAlign','center','verticalALign','top')
        text2(1/4.7,-0.37/yExpand,'NREM',ax,'horizontalAlign','center','verticalALign','top')
        text2(3.5/4.7,-0.25/yExpand,'Post-cond.',ax,'horizontalAlign','center','verticalALign','top')
        text2(3.5/4.7,-0.37/yExpand,'NREM',ax,'horizontalAlign','center','verticalALign','top')
        
    end

end
end

function panel_02_03(x,y)
width=16;
height=10;
xGap=11;
yGap=10;
%%
coact=poolVar('coactCompCell.mat');

info=poolVar('okUnit.cellinfo.mat');
rip=poolVar('ripFrMod.mat');
spdl=poolVar('spdlFrMod.mat');
delta=poolVar('deltaFrMod.mat');
hfo=poolVar('hfoFrMod.mat');

ratList=fieldnames(coact);


tempSes=2;
sigHC=3;

%% get partners
partner={};

for ratIdx=1:length(ratList)
    rat=ratList{ratIdx};
    
    
    region=relabel_region(coact.(rat).region,'minCellNum',0);
    
    temp=cell(size(info.(rat).channel));
    
    [cIdx,rIdx]=find(coact.(rat).ica(tempSes).homecage(sigHC).nrem);
    cList=unique(cIdx)';
    for cc=cList
        temp{cc}={region{ rIdx(cIdx==cc)}};
    end
    
    partner=[partner,temp];
    
end

%% get cellinfo
cellType=[];
reg=[];
for ratIdx=1:length(ratList)
    rat=ratList{ratIdx};
    
    cellType=[cellType;info.(rat).cellType.type'];
    reg=[reg,info.(rat).region];
end

[reg,regList]=relabel_region(reg,'minCellNum',0);
reg=reg';

%% get param in HC

hc.nSwrPart=[];
hc.nSwrGain=[];

for ratIdx=1:length(ratList)
    rat=ratList{ratIdx};
    
    if isfield(rip,rat)
        hc.nSwrPart=[hc.nSwrPart;cat(1,rip.(rat).hcNrem(2:3).participation)'*100];
        hc.nSwrGain=[hc.nSwrGain;cat(1,rip.(rat).hcNrem(2:3).gain)'*100];
        
    else
        temp=nan(size(info.(rat).channel,2),2);
        hc.nSwrPart=[hc.nSwrPart;temp];
        hc.nSwrGain=[hc.nSwrGain;temp];
    end
end
%%
hc.hfoPart=[];
hc.hfoGain=[];

for ratIdx=1:length(ratList)
    rat=ratList{ratIdx};
    
    hc.hfoPart=[hc.hfoPart;cat(1,hfo.(rat).hcNrem(2:3).participation)'*100];
    hc.hfoGain=[hc.hfoGain;cat(1,hfo.(rat).hcNrem(2:3).gain)'*100];
end
%%

hc.spdlPart=[];
hc.spdlGain=[];
for ratIdx=1:length(ratList)
    rat=ratList{ratIdx};
    
    hc.spdlPart=[hc.spdlPart;cat(1,spdl.(rat).pfc(2:3).participation)'*100];
    hc.spdlGain=[hc.spdlGain;cat(1,spdl.(rat).pfc(2:3).gain)'*100];
end

hc.preDelPart=[];
hc.preDelGain=[];
hc.postDelPart=[];
hc.postDelGain=[];
hc.delSpkPart=[];
for ratIdx=1:length(ratList)
    rat=ratList{ratIdx};
    
    hc.delSpkPart=[hc.delSpkPart;cat(1,delta.(rat).pfc.deltaSpike(2:3).participation)'*100];
    hc.preDelPart=[hc.preDelPart;cat(1,delta.(rat).pfc.preDelta(2:3).participation)'*100];
    hc.postDelPart=[hc.postDelPart;cat(1,delta.(rat).pfc.postDelta(2:3).participation)'*100];
    hc.preDelGain=[hc.preDelGain;cat(1,delta.(rat).pfc.preDelta(2:3).gain)'*100];
    hc.postDelGain=[hc.postDelGain;cat(1,delta.(rat).pfc.postDelta(2:3).gain)'*100];
end

%%
colTemp=setCoactColor();

col=[colTemp.cellType.inh
    colTemp.cellType.nc
    colTemp.cellType.ex];
%%
CTname={'Inhibitory','Not classified','Excitatory'};

mesList={'nSwrGain','hfoGain','delSpkPart'};

mesNameList.nSwrPart={'SWR' 'participation (%)'};
mesNameList.nSwrGain={'FR modulation' 'by SWRs (%)'};
mesNameList.wSwrPart={'Wake SWR' 'participation (%)'};
mesNameList.wSwrGain={'Wake SWR' 'gain (%)'};
mesNameList.hfoPart={'HFO' 'participation (%)'};
mesNameList.hfoGain={'FR modulation' 'by HFOs (%)'};
mesNameList.spdlPart={'spindle' 'participation (%)'};
mesNameList.spdlGain={'spindle' 'gain (%)'};
mesNameList.delSpkPart={'Delta spike' 'participation (%)'};
mesNameList.preDelPart={'Delta onset' 'participation (%)'};
mesNameList.postDelPart={'Delta offset' 'participation (%)'};
mesNameList.preDelGain={'Delta onset' 'gain (%)'};
mesNameList.postDelGain={'Delta offset' 'gain (%)'};

targetReg={'PrL L5'};
partnerList={'PrL L5','BLA','vCA1'};

yRange.PrLL5.nSwrGain=[0,200];
yRange.PrLL5.hfoGain=[0,250];
yRange.PrLL5.delSpkPart=[0,4];

cellLeg={};
for n=1:3
    cellLeg{n}=sprintf('\\color[rgb]{%f %f %f}%s',col(n,:),CTname{n})
end

for regIdx=1:length(targetReg)
    
    target=find(strcmp(reg,targetReg{regIdx}));
    
    patList=partnerList;
    patList(strcmp(patList,targetReg{regIdx}))=[];
    cList={};
    
    for n=1:length(patList)
        cList{n}=target(cellfun(@(x) any(strcmp(patList{n},x)),partner(target)));
    end
    cName=cellfun(@(x) ['' x], patList,'UniformOutput',false);
    
    
    for mesIdx=1:length(mesList)
        mesName=mesList{mesIdx};
        subplotInMM(x+(width+xGap)*(mesIdx-1)+(mesIdx==3)*1,y,width,height)
        hold on
            
        fprintf('\n%s\n',[mesNameList.(mesName){:}])
            
        for CTidx=0:1
            fprintf('%s ',CTname{3-n*CTidx})
            avg=[];
            err=[];
            p=[];
            for n=1:length(cName)
                val=hc.(mesName)(cList{n},:);
                subCT=cellType(cList{n});
                subset=val(subCT==(1-2*CTidx),:);
                
                avg(n,:)=nanmean(subset,1);
                if size(subset,1)>1
                    err(n,:)=nanste(subset,[],1);
                    if any(all(~isnan(subset),2))
                        p(n)=signrank(subset(:,1),subset(:,2));
                    else
                        p(n)=nan;
                    end
                else
                    err(n,:)=[0,0];
                    p(n)=nan;
                end
            end
            hold on
            
            xVal{1}=(0:length(cName)-1)*(length(cName)+0.5)+1-0.2+CTidx;
            plot(xVal{1}+[0;0],(avg(:,1)+err(:,1).*[0,1])','Color',col(3-2*CTidx,:))
            bar(xVal{1},avg(:,1),0.1,'EdgeColor',col(3-2*CTidx,:),'facecolor',0.99*[1,1,1])
            
            xVal{2}=(0:length(cName)-1)*(length(cName)+0.5)+1+0.2+CTidx;
            plot(xVal{2}+[0;0],(avg(:,2)+err(:,2).*[0,1])','Color',col(3-2*CTidx,:))
            bar(xVal{2},avg(:,2),0.1,'EdgeColor',col(3-2*CTidx,:),'facecolor',col(3-2*CTidx,:))
            ax=fixAxis;
            for n=1:length(p)
                fprintf('%s:p=%f ',cName{n},p(n))
                if p(n)<0.001
                    sig='***';
                elseif p(n)<0.01
                    sig='**';
                elseif p(n)<0.05
                    sig='*';
                else
                    continue
                end
                sigPosY=max((avg(n,:)+err(n,:)))+diff(yRange.(strrep(targetReg{regIdx},' ','')).(mesName))/20;
                plot((length(cName)+0.5)*(n-1)+CTidx+1+0.2*[-1,1],sigPosY+[0,0],'k-')
                text((length(cName)+0.5)*(n-1)+CTidx+1,sigPosY,sig,'HorizontalAlignment','center','fontsize',7)
            end
            fprintf('\n')
  
        end
        
        ylabel(mesNameList.(mesName))
        set(gca,'XTick',[])
        xlim([0.5,length(cName)*2+0.5+0.5])
        ylim(yRange.(strrep(targetReg{regIdx},' ','')).(mesName))
        ax=fixAxis;

        for cIdx=1:length(cName)
            text((length(cName)+0.5)*(cIdx-1)+1+0.5,ax(3)-diff(ax(3:4))*0.075,...
                cName{cIdx},'horizontalAlign','center','fontsize',5,'verticalAlign','top')
        end
        text(mean([0.5,length(cName)*2+0.5+0.5]),ax(3)-diff(ax(3:4))*0.25,'Coupled region',...
            'horizontalAlign','center','fontsize',5,'verticalAlign','top')
        
        title(['Cells in ' targetReg{regIdx}],'fontsize',5,'fontweight','normal')

        
    end    
    for idx=3
        subplotInMM(x+(width+xGap)*(idx-1)+width+(idx==3)*1+0.5,y,5,height)
        xlim([0,5])
        ylim([0,height])
        text(0,height,cellLeg{3},'verticalAlign','middle','fontsize',5)
        text(0,height-2,cellLeg{1},'verticalAlign','middle','fontsize',5)

        rectangle('Position',[0,height-5.25,2,1],'linestyle','-','EdgeColor','k')
        text(2.5,height-4,'Pre-','verticalAlign','middle','fontsize',5)    
        text(2.5,height-5.5,'cond.','verticalAlign','middle','fontsize',5)    

        rectangle('Position',[0,height-8.75,2,1],'linestyle','-','EdgeColor','k','FaceColor','k')
        text(2.5,height-7.5,'Post-','verticalAlign','middle','fontsize',5)    
        text(2.5,height-9,'cond.','verticalAlign','middle','fontsize',5)   

        axis off
    end
    
end

%%


end

%%
function panel_04(x,y)
width=16;
height=10;

%%
coact=poolVar('coactCompCell.mat');
info=poolVar('okUnit.cellinfo.mat');
cue=poolVar('cueTrigSpk.mat');

ratList=fieldnames(coact);


tempSes=2;
sigHC=3;
%% get partners
partner={};

for ratIdx=1:length(ratList)
    rat=ratList{ratIdx};
        region=relabel_region(coact.(rat).region,'minCellNum',0);        
        temp=cell(size(info.(rat).channel));
        
        [cIdx,rIdx]=find(coact.(rat).ica(tempSes).homecage(sigHC).nrem);
        cList=unique(cIdx)';
        for cc=cList
            temp{cc}={region{ rIdx(cIdx==cc)}};
        end
        
        partner=[partner,temp];
end

%% get cellinfo
cellType=[];
reg=[];
for ratIdx=1:length(ratList)
    rat=ratList{ratIdx};
    
    cellType=[cellType;info.(rat).cellType.type'];
    reg=[reg,info.(rat).region];
end

[reg,regList]=relabel_region(reg,'minCellNum',0);
reg=reg';

%%
tBin=cue.(ratList{1}).tBin;

tWinShock=2;
tWinCue=2;

baseBin=find(tBin>-tWinCue&tBin<0);
onsetBin=find(tBin>0&tBin<tWinCue);

befShockBin=find(tBin>30-tWinShock&tBin<30);
aftShockBin=find(tBin>32&tBin<32+tWinShock);

for sesIdx=1:5
    cueFR{sesIdx}=[];
end
shockFR=[];
for ratIdx=1:length(ratList)
    rat=ratList{ratIdx};
    for sesIdx=1:5
        cueFR{sesIdx}=[cueFR{sesIdx};
            [mean(cue.(rat).FR{sesIdx}(baseBin,:));mean(cue.(rat).FR{sesIdx}(onsetBin,:))]'];
    end
    shockFR=[shockFR;
        [mean(cue.(rat).FR{2}(befShockBin,:));mean(cue.(rat).FR{2}(aftShockBin,:))]'];
end


mi=@(x) diff(x,1,2)./sum(x,2);
cellRegList={'PrL L5','vCA1','BLA'};
chName={'Base','Cond','C-Tone','C-Ext','R-Ext'};

for cellRegIdx=1:length(cellRegList)
    cellReg=cellRegList{cellRegIdx};
    
    partnerRegList={};
    for partnerIdx=1:length(cellRegList)
        if ~strcmp(cellReg,cellRegList{partnerIdx})
            partnerRegList{end+1}=cellRegList{partnerIdx}
        end
    end
    
    for partnerIdx=1:2
        partnerReg=partnerRegList{partnerIdx};
        cIdx=find(strcmp(reg,cellReg)&cellfun(@(x) any(strcmp(x,partnerReg)),partner)');
        nonPartnerReg=partnerRegList(~ismember(partnerRegList,partnerReg));
        
        for sesIdx=1:6
            if sesIdx==6
                ciSub(cellRegIdx,sesIdx).ci{partnerIdx}=(mi(shockFR(cIdx,:)));
                ciSub(cellRegIdx,sesIdx).trigName='Shock';
            else
                ciSub(cellRegIdx,sesIdx).ci{partnerIdx}=(mi(cueFR{sesIdx}(cIdx,:)));
                ciSub(cellRegIdx,sesIdx).trigName=['Cue onset in ' chName{sesIdx}];
            end
            ciSub(cellRegIdx,sesIdx).cellReg=cellReg;
            ciSub(cellRegIdx,sesIdx).partner{partnerIdx}=partnerReg;
            ciSub(cellRegIdx,sesIdx).cType{partnerIdx}=cellType(cIdx);
        end
    end
end

%%
colTemp=setCoactColor();

col=[
    colTemp.cellType.ex
    colTemp.cellType.inh
    ];
CTname={'Excitatory','Inhibitory'};
cellLeg={};
for n=1:2
    cellLeg{n}=sprintf('\\color[rgb]{%f %f %f}%s',col(n,:),CTname{n})
end

ylimList=[-0.6,0.5;
          -1,0.5;
          -0.5,0.4];

subplotInMM(x,y,width,height)
regIdx=1;
sesIdx=6;
        hold on
        for cType=1:2
            data={};
            for patIdx=1:2
                data{patIdx}=ciSub(regIdx,sesIdx).ci{patIdx}(ciSub(regIdx,sesIdx).cType{patIdx}==3-2*cType);
            end
            avg=cellfun(@nanmean,data);
            err=cellfun(@nanste,data);
            
            for n=1:length(data)
                if isempty(data{n})
                    pZero(n)=1;
                else
                    pZero(n)=signrank(data{n});
                end
            end
            
            if all(cellfun(@(x) sum(~isnan(x)),data)>0)
                p=ranksum(data{:});
            else
                p=1;
            end
            
            x=(1:2)+(cType-1)*3;
            
            bar(x,avg,'LineStyle','none','FaceColor',col(cType,:))
            errBar=avg+(2*(avg>0)-1).*[0,0;err];
            
            for patIdx=1:2
                plot(x(patIdx)+[0,0],errBar(:,patIdx),'-','color',col(cType,:))
            end
            
            ylim(ylimList(regIdx,:))
                        
            if p<0.001
                sigTxt='***';
            elseif p<0.01
                sigTxt='**';
            elseif p<0.05
                sigTxt='*';
            else
                sigTxt='';
            end
            if ~isempty(sigTxt)
                
                plot(x,max([0;errBar(:)])+diff(ylimList(regIdx,:))/20+[0,0],'k-')
                text(mean(x),max([0;errBar(:)])+diff(ylimList(regIdx,:))/20,sigTxt,'HorizontalAlignment','center','FontSize',7)
            end
            
        end
        xlim([0,6])
        ax=fixAxis;

        xTickPos=[1,2,4,5];
        set(gca,'XTick',xTickPos,'XTickLabel',[])
        for n=1:length(xTickPos)
            text(xTickPos(n),ax(3)-diff(ax(3:4))*0.025,ciSub(regIdx,sesIdx).partner{mod(n-1,2)+1},...
                'verticalAlign','top','horizontalAlign','right','rotation',35,'fontsize',5)
        end
        ylabel({'Shock' 'modulation index'})
        text2(0.5,-0.4,'Coupled region',ax,'verticalAlign','top','horizontalAlign','center','fontsize',5)

        text2(0.5,0.4,cellLeg{1},ax,'verticalAlign','top','fontsize',5)
        text2(0.5,0.25,cellLeg{2},ax,'verticalAlign','top','fontsize',5)
        title(['Cells in ' ciSub(regIdx,sesIdx).cellReg],'fontsize',5,'fontweight','normal')

end

%%
function [sigPos,sigTxt]=findSig(pVal,pos)
    [sortedPos,posIdx]=sort(pos);
    posIdx(posIdx)=1:length(posIdx)
    
    sigPos=[];
    sigTxt={};
    tempPos=-1;
    sCnt=0;
    for n=1:3
        switch n
            case 1
                n=find(posIdx==1);
                m=find(posIdx==2);
            case 2
                n=find(posIdx==2);
                m=find(posIdx==3);
            case 3        
                n=find(posIdx==1);
                m=find(posIdx==3);
                if length(sigTxt)>0
                    tempPos=tempPos+2;
                end
        end
        idx=find((pVal(:,1)==n & pVal(:,2)==m) |  (pVal(:,1)==m & pVal(:,2)==n))
        if (pVal(idx,3)<0.05);
            sCnt=sCnt+1;
            sigPos(sCnt,:)=[tempPos,tempPos,pos(pVal(idx,1:2))'];
            sigTxt{sCnt}=getSigTxt(pVal(idx,3));
        end
    end    
end

function sigTxt=getSigTxt(p)
    if p<0.001
        sigTxt='***';
    elseif p<0.01
        sigTxt='**';
    elseif p<0.05
        sigTxt='*';
    else 
        sigTxt='';
    end
end

function [sigPos,sigTxt]=findSigPos(pVal)
    [~,order]=sort(diff(pVal(:,1:2),1,2));
    empty=true(1,2*max(max(pVal(:,1:2))));
    sCnt=0;
    sigPos=[];
    sigTxt={};
    for n=order'
        if pVal(n,3)>=0.05
            continue
        end
         sCnt=sCnt+1;
         level=1;
         
         while ~all(empty(level,2*min(pVal(n,1:2)):2*max(pVal(n,1:2))-1))
             level=level+1;
             if size(empty,1)<level
                 empty(level,:)=true(1,2*max(max(pVal(:,1:2))));
                 break
             end
         end
         empty(level,2*min(pVal(n,1:2)):2*max(pVal(n,1:2))-1)=false;
         sigPos(sCnt,:)=[pVal(n,1:2),level*[1,1]];
         sigTxt{sCnt}=getSigTxt(pVal(n,3));
    end
end
