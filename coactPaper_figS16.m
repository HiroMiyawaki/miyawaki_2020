function coactPaper_figS16()
labelWeight=true;
labelCase=true;
labelSize=9;
letGapX=5;
letGapY=5;

close all
fh=initFig('height',13);



x=11;y=8;
panel_01(x,y);
panelLetter2(x-letGapX,y-letGapY,alphabet(1,labelCase),'fontSize',labelSize,'isBold',labelWeight)
drawnow();

x=11;y=8+35;
panel_02(x,y)
panelLetter2(x-letGapX,y-letGapY,alphabet(2,labelCase),'fontSize',labelSize,'isBold',labelWeight)
drawnow();

x=11+65;y=8+35;
panel_03(x,y);
panelLetter2(x-letGapX,y-letGapY,alphabet(3,labelCase),'fontSize',labelSize,'isBold',labelWeight)
drawnow();

print(fh,'figS16.pdf','-dpdf','-painters','-r300')

end
%%
%%

%%
function panel_01(x,y)

nCellWidth=18;
nCellHeight=15;
nCellGapXintra=6.5;

coact=poolVar('coactCompCell.mat');
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

%%
colTemp=setCoactColor();

piCol=[colTemp.cellType.inh;
       colTemp.cellType.nc;
       colTemp.cellType.ex];

targetReg={'PrL L5','BLA','vCA1'};


eiLeg{1}=sprintf('\\color[rgb]{%f %f %f}%s',colTemp.cellType.ex, 'Excitatory');
eiLeg{2}=sprintf('\\color[rgb]{%f %f %f}%s',colTemp.cellType.inh, 'Inhibitory');
eiLeg{3}=sprintf('\\color[rgb]{%f %f %f}%s',colTemp.cellType.nc, 'Not classified');


for n=1:2
    yTickPos.beh.PrLL5{n}=[];
    yTickPos.beh.vCA1{n}=[];
    yTickPos.beh.BLA{n}=[];

    yTickPos.nrem.PrLL5{n}=[];
    yTickPos.nrem.vCA1{n}=[];
    yTickPos.nrem.BLA{n}=[];
end

yTickPos.beh.BLA{1}=[0.5,1,2,5];
yTickPos.beh.vCA1{1}=[0.5,1,2,5];

yTickPos.nrem.vCA1{2}=[2,5,10,20];
%%
for targetIdx=1:length(targetReg)
    yTop=y;
    target=find(strcmp(reg,targetReg{targetIdx}));
    partnerList=targetReg;
    partnerList(strcmp(partnerList,targetReg{targetIdx}))=[]
    
    eiRatio=zeros(length(partnerList)+1,3);

    partnerName={};
    partnerNameCore={};
    pairCol=[];
    for n=1:length(partnerList)+1
        if n>length(partnerList)
            id=target(cellfun(@isempty,partner(target)));
            partnerName{n}='Others';
            partnerNameCore{n}='Others';
            pairCol(n,:)=colTemp.region.others;
        else
            id=target(cellfun(@(x) any(strcmp(x,partnerList{n})), partner(target)));
            partnerName{n}=['Coupled with ' partnerList{n}];
            partnerNameCore{n}=partnerList{n};
            pairCol(n,:)=colTemp.pair.(strrep(strrep([targetReg{targetIdx} partnerList{n}],' ',''),'/',''));
        end
        pairLeg{n+1}=sprintf('\\color[rgb]{%f %f %f}%s',pairCol(n,:),partnerName{n});
        eiId{1}=id(cellType(id)==1);
        eiId{2}=id(cellType(id)==-1);
        
        cnt=histcounts(cellType(id),-1.5:1.5);
        eiRatio(n,:)=cnt/sum(cnt)*100;

        fprintf('in %s, %s cells : n=%d\n',targetReg{targetIdx},partnerNameCore{n},sum(cnt))

    end

    regRatio=zeros(2,3);
    for eiIdx=1:2
        id=target(cellType(target)==3-eiIdx*2);
        cnt=histcounts(cellfun(@length,partner(id)),-0.5:2.5);
        regRatio(eiIdx,:)=cnt/sum(cnt)*100;
    end

    % # partner regions
    xShift=(targetIdx-1)*(nCellWidth+nCellGapXintra);
    yShift=0;
    subplotInMM(x+xShift,yTop+yShift,nCellWidth,nCellHeight)

    bar(0:2,regRatio','linestyle','none')
    colormap(gca,piCol([3,1],:))
    box off
    xlabel('# partner region')
    if targetIdx==1
        ylabel({'Fraction of' 'cells (%)'})
    end
    xlim([-0.5,2.5])
    title(targetReg{targetIdx},'fontsize',5,'fontweight','normal')
    ax=fixAxis;
    text2(0.7,1,eiLeg(1:2),ax,'verticalAlign','top')
end
end
function panel_02(x,y)

nCellHeight=15;

behWidth=21;
behHeight=16;
behGapXintra=4;

yGapIntraTop=14;

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

%%
hcIdx=[1,3,5,10,12];
behList={'wake','nrem','rem'};
for behIdx=1:3
    beh=behList{behIdx};
    behFR.(beh)=[];
    for ratIdx=1:length(ratList)
        rat=ratList{ratIdx};
        
        temp=meanFR.(rat).noDiv.Hz.(beh)(hcIdx,:);
        temp(isnan(temp))=0;
        dur=meanFR.(rat).noDiv.duration.(beh)(hcIdx);
        
        behFR.(beh)=[behFR.(beh),sum(temp.*dur',1)/sum(dur)];
    end
end
%%
colTemp=setCoactColor()

targetReg={'PrL L5','BLA','vCA1'};
behList={'wake','nrem','rem'};
cellTypeList={'excitatory cells','inhibitory cells'};


for n=1:2
    yTickPos.beh.PrLL5{n}=[];
    yTickPos.beh.vCA1{n}=[];
    yTickPos.beh.BLA{n}=[];

    yRange.beh.PrLL5{n}=[];
    yRange.beh.vCA1{n}=[];
    yRange.beh.BLA{n}=[];
end

yTickPos.beh.BLA{1}=[0.5,1,2,5];
yTickPos.beh.vCA1{1}=[0.5,1,2,5];

yRange.beh.PrLL5{1}=[1,4.2];
%%
for targetIdx=1:length(targetReg)
    yTop=y;
    target=find(strcmp(reg,targetReg{targetIdx}));
    partnerList=targetReg;
    partnerList(strcmp(partnerList,targetReg{targetIdx}))=[]
    
    eiFR.mean=zeros(length(partnerList)+1,length(behList),2);
    eiFR.ste=zeros(length(partnerList)+1,length(behList),2);
    eiFR.raw=cell(length(partnerList)+1,length(behList),2);
    
    partnerName={};
    partnerNameCore={};
    pairCol=[];
    for n=1:length(partnerList)+1
        if n>length(partnerList)
            id=target(cellfun(@isempty,partner(target)));
            partnerName{n}='Others';
            partnerNameCore{n}='Others';
            pairCol(n,:)=colTemp.region.others;
            pairLeg{2*n-1}=sprintf('\\color[rgb]{%f %f %f} %s', pairCol(n,:),partnerNameCore{n});
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
                
        for eiIdx=1:2           
            for behIdx=1:3
                beh=behList{behIdx};
                eiFR.mean(n,behIdx,eiIdx)=nanmean(log10(behFR.(beh)(eiId{eiIdx})));
                eiFR.ste(n,behIdx,eiIdx)=nanste(log10(behFR.(beh)(eiId{eiIdx})));
                eiFR.raw{n,behIdx,eiIdx}=(behFR.(beh)(eiId{eiIdx}));
            end            
        end
    end
    
    eiFR.p.all=ones(3,2);
    eiFR.p.each=ones(3,2,3,3);
    for behIdx=1:3
        for eiIdx=1:2
            temp=cellfun(@(x) ones(size(x)),eiFR.raw(:,behIdx,eiIdx),'UniformOutput',false);
            grp=[];
            for regIdx=1:size(temp,1)
                grp=[grp,regIdx*temp{regIdx}]
            end
            val=cat(2,eiFR.raw{:,behIdx,eiIdx});
            [p,~,s]=kruskalwallis(val,grp,'off');
            eiFR.p.all(behIdx,eiIdx)=p;

            idx=0;
            nGrp=3;
            for n=1:nGrp-1
                for m=n+1:nGrp
                    idx=idx+1;
                    p=ranksum(eiFR.raw{n,behIdx,eiIdx},eiFR.raw{m,behIdx,eiIdx});
                    p=p*nGrp*(nGrp-1)/2;
                    eiFR.p.each(behIdx,eiIdx,idx,:)=[n,m,p];
                end
            end

        end
    end

    % FR for each state
    yShift=(targetIdx-1)*(nCellHeight+yGapIntraTop);
    xShift=0;
    for eiIdx=1:2
        subplotInMM(x+xShift+(behWidth+behGapXintra)*(eiIdx-1),yTop+yShift,behWidth,behHeight)

        hold on
        posPool=[];
        poolAvg=[];
        for n=1:length(partnerName)
            avg=10.^eiFR.mean(n,:,eiIdx);
            pos=10.^(eiFR.mean(n,:,eiIdx)+eiFR.ste(n,:,eiIdx))-avg;
            neg=avg-10.^(eiFR.mean(n,:,eiIdx)-eiFR.ste(n,:,eiIdx));
            errorbar((1:3)+(n-2)*0.1,avg,neg,pos,...
                'linestyle','none','color',pairCol(n,:),'CapSize',0,'Marker','.','MarkerSize',4,'linewidth',0.5)
            posPool(n,:)=avg+pos;
            poolAvg(n,:)=avg;
        end
        for behIdx=1:3
            if eiFR.p.all(behIdx,eiIdx)<0.001
                sigTxt='###';
            elseif eiFR.p.all(behIdx,eiIdx)<0.01
                sigTxt='##';
            elseif eiFR.p.all(behIdx,eiIdx)<0.05
                sigTxt='#';
            else
                sigTxt='';
            end
            
            if ~isempty(sigTxt)
                text(behIdx+0.02,max(posPool(:,behIdx))*1.15,sigTxt,'horizontalAlign','center','fontsize',5,'verticalAlign','middle')
                
                if any(eiFR.p.each(behIdx,eiIdx,:,3)<0.05)
                    [sigPos,sigTxt]=findSig(squeeze(eiFR.p.each(behIdx,eiIdx,:,:)),poolAvg(:,behIdx));           
                    
                    for sigIdx=1:size(sigPos,1)
                        if sigPos(sigIdx,1)<0
                            sigX=[-0.18,-0.22,-0.22,-0.18];
                            sigTxtX=-0.23;
                            vAlign='middle';
                        else
                            sigX=[0.18,0.22,0.22,0.18];
                            sigTxtX=0.22;
                            vAlign='top';
                        end
                        sigY=sort(sigPos(sigIdx,3:4)).*[1.00,1.00];
                        
                        plot(behIdx+sigX,sigY([1,1,2,2]),'k-','linewidth',0.5)
                        text(behIdx+sigTxtX,geomean(sigPos(sigIdx,3:4)),sigTxt{sigIdx},'fontsize',6,...
                            'Rotation',90,'VerticalAlignment',vAlign,'HorizontalAlignment','center')
                    end
                    
                end
                
            
            end
        end

        
        set(gca,'YScale','log')
        set(gca,'XTick',1:3,'XTickLabel',upper(behList),'XTickLabelRotation',20)
        axis tight
        xlim([0.5,3.5])
        ax=fixAxis;
        if isempty(yRange.beh.(strrep(targetReg{targetIdx},' ','')){eiIdx})
            ylim(exp(log(ax(3:4))+diff(log(ax(3:4)))*[-1,1]/10))
        else
            ylim(yRange.beh.(strrep(targetReg{targetIdx},' ','')){eiIdx});
        end
        ax=fixAxis;

        tempTick=yTickPos.beh.(strrep(targetReg{targetIdx},' ','')){eiIdx};
        if ~isempty(tempTick)
            set(gca,'YTick',tempTick)
        end   
        
        title({targetReg{targetIdx},cellTypeList{eiIdx}},'fontsize',5,'fontweight','normal')
        if eiIdx==1
            ylabel('Firing rate (Hz)')
        else
        for n=1:5
            text2(1,1-0.15*(n-1),pairLeg{n},ax,'verticalALign','top')
        end
        end
    end
end
end
function panel_03(x,y)
yGapIntraBottom=14;
nremFRwidth=28;
nremFRheigth=15;
nremFRGapX=13;

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
   
    for eiIdx=2
        xShift=0;
        yShift=(targetIdx-1)*yGapIntraBottom+yExpandSum(eiIdx)*nremFRheigth;
        if targetIdx==1
            yExpand=1;
        else
            yExpand=1;
        end
        
        subplotInMM(x+xShift+(nremFRwidth+nremFRGapX)*(eiIdx-1)*0,yTop+yShift,nremFRwidth,nremFRheigth*yExpand)
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
        
            ylabel('Firing rate (Hz)')
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



