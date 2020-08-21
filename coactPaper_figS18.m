function coactPaper_figS18()
labelWeight=true;
labelCase=true;
labelSize=9;
letGapX=9;
letGapY=3;

close all
fh=initFig('height',8);


x=11;y=5;
panel_01(x,y);
panelLetter2(x-letGapX,y-letGapY-2,alphabet(1,labelCase),'fontSize',labelSize,'isBold',labelWeight)
drawnow();

x=11+28;y=5;
panel_02(x,y);
panelLetter2(x-letGapX,y-letGapY,alphabet(2,labelCase),'fontSize',labelSize,'isBold',labelWeight)
drawnow();

print(fh,'figS18.pdf','-dpdf','-painters','-r300')

end

%%
function panel_01(x,y)
width=16;
height=14;
yGap=10;
xGap=9;
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
    0,0.4];


for regIdx=2:3;
    sesIdx=6;
    subplotInMM(x,y+(height+yGap)*(regIdx-2),width,height)
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
        
        xVal=(1:2)+(cType-1)*3;
        
        bar(xVal,avg,'LineStyle','none','FaceColor',col(cType,:))
        errBar=avg+(2*(avg>0)-1).*[0,0;err];
        
        for patIdx=1:2
            plot(xVal(patIdx)+[0,0],errBar(:,patIdx),'-','color',col(cType,:))
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
            
            plot(xVal,max([0;errBar(:)])+diff(ylimList(regIdx,:))/20+[0,0],'k-')
            text(mean(xVal)+0.0375,max([0;errBar(:)])+diff(ylimList(regIdx,:))/20,sigTxt,'HorizontalAlignment','center','FontSize',6)
        end
        
    end
    xlim([0,6])
    ax=fixAxis;
    xTickPos=[1,2,4,5];
    set(gca,'XTick',xTickPos,'XTickLabel',[])
    for n=1:length(xTickPos)
        text(xTickPos(n),ax(3)-diff(ax(3:4))*0.025,ciSub(regIdx,sesIdx).partner{mod(n-1,2)+1},...
            'verticalAlign','top','horizontalAlign','right','rotation',40,'fontsize',5)
    end
    
    if regIdx==3
        text2(0.5,-0.4,'Coupled region',ax,'verticalAlign','top','horizontalAlign','center','fontsize',5)
    end
    ylabel({'Shock' 'modulation index'},'fontsize',5)
    title(['Cells in ' ciSub(regIdx,sesIdx).cellReg],'fontsize',5,'fontweight','normal')
    
end

end
%%
function panel_02(x,y)
width=16;
height=14;
yGap=10;
xGap=11;
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
chName={'baseline','conditioning','cue retention','C-Ext','R-Ext'};

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

ylimList=[-0.3,0.5;
    -0.2,0.5;
    -0.5,0.4];

rowOrder=[3,1,2];
for regIdx=1:3;
    for sesIdx=1:3;
        subplotInMM(x+(width+xGap)*(sesIdx-1),y+(height+yGap)*(rowOrder(regIdx)-1),width,height)
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
            
            xVal=(1:2)+(cType-1)*3;
            
            bar(xVal,avg,'LineStyle','none','FaceColor',col(cType,:))
            errBar=avg+(2*(avg>0)-1).*[0,0;err];
            
            for patIdx=1:2
                plot(xVal(patIdx)+[0,0],errBar(:,patIdx),'-','color',col(cType,:))

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
                
                plot(xVal,max([0;errBar(:)])+diff(ylimList(regIdx,:))/20+[0,0],'k-')
                text(mean(xVal),max([0;errBar(:)])+diff(ylimList(regIdx,:))/20,sigTxt,'HorizontalAlignment','center','FontSize',7)
            end
            
        end
        xlim([0,6])
        ax=fixAxis;
 
        xTickPos=[1,2,4,5];
        set(gca,'XTick',xTickPos,'XTickLabel',[])
        for n=1:length(xTickPos)
            text(xTickPos(n),ax(3)-diff(ax(3:4))*0.025,ciSub(regIdx,sesIdx).partner{mod(n-1,2)+1},...
                'verticalAlign','top','horizontalAlign','right','rotation',40,'fontsize',5)
        end
        if  rowOrder(regIdx)==2
            ylabel(sprintf('Cue modulation index in %s session',chName{sesIdx}))
        end
        if rowOrder(regIdx)==3
            text2(0.5,-0.4,'Coupled region',ax,'verticalAlign','top','horizontalAlign','center','fontsize',5)
        end
        if rowOrder(regIdx)==1&sesIdx==3
            text2(1.1,1,cellLeg,ax,'verticalAlign','top','fontsize',5)
        end
        title(['Cells in ' ciSub(regIdx,sesIdx).cellReg],'fontsize',5,'fontweight','normal')
        
    end
end


end



