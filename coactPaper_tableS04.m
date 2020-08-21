function coactPaper_tableS04()
labelWeight=true;
labelCase=true;
labelSize=9;
letGapX=5;
letGapY=5;

close all
fh=initFig('height',3);


x=10;y=2;
panel_01(x,y);
drawnow();

print(fh,'tableS04.pdf','-dpdf','-painters','-r300')

end
%%
%%

function panel_01(x,y)
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
targetReg={'PrL L5','BLA','vCA1'};

tableGapX=2;
tableWidth=100;

nCellHeight=10;

yGapIntraTop=10;

tableHight=nCellHeight*2+yGapIntraTop;

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
yShift=0;
xShift=0;
subplotInMM(x+xShift,y+yShift,tableWidth,tableHight,tableWidth)
set(gca,'YDir','reverse')
xlim([0,tableWidth])
ylim([0,tableHight])
hold on
xMargin=2;
yMargin=0;
cellWidth=(tableWidth-xMargin*2)/4.7
cellHight=(tableHight-yMargin*2)/4.5;
lineGap=0
rowIdx=0;
plot([0,tableWidth],(rowIdx+lineGap)*cellHight+yMargin+[0,0],'k-','LineWidth',0.5)
rowIdx=rowIdx+1;
for colIdx=1:4
    if colIdx<=length(targetReg)
    text((colIdx)*cellWidth+xMargin,(rowIdx-0.5)*cellHight+yMargin,{'Coupled ' ['with ' targetReg{colIdx}]},'fontsize',5,...
        'horizontalAlign','left','verticalALign','middle')
    else
    text((colIdx)*cellWidth+xMargin,(rowIdx-0.5)*cellHight+yMargin,'Other cells','fontsize',5,...
        'horizontalAlign','left','verticalALign','middle')
    end
end
plot([0,tableWidth],(rowIdx+lineGap)*cellHight+yMargin+[0,0],'k-','LineWidth',0.5)
for idx=1:length(targetReg)
    rowIdx=rowIdx+1;        
    target=find(strcmp(reg,targetReg{idx}));

    
    for colIdx=0:4;
        interpreter='tex';
        if colIdx==0
            txt={'Cells' ['in ' targetReg{idx}]};
        elseif colIdx==idx
            txt='N.A.';
        else
            if colIdx==4
                id=target(cellfun(@isempty,partner(target)));            
            else
                id=target(cellfun(@(x) any(strcmp(x,targetReg{colIdx})), partner(target)));
            end
            cnt=histcounts(cellType(id),-1.5:1.5);
            txt=sprintf('%d / %d / %d',cnt(3),cnt(1),cnt(2));
        end
        text((colIdx)*cellWidth+xMargin,(rowIdx-0.5)*cellHight+yMargin,txt,'fontsize',5,...
            'horizontalAlign','left','verticalALign','middle','Interpreter',interpreter)
    end
    
end
plot([0,tableWidth],(rowIdx+lineGap)*cellHight+yMargin+[0,0],'k-','LineWidth',0.5)
axis off

end





