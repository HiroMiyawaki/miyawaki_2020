function coactPaper_figS14()
labelWeight=true;
labelCase=true;
labelSize=9;
letGapX=5;
letGapY=6;

close all
fh=initFig('height',3);



x=15+23;y=5;
panel_01(x,y);
drawnow();

print(fh,'figS14.pdf','-dpdf','-painters','-r300')
%
end

%%
function panel_01(x,y)
width=25;
yGapInter=12;
height=12;
sigCue=poolVar('icaReacZNCCGchamberCue_sig.mat');
sigChamber=poolVar('icaReacZNCCGchamber_sig.mat');
sigHomecage=poolVar('icaReacZNCCG_sig.mat');
xGap=6;

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


targetPair={'BLA','PrL L2/3';
    'vCA1','BLA'};

colTempate=setCoactColor();
col=[
    colTempate.pair.BLAPrLL23;
    0.5*[1,1,1];
    colTempate.pair.vCA1BLA;
    0.5*[1,1,1];
    ];
yRange=[0,40;
    0,60];
for pairIdx=1:size(targetPair,1)
    subplotInMM(x+(width+xGap)*(pairIdx-1),y,width,height)
    
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
    if pairIdx==1
        ylabel({'Fraction of pairs with' 'significant peak (%)'})
    end
end
end
