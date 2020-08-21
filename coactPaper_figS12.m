function coactPaper_figS12()
labelWeight=true;
labelCase=true;
labelSize=9;
letGapX=3;
letGapY=6;

close all
fh=initFig('height',3,'nColumn',2);

x=15+27;y=5;
panel_05(x,y);
drawnow();

print(fh,'figS12.pdf','-dpdf','-painters','-r300')

end

%%
function panel_05(x,y)
width=18;
wGap=15;
height=12;

hfoDrop=poolVar('icaReacCCG_dropHFO_baseCond_sh.mat');
sig=poolVar('icaReacZNCCG_sig.mat');
hfoSig=poolVar('icaReacZNCCG_exHFObaseCond_sig.mat');

swrDrop=poolVar('icaReacCCG_dropSWR_sh.mat');
swrSig=poolVar('icaReacCCG_exSWR_sig.mat');


ratList=fieldnames(swrDrop);
%%
tempIdx=2; %conditioning
targetIdx=3; %homecage3: after conditioning
dropTypeList={'exRipple','exHFO','exSpindle'};


isSigAll=[];
regAll={};
for dIdx=1:length(dropTypeList)
    dropType=dropTypeList{dIdx};
    pAll.(dropType)=[];
    isExSigAll.(dropType)=[];
end

for ratIdx=1:length(ratList)
    rat=ratList{ratIdx};
    
    for dIdx=1:length(dropTypeList)
        dropType=dropTypeList{dIdx};
        if strcmp(dropType,'exHFO')
            shPeak=hfoDrop.(rat)(tempIdx).(dropType).shuffle.peak(:,:,targetIdx);
            realPeak=max(hfoDrop.(rat)(tempIdx).(dropType).real.ccg(:,:,targetIdx),[],2);
            isExSigAll.(dropType)=[isExSigAll.(dropType);hfoSig.(rat)(tempIdx).(dropType).significance(:,targetIdx)>0];
        else
            shPeak=swrDrop.(rat)(tempIdx).(dropType).shuffle.peak(:,:,targetIdx);
            realPeak=max(swrDrop.(rat)(tempIdx).(dropType).real.ccg(:,:,targetIdx),[],2);
            isExSigAll.(dropType)=[isExSigAll.(dropType);swrSig.(rat)(tempIdx).(dropType).significance(:,targetIdx)>0];
        end
        
        temp=zeros(size(realPeak));
        for n=1:size(shPeak,1)
            temp(n)=sum(shPeak(n,:)<realPeak(n))/500;
        end
        pAll.(dropType)=[pAll.(dropType);temp];
    end
    isSigAll=[isSigAll;sig.(rat)(tempIdx).nrem.significance(:,targetIdx)>0];
    regAll=[regAll;hfoDrop.(rat)(tempIdx).region(hfoDrop.(rat)(tempIdx).pairID)];
end

[pairListAll,~,pairIDall]=uniqueCellRows(regAll);
interRegAll=ismember(pairIDall,find(~strcmp(pairListAll(:,1),pairListAll(:,2))));

isSig=isSigAll(interRegAll);
reg=regAll(interRegAll,:);
for dIdx=1:length(dropTypeList)
    dropType=dropTypeList{dIdx};
    p.(dropType)=pAll.(dropType)(interRegAll);
    isExSig.(dropType)=isExSigAll.(dropType)(interRegAll);
    
end

[pairList,~,pairID]=uniqueCellRows(reg);
%%
dropName={'SWR','HFO','Spindle'};
targetPairList={'vCA1' 'BLA' };

temp=setCoactColor;
col=temp.pair.vCA1BLA;

legTxt={};
for n=1:size(targetPairList,1)
    legTxt{n}=sprintf('\\color[rgb]{%f %f %f}%s - %s',col(n,:), targetPairList{n,:})
end

for pairIdx=1:size(targetPairList,1)
    targePairtIdx=find(strcmp(pairList(:,1),targetPairList{pairIdx,1})& strcmp(pairList(:,2),targetPairList{pairIdx,2}));
    tarSig=isSig(pairID==targePairtIdx,:)
    for dIdx=1:length(dropTypeList)
        dropType=dropTypeList{dIdx};
        xx=isExSig.(dropType)(pairID==targePairtIdx)==1;
        yy=tarSig==1;
        ct=[sum(xx&yy)  sum(xx&~yy)
            sum(~xx&yy) sum(~xx&~yy)];
        
        frac=(ct(1,:)./sum(ct,1)*100)
        exEvent(pairIdx,dIdx)=100-frac(1);
        
        tarP=p.(dropType)(pairID==targePairtIdx);
        xx=tarP<(0.01/2);
        yy=tarSig==1;
        ct=[sum(xx&yy)  sum(xx&~yy)
            sum(~xx&yy) sum(~xx&~yy)];
        
        num=sum(ct,1);
        frac=(ct(1,:)./sum(ct,1)*100);
        jitEvent(pairIdx,dIdx)=frac(1);
        
    end
end

for n=1:2
        subplotInMM(x+(width+wGap)*(n-1),y,width,height)
    hold on
    if n==1
        val=jitEvent;
        yTxt={'Pair with' 'sig. drop (%)'};
        xTxt='Excluded Events';
    else
        val=exEvent;
        yTxt={'Pair lost' 'sig. peak (%)'};
        xTxt='Excluded Events';
    end
    for m=1
        bar((1:3)+0*(2*m-3),val(m,:),0.2,'FaceColor',col(m,:),'linestyle','none')
    end
    set(gca,'xtick',1:3,'XTickLabel',dropName)
    xlim([0.5,3.5])
    ylim([0,100])
    ylabel(yTxt,'fontsize',5,'fontweight','normal')
    xlabel(xTxt,'fontsize',5,'fontweight','normal')
    ax=fixAxis;
    text2(0.5,1.15,legTxt{1},ax,'fontsize',5)
end



end