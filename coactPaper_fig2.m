function coactPaper_fig2()
labelWeight=true;
labelCase=true;
labelSize=9;
letGapX=5;
letGapY=6;

close all
fh=initFig('height',7);

x=8;y=6;
panel_01(x,y);
panelLetter2(x-letGapX-1,y-letGapY+1,alphabet(1,labelCase),'fontSize',labelSize,'isBold',labelWeight)
drawnow();

x=66;y=6;
panel_02(x,y);
panelLetter2(x-letGapX-4,y-letGapY+1,alphabet(2,labelCase),'fontSize',labelSize,'isBold',labelWeight)
drawnow();

x=66;y=48.5+2.8;
panel_03(x,y)
panelLetter2(x-letGapX-4,y-letGapY+2,alphabet(3,labelCase),'fontSize',labelSize,'isBold',labelWeight)
panelLetter2(x-letGapX-4+30,y-letGapY+2,alphabet(4,labelCase),'fontSize',labelSize,'isBold',labelWeight)
drawnow();


print(fh,'fig02.pdf','-dpdf','-painters','-r300')

end
%%
function panel_01(x,y)
icaWavelet=poolVar('icaCoactTrigWaveletHT.mat');
base=poolVar('basicMetaData.mat','base');

ratList=fieldnames(icaWavelet);

%%
tRange=420*[-1,1];
width=16;
hight=14;
xMargin=8;
yMargin=5;
cLim=0.3*[-1,1];
%%
tBinWavelet=(-(size(icaWavelet.(ratList{1}).nrem.wavelet,2)-1)/2:(size(icaWavelet.(ratList{1}).nrem.wavelet,2)-1)/2)/1.25;
reactName='ICA reactivation';

reg={};
wavelet=[];
sigLebel=[];
pID=[];
animal=[];
col=setCoactColor();
for ratIdx=1:length(ratList)
    rat=ratList{ratIdx};
    reg=[reg;icaWavelet.(rat).nrem.region];
    animal=[animal,ratIdx*ones(1,size(icaWavelet.(rat).nrem.wavelet,3))];
    
    ch=icaWavelet.(rat).param.Ch;
    chNameList.(rat)=relabel_region(base.(rat).Ch.names,'minCellNum',0);
    
    wavelet=cat(3,wavelet,icaWavelet.(rat).nrem.wavelet);
    
    pID=[pID,icaWavelet.(rat).nrem.pairID'];
    sigLebel=[sigLebel;icaWavelet.(rat).nrem.sigLevel];
end
%%
fratio=mean(-diff(log2(icaWavelet.(ratList{1}).f)));
fMax=max(icaWavelet.(ratList{1}).f);

yTickLabel=[1,3,10,30,100,300];
yTick=length(icaWavelet.(ratList{1}).f)-(log2(fMax)-log2(yTickLabel))/fratio;

swrFreqPos=length(icaWavelet.(ratList{1}).f)-(log2(fMax)-log2([110,130,150,170]))/fratio;


debug=false;

%%
[regList,~,regID]=unique(reg);
regID=reshape(regID,size(reg));

[pairList,~,pairID]=unique(regID,'rows');

cnt=histcounts(pairID,0.5:length(pairList)+0.5);
[~,order]=sort(cnt,'descend');

for typeIdx=1:2%length(order)
    target=find(pairID==order(typeIdx));
    target(sigLebel(target)~=1)=[];
    
    if isempty(target)
        continue
    end
    
    nCol=3;
    
    
    
    probeOrder=[1,3,2];
    for n=1:3
        temp={};
        for idx=1:length(target)
            rat=ratList{animal(target(idx))};
            temp{idx}=chNameList.(rat){icaWavelet.(rat).param.Ch(probeOrder(n))};
        end
        probeName{n}=join(unique(temp),'/');
    end
    
    for n=1:3
        
        subplotInMM(x+(width+xMargin)*(mod(typeIdx-1,nCol)),y+(hight+yMargin)*(n-1),width,hight)
        hold on
        imagesc(tBinWavelet,[],mean(flipud(wavelet(:,:,target,probeOrder(n))),3))
        set(gca,'YTick',yTick,'YTickLabel',yTickLabel)
        set(gca,'xtick',-400:400:400)
        colormap(gca,col.wavelet.map)
        axis tight

        set(gca,'CLim',cLim)
        xlim(tRange)
        ax=fixAxis;
        hold on
        plot([0,0],ax(3:4),'-','color',0.7*[1,1,1],'linewidth',0.25)
        if debug
            for ii=1:length(swrFreqPos)
                plot(ax(1:2),swrFreqPos(ii)+[0,0],'k-','linewidth',0.25)
            end
        end
        ylabel([probeName{n}{:} ' (Hz)'],'fontsize',5,'fontweight','normal')

        box off
        if (typeIdx==1 && n==2) ||  (typeIdx==2 && n==1)
            fPos=get(gcf,'paperPosition');
            fPos=fPos*10;
            scale=fPos(3:4);

            xMM=x+(width+xMargin)*(mod(typeIdx-1,nCol))+width/2+[-1,0]*1.75-1;
            yMM=fliplr(y+(hight+yMargin)*(n-1)+[0,1]*1.75+2.5);
            if typeIdx==2
                lType='none';
            else
                lType='-';
            end
            annotation('arrow',xMM/scale(1) ,1-yMM/scale(2),'color','w','HeadWidth',4,'HeadLength',4,...
                'LineWidth',1,'LineStyle',lType)
        end
        if n==1
            title(join(reg(target(1),:),' - '),'fontsize',5,'fontweight','normal');
        end
        if n==3
            xlabel({'Time from ' 'coactivation peak (ms)'},'fontsize',5,'fontweight','normal')
        end
        
    end
    
end
subplotInMM(x+(width+xMargin)*2-xMargin+1.5,y,1,hight+(hight+yMargin)*2)
imagescXY([0,1],cLim,linspace(cLim(1),cLim(2),size(col.coact.map,1)))
colormap(col.wavelet.map)
xlim([0,1])
box off
set(gca,'XTick',[],'YAxisLocation','right')
ax=fixAxis;
text2(6.5,0.5,'Power (z)',ax,'fontsize',5,'rotation',-90,'horizontalALign','center')


end

%%
function panel_02(x,y)

width=18;
wGap=9;
hGap=8;
height=(33-hGap*2)/3;

nremPETH=poolVar('evtTrigIcaCoact-postNREM.mat');

cLim=[-1,8];
xLim.SWR=420*[-1,1];
xLim.HFO=420*[-1,1];
xLim.spindle=2*[-1,1];

ratList=fieldnames(nremPETH);
yTick={0:15:40,0:5:10}
%%
rate.SWR=[];
rate.spindle=[];
rate.HFO=[];
sig=[];
reg={};

for ratIdx=1:length(ratList)
    rat=ratList{ratIdx};
    target=find(nremPETH.(rat).sigLevel==1);

    rate.SWR=cat(2,rate.SWR,nremPETH.(rat).swr.avg(target,:)');
    rate.HFO=cat(2,rate.HFO,nremPETH.(rat).hfo.avg(target,:)');
    rate.spindle=cat(2,rate.spindle,nremPETH.(rat).spindle.avg(target,:)');

    reg=cat(1,reg,nremPETH.(rat).region(target,:));
end
tBinSize=0.02;
tWin=nremPETH.(ratList{1}).param.tWin
nWin=ceil(tWin/tBinSize);
tBin=(-nWin:nWin)*tBinSize*1e3;
%%
smSigma=0.02;
smBin=0:tBinSize:smSigma*4;
smBin=[-fliplr(smBin),smBin(2:end)];
smCore=normpdf(smBin,0,smSigma);
smCore=smCore/sum(smCore);
%%
trigList={'SWR','HFO','spindle'};
for trigIdx=1:3
    rate.(trigList{trigIdx})=zscore(Filter0(smCore,squeeze(rate.(trigList{trigIdx}))),[],1);
end
%%
pairList={'BLA' ,'PrL L5' ;
    'vCA1' 'PrL L5'};

sortBin=(size(rate.HFO,1)+1)/2+(-1:1);
col=setCoactColor;
for pairIdx=1:2
    target=find((strcmp(reg(:,1),pairList{pairIdx,1})&strcmp(reg(:,2),pairList{pairIdx,2}))|...
        (strcmp(reg(:,1),pairList{pairIdx,2})&strcmp(reg(:,2),pairList{pairIdx,1})));
    
    if pairIdx==1
        sortSignal=rate.HFO(:,target);
    else
        sortSignal=rate.SWR(:,target);
    end
    [~,order]=sort(mean(sortSignal(sortBin,:),1));
    
    
    for trigIdx=1:3
        subplotInMM(x+(wGap+width)*(pairIdx-1),y+(hGap+height)*(trigIdx-1),width,height)
        signal=rate.(trigList{trigIdx})(:,target);
        signal=signal(:,order);
        if strcmp((trigList{trigIdx}),'spindle')
            imagescXY(tBin/1e3,[],signal)
            xlabel(sprintf('Time from %s peak (s)',trigList{trigIdx}))
            set(gca,'xtick',-2:2:2)
        else
            imagescXY(tBin,[],signal)
            xlabel(sprintf('Time from %s peak (ms)',trigList{trigIdx}))
            set(gca,'xtick',-400:400:400)
        end
        set(gca,'YTick',yTick{pairIdx});
        ylabel('# pair')
        box off
        colormap(gca,col.coact.map)
        set(gca,'clim',cLim)
        ax=fixAxis;

        xlim(xLim.(trigList{trigIdx}))
        if trigIdx==1
            title(join(pairList(pairIdx,:),' - '),'fontsize',5,'fontweight','normal')
        end
 
        
    end
end
subplotInMM(x+(wGap+width)*2-wGap+2,y,1,(hGap+height)*3-hGap)
imagescXY([0,1],cLim,linspace(cLim(1),cLim(2),size(col.coact.map,1)))
set(gca,'clim',cLim)
colormap(gca,col.coact.map)
set(gca,'XTick',[],'YAxisLocation','right')
xlim([0,1])
text(4.5,mean(cLim),'Coactivation strength (z)','horizontalAlign','center','rotation',-90,'fontsize',5)
box off

end
%%

function panel_03(x,y)
width=16;
wGap=15.5;
height=9;

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
targetPairList={'BLA' ,'PrL L5' ;
    'vCA1' 'PrL L5'};

temp=setCoactColor;
col=[temp.pair.BLAPrLL5;
    temp.pair.vCA1PrLL5];

legTxt={};
for n=1:2
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
        xTxt='Excluded events';
    else
        val=exEvent;
        yTxt={'Pair lost' 'sig. peak (%)'};
        xTxt='Excluded events';
    end
    for m=1:2
        bar((1:3)+0.15*(2*m-3),val(m,:),0.2,'FaceColor',col(m,:),'linestyle','none')
    end
    set(gca,'xtick',1:3,'XTickLabel',dropName)
    xlim([0.5,3.5])
    ylim([0,100])
    ylabel(yTxt,'fontsize',5,'fontweight','normal')
    xlabel(xTxt,'fontsize',5,'fontweight','normal')
    ax=fixAxis;
    text2(0.5,1.15,legTxt{1},ax,'fontsize',5)
    text2(0.5,1,legTxt{2},ax,'fontsize',5)
end



end

