function coactPaper_fig4()
labelWeight=true;
labelCase=true;
labelSize=9;
letGapX=5;
letGapY=4;

close all
fh=initFig4('height',7.5);


x=15;y=6;
panel_01(x,y);
panelLetter2(x-letGapX-2,y-letGapY,alphabet(1,labelCase),'fontSize',labelSize,'isBold',labelWeight)
drawnow();

x=15+54;y=6;
panel_02(x,y);
panelLetter2(x-letGapX-2,y-letGapY,alphabet(2,labelCase),'fontSize',labelSize,'isBold',labelWeight)
drawnow();


print(fh,'fig04.pdf','-dpdf','-painters','-r300')
end



function panel_01(x,y)
base=poolVar('basicMetaData.mat','base');
wave=poolVar('icaCoactTrigWaveletCueRet.mat');
sigCue=poolVar('icaReacZNCCGchamberCue_sig.mat');
sig=poolVar('icaReacZNCCG_sig.mat');

ratList=fieldnames(wave);

%%
tRange=420*[-1,1];
width=19;
hight=14;
xMargin=5;
yMargin=5;
cLim=2*[-1,1];
%%
tBinWavelet=(-(size(wave.(ratList{1}).wavelet,2)-1)/2:(size(wave.(ratList{1}).wavelet,2)-1)/2)/1.25;

reg={};
wavelet=[];
sigLebel=[];
sigQ=[];
pID=[];
animal=[];
col=setCoactColor();
for ratIdx=1:length(ratList)
    rat=ratList{ratIdx};
    reg=[reg;wave.(rat).region];
    animal=[animal,ratIdx*ones(1,size(wave.(rat).wavelet,3))];
    
    ch=wave.(rat).param.Ch;
    chNameList.(rat)=relabel_region(base.(rat).Ch.names,'minCellNum',0);
    
    wavelet=cat(3,wavelet,wave.(rat).wavelet);
    
    pID=[pID,wave.(rat).pairID'];
    sigLebel=[sigLebel;wave.(rat).sigLevel];

    temp=sig.(rat)(2).region(sig.(rat)(2).pairID);
    
    across=find(cellfun(@(x,y) ~strcmp(x,y),temp(:,1),temp(:,2)) & ...
        (sig.(rat)(2).nrem.significance5(:,3)>0)); 
    
    sigQ=[sigQ;sigCue.(rat)(2).significance(across)];
    
end
%%
fratio=mean(-diff(log2(wave.(ratList{1}).f)));
fMax=max(wave.(ratList{1}).f);

yTickLabel=[1,3,10,30,100,300];
yTick=length(wave.(ratList{1}).f)-(log2(fMax)-log2(yTickLabel))/fratio;


%%
[regList,~,regID]=unique(reg);
regID=reshape(regID,size(reg));

[pairList,~,pairID]=unique(regID,'rows');

cnt=histcounts(pairID,0.5:length(pairList)+0.5);
[~,order]=sort(cnt,'descend');

for typeIdx=1:2%length(order)
    target=find(pairID==order(typeIdx));
    target(sigLebel(target)~=1)=[];
    
    
    target=target(sigQ(target)==1);
    
    if isempty(target)
        continue
    end
    
    nCol=3;
    
    
    
    probeOrder=[1,3,2];
    for n=1:3
        temp={};
        for idx=1:length(target)
            rat=ratList{animal(target(idx))};
            temp{idx}=chNameList.(rat){wave.(rat).param.Ch(probeOrder(n))};
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
        plot([0,0],ax(3:4),'-','color',0.7*[1,1,1],'linewidth',0.25)
        hold on
        if typeIdx==1
            ylabel([probeName{n}{:} ' (Hz)'],'fontsize',5,'fontweight','normal')
        end
        box off
        if n==2 ||  (typeIdx==2 && n==1) || n==3
            fPos=get(gcf,'paperPosition');
            fPos=fPos*10;
            scale=fPos(3:4);
            if n==3
                xMM=x+(width+xMargin)*(mod(typeIdx-1,nCol))+width/2+[-1,0]*1.75-1;
                yMM=fliplr(y+(hight+yMargin)*(n-1)+[0,1]*1.75+3);
            else
                xMM=x+(width+xMargin)*(mod(typeIdx-1,nCol))+width/2+[-1,0]*1.75-1;
                yMM=fliplr(y+(hight+yMargin)*(n-1)+[0,1]*1.75+2.5);
            end
                
            if typeIdx==2
                lType='-';
            else
                lType='-';
            end
            annotation('arrow',xMM/scale(1) ,1-yMM/scale(2),'color','k','HeadWidth',4,'HeadLength',4,...
                'LineWidth',1,'LineStyle',lType)
            
            if typeIdx==2&&n==3
            xMM=x+(width+xMargin)*(mod(typeIdx-1,nCol))+width/2+[1,0]*1.75+1.25;
            yMM=fliplr(y+(hight+yMargin)*(n-1)+[0,1]*1.75+5.25);
                lType='none';
            annotation('arrow',xMM/scale(1) ,1-yMM/scale(2),'color','k','HeadWidth',4,'HeadLength',4,...
                'LineWidth',1,'LineStyle',lType)
            end
            
        end
        if n==1
            title(join(reg(target(1),:),' - '),'fontsize',5,'fontweight','normal');
        end
        if n==3 && typeIdx==1
            textInMM(x+width+xMargin/2,y+hight*3+yMargin*2+5,'Time from coactivation peak (ms)',...
                'horizontalALign','center','verticalALign','baseline')
        end
        
    end
    
end

subplotInMM(x,y+hight+(hight+yMargin)*2+7,width*2+xMargin,1)
imagescXY(cLim,[0,1],linspace(cLim(1),cLim(2),size(col.coact.map,1))')
colormap(col.wavelet.map)
xlim(cLim)
box off
set(gca,'YTick',[],'YAxisLocation','right')
xlabel('Power (z)','fontsize',5)

end


function panel_02(x,y)

width=19;
wGap=4;
hGap=7;
totalH=14*3+5*2;
height=(totalH-hGap*3)/4;
%%
qPETH=poolVar('evtTrigIcaCoact-cueRet.mat');
gPETH=poolVar('gammaTrigIcaCoact-cueRet.mat');

sigCue=poolVar('icaReacZNCCGchamberCue_sig.mat');

cLim=[-1,7];

ratList=fieldnames(qPETH);
yTick={0:3:15,0:2:10}
%%
rate.SWR=[];
rate.aHFO=[];
rate.slowGamma=[];
rate.cRipple=[];

sigQ=[];
reg={};

for ratIdx=1:length(ratList)
    rat=ratList{ratIdx};
    target=find(qPETH.(rat).sigLevel==1);
    
    rate.SWR=cat(2,rate.SWR,qPETH.(rat).swr.avg(target,:)');
    rate.aHFO=cat(2,rate.aHFO,qPETH.(rat).hfo.avg(target,:)');
    rate.slowGamma=cat(2,rate.slowGamma,gPETH.(rat).peth(2).avg(target,:)');
    rate.cRipple=cat(2,rate.cRipple,gPETH.(rat).peth(6).avg(target,:)');
    
    gPETH.(rat).peth(2).param.freqRange
    gPETH.(rat).peth(6).param.freqRange
    
    
    reg=cat(1,reg,qPETH.(rat).region(target,:));
    
    reacID=qPETH.(rat).reacID(target,:);
    for rIdx=1:size(reacID,1)
        idx=find((sigCue.(rat)(2).pairID(:,1)==reacID(rIdx,1) & sigCue.(rat)(2).pairID(:,2)==reacID(rIdx,2)) | ...
            (sigCue.(rat)(2).pairID(:,1)==reacID(rIdx,2) & sigCue.(rat)(2).pairID(:,2)==reacID(rIdx,1)));
        sigQ(end+1)=sigCue.(rat)(2).significance(idx);
        
    end
    
end


%%
tBinSize=0.02;
tWin=qPETH.(ratList{1}).param.tWin
nWin=ceil(tWin/tBinSize);
tBin=(-nWin:nWin)*tBinSize*1e3;
%%
smSigma=0.02;
smBin=0:tBinSize:smSigma*4;
smBin=[-fliplr(smBin),smBin(2:end)];
smCore=normpdf(smBin,0,smSigma);
smCore=smCore/sum(smCore);
%%
trigList={'SWR','aHFO','slowGamma','cRipple'};

trigName.SWR='SWR peak';
trigName.aHFO='aHFO peak';
trigName.slowGamma='PrL slow gamma peak';
trigName.cRipple='PrL cortical ripple peak';

%%
pairList={'BLA' ,'PrL L5' ;
    'vCA1' 'PrL L5'};

sortBin=(size(rate.aHFO,1)+1)/2+(-3:3);
col=setCoactColor;
for pairIdx=1:2
    target=find((strcmp(reg(:,1),pairList{pairIdx,1})&strcmp(reg(:,2),pairList{pairIdx,2}))|...
        (strcmp(reg(:,1),pairList{pairIdx,2})&strcmp(reg(:,2),pairList{pairIdx,1})));
    target=target(sigQ(target)==1);

    
    if pairIdx==1
        sortSignal=rate.aHFO(:,target);
    else
        sortSignal=rate.SWR(:,target);
    end
    sortSignal=Filter0(smCore,sortSignal)';
    sortSignal=zscore(sortSignal,[],2);
    [~,order]=sort(mean(sortSignal(:,sortBin),2),'descend');
    
    
    for trigIdx=1:4
        
        subplotInMM(x+(width+wGap)*(pairIdx-1),y+(hGap+height)*(trigIdx-1),width,height)

        signal=rate.(trigList{trigIdx})(:,target);
        signal=signal(:,order);
        
        signal=Filter0(smCore,signal)';
        signal=zscore(signal,[],2);
        imagesc(tBin,[],signal)
        xlim(440*[-1,1])        
        set(gca,'clim',[-1,7])
       set(gca,'xtick',-400:400:400)
       if pairIdx==1
           textInMM(x+width+wGap/2,y+(hGap+height)*(trigIdx-1)+height+5,...
               sprintf('Time from %s (ms)',trigName.(trigList{trigIdx})),...
               'horizontalAlign','center','verticalAlign','baseline')
       end
       
       if pairIdx==1 && trigIdx==4
           mid=(size(signal,2)+1)/2;
           [~,mxPos]=max(signal(:,mid+(-5:5)),[],2);
           gap=(mxPos-5)*tBinSize*1e3;
           mean(gap)
           ste(gap)
           signrank(gap)
           fprintf('%s triggered %s-%s coactivation; \\deltat = %f +/- %f, p=%f\n',...
                        trigList{trigIdx}, pairList{pairIdx,:},mean(gap), ste(gap),signrank(gap))
           
       end

        if pairIdx==1
            ylabel('# pair')
        end
        box off
        
        colormap(gca,col.coact.map)
        if trigIdx==1
            title(join(pairList(pairIdx,:),' - '),'fontsize',5,'fontweight','normal')
        end
    end
end
subplotInMM(x,y+totalH+7,width*2+wGap,1)
imagescXY(cLim,[0,1],linspace(cLim(1),cLim(2),size(col.coact.map,1))')
set(gca,'clim',cLim)
colormap(gca,col.coact.map)
set(gca,'YTick',[])
xlim(cLim)
xlabel('Coactivation strength (z)','fontsize',5)
box off


end


