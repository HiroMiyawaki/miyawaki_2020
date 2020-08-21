function coactPaper_figS15()
labelWeight=true;
labelCase=true;
labelSize=9;
letGapX=5;
letGapY=6;

close all
fh=initFig('height',3,'nColumn',1);

x=10;y=5;
panel_01(x,y);
drawnow();

print(fh,'figS15.pdf','-dpdf','-painters','-r300')
%
end
%%

function panel_01(x,y)

width=18-8/3;
wGap=4;
hGap=7;
height=9;
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
rate.fastGamma=[];

sigQ=[];
reg={};

for ratIdx=1:length(ratList)
    rat=ratList{ratIdx};
    target=find(qPETH.(rat).sigLevel==1);
    
    rate.SWR=cat(2,rate.SWR,qPETH.(rat).swr.avg(target,:)');
    rate.aHFO=cat(2,rate.aHFO,qPETH.(rat).hfo.avg(target,:)');
    rate.fastGamma=cat(2,rate.fastGamma,gPETH.(rat).peth(3).avg(target,:)');
    
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
tWin=qPETH.(ratList{1}).param.tWin;
nWin=ceil(tWin/tBinSize);
tBin=(-nWin:nWin)*tBinSize*1e3;
%%
smSigma=0.02;
smBin=0:tBinSize:smSigma*4;
smBin=[-fliplr(smBin),smBin(2:end)];
smCore=normpdf(smBin,0,smSigma);
smCore=smCore/sum(smCore);
%%
trigList={'fastGamma'};

trigName.SWR='SWR peak';
trigName.aHFO='aHFO peak';
trigName.slowGamma='PrL slow gamma peak';
trigName.fastGamma='PrL fast gamma peak';
trigName.cRipple='PrL ripple peak';

%%
pairList={'BLA' ,'PrL L5' ;
    'vCA1' 'PrL L5'};

sortBin=(size(rate.(trigList{1}),1)+1)/2+(-3:3);
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
    
    
    for trigIdx=1

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

subplotInMM(x,y+height+7,width*2+wGap,1)
imagescXY(cLim,[0,1],linspace(cLim(1),cLim(2),size(col.coact.map,1))')
set(gca,'clim',cLim)
colormap(gca,col.coact.map)
set(gca,'YTick',[])
xlim(cLim)
xlabel('Coactivation strength (z)','fontsize',5)
box off


end

