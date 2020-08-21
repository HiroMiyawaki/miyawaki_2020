function coactPaper_fig1()
labelWeight=true;
labelCase=true;
labelSize=9;
letGapX=5;
letGapY=6;
% 
close all
fh=initFig('height',9);

x=13;y=7;
panel_01(x,y);
panelLetter2(x-letGapX-7,y-letGapY,alphabet(1,labelCase),'fontSize',labelSize,'isBold',labelWeight)
drawnow();

x=60;y=7;
panel_02(x,y);
panelLetter2(x-letGapX,y-letGapY,alphabet(2,labelCase),'fontSize',labelSize,'isBold',labelWeight)
drawnow();

x=6;y=40;
panel_03(x,y);
panelLetter2(x-letGapX,y-letGapY,alphabet(3,labelCase),'fontSize',labelSize,'isBold',labelWeight)
drawnow();

x=56;y=40;
panel_04(x,y);
panelLetter2(x-letGapX,y-letGapY,alphabet(4,labelCase),'fontSize',labelSize,'isBold',labelWeight)
drawnow();

x=58;y=64;
panel_05(x,y);
panelLetter2(x-letGapX-2,y-letGapY,alphabet(5,labelCase),'fontSize',labelSize,'isBold',labelWeight)
drawnow();

x=89;y=64;
panel_06(x,y);
panelLetter2(x-letGapX-2,y-letGapY,alphabet(6,labelCase),'fontSize',labelSize,'isBold',labelWeight)
drawnow();


print(fh,'fig01.pdf','-dpdf','-painters','-r300')

end
%%
function panel_01(x,y)
    load('~/data/Fear/triple/hoegaarden181115/hoegaarden181115.basicMetaData.mat')
    load([basicMetaData.AnalysesName '-icaReacZNCCG.mat'])
    load([basicMetaData.AnalysesName '-icaReacCCG_sig.mat'])

    col=setCoactColor;
    
    totalHeigh=18;
    width=15;

    gapY=1;
    gapX=3.5;
    totalY=0;
    smSigma=20; %in ms
    cLim=0.01*[-1,1];

    reg=icaReacZNCCG(2).region(icaReacZNCCG(2).pairID);
    peakVal=icaReacCCG_sig(2).nrem.peakValue(:,[2,3]);
    sig=icaReacCCG_sig(2).nrem.significance(:,[2,3]);
    ccgVal=icaReacZNCCG(2).nrem.real.ccg(:,:,2:3);

    tBinSize=icaReacZNCCG(2).tBinSize*1e3;
    nShowBin=21;
    cBin=(size(ccgVal,2)+1)/2;
     
    tBin=(-nShowBin:nShowBin)*tBinSize;

    nSm=ceil(smSigma/tBinSize);
    xSm=(-nSm*4:nSm*4)*tBinSize;
    smCore=normpdf(xSm,0,smSigma);
    smCore=smCore/sum(smCore);

    for n=1:2
        ccgVal(:,:,n)=Filter0(smCore,ccgVal(:,:,n));
    end
    ccgVal=ccgVal(:,cBin+(-nShowBin:nShowBin),:);
    nPair=0;
   for n=1:3
        switch n
            case 1
                target={'BLA','PrL L5'};
            case 2
                target={'vCA1','PrL L5'};
            case 3
                target={'vCA1','BLA'};
            otherwise
        end
        nPair=nPair+sum(strcmp(reg(:,1),target{1})&strcmp(reg(:,2),target{2}));
   end
   
   eachHight=(totalHeigh-gapY*2)/nPair;
    for n=1:3
        switch n
            case 1
                target={'BLA','PrL L5'};
            case 2
                target={'vCA1','PrL L5'};
            case 3
                target={'vCA1','BLA'};
            otherwise
                continue
        end

        idx=find(strcmp(reg(:,1),target{1})&strcmp(reg(:,2),target{2}));
        subSig=sig(idx,:);
        subPeak=peakVal(idx,:);
        fprintf('%s-%s, n=%d\n',target{:},length(idx));

        [~,order]=sort(mean(ccgVal(idx,nShowBin+1+(-3:3),2),2),'descend');
        idx=idx(order);
        subSig=subSig(order,:);

        height=length(idx)*eachHight;% cm
        for m=0:1
            subplotInMM(x+(width+gapX)*m,y+totalY,width,height)
            imagesc(tBin,1:length(idx),ccgVal(idx,:,1+m))
            box off
            set(gca,'ytick',[])
            if n~=3
                set(gca,'xtick',[])
            else
                set(gca,'xtick',-400:400:400)
            end
            xlim(400*[-1,1])
            set(gca,'clim',cLim)
            ax=fixAxis;
            if n==1
                if m==0
                    title({'Pre-conditioning' 'NREM'},'fontweight','normal','fontsize',5)
                else
                    title({'Post-conditioning' 'NREM'},'fontweight','normal','fontsize',5)
                end
            end
            if n==3
                xlabel('\Deltatime (ms)','fontsize',5)
            end
            if m==0
               text2(-0.05,0.5,join(target, ' - '),ax,'fontsize',5,'horizontalAlign','right')
            end
            colormap(gca,col.coact.map)
        end
        subplotInMM(x+width+0.5,y+totalY,2.5,height)
        imagesc([subSig(:,1),-2*ones(size(subSig(:,1))),subSig(:,2)])
        set(gca,'clim',[-2,1])
        colormap(gca,[1,1,1;flipud(col.pVal)])
        box off
        axis off
        

        totalY=totalY+height+gapY;
    end

        subplotInMM(x+width*2+gapX+0.5,y,1,totalY-gapY)
        imagescXY([],cLim,linspace(cLim(1),cLim(2),512));
        set(gca,'clim',cLim)
        colormap(gca,col.coact.map)
        box off
        set(gca,'XTick',[])
        set(gca,'YAxisLocation','right')
        set(gca,'YTick',[cLim(1),0,cLim(2)])
        ax=fixAxis;
        text2(7,0.5,'Correlation',ax,'horizontalALign','center','Rotation',-90)
end

%%
function panel_02(x,y)
    load('~/data/Fear/triple/hoegaarden181115/hoegaarden181115.basicMetaData.mat')
    load([basicMetaData.AnalysesName '-icaReacZNCCG.mat'])
    load([basicMetaData.AnalysesName '-icaReacCCG_sig.mat'])
    width=54;
    height=18;
    subplotInMM(x,y,width,height,true)

    col=setCoactColor();
    
    preHC=[1,2,3,3,3,3,4];
    chColAcross=flipud(col.pVal);
    changeName={'','Kept','Gained'};
    legChange={};
    for n=1:length(changeName)
        legChange{n}=sprintf('\\color[rgb]{%f %f %f}%s',chColAcross(length(changeName)+1-n,:),changeName{length(changeName)+1-n});
    end
        

    sig=icaReacCCG_sig(2).nrem.significance;

    [regList,~,regIdx]=unique(icaReacCCG_sig(2).region);
    regIdx=reshape(regIdx,size(icaReacCCG_sig(2).region));
    nID=length(regIdx);
            
            
    labelPos=arrayfun(@(x) mean(find(regIdx==x)),1:length(regList));          

    across=find(diff(regIdx(icaReacCCG_sig(2).pairID),1,2)~=0);
                
    hold on

    sigLevel=sig(:,2:3)==1;
    
    typeName={'Pre-conditioning', 'Post-conditioning', 'Difference'};

    ensLeg='Ensembles in '
    for n=1:length(labelPos)
        nameID=strrep(regList{n},' ', '');
        nameID=strrep(nameID,'/', '');
        
        ensLeg= [ensLeg, sprintf('\\color[rgb]{%f %f %f}%s ',col.region.(nameID),regList{n})];
        if n<length(labelPos)
            ensLeg=[ensLeg,sprintf('\\color[rgb]{0 0 0}/ ')];
        end
    end
    
    for typeIdx=1:3
        centerX=3*(typeIdx-1);
                
        if typeIdx~=3
            for idx=across'
                if ~sigLevel(idx,typeIdx)
                    continue
                end
                plot(cos(2*pi/nID*icaReacCCG_sig(2).pairID(idx,:))+centerX,...
                    sin(2*pi/nID*icaReacCCG_sig(2).pairID(idx,:)),'-','color','k','LineWidth',0.5)

            end
        else
            for idx=across'
                if ~sigLevel(idx,1) && ~sigLevel(idx,2)
                    continue
                end
                plot(cos(2*pi/nID*icaReacCCG_sig(2).pairID(idx,:))+centerX,...
                    sin(2*pi/nID*icaReacCCG_sig(2).pairID(idx,:)),'-','color',chColAcross(diff(sigLevel(idx,:))+2,:),'LineWidth',0.5)

            end
        end

        for n=1:length(labelPos)
            nameID=strrep(regList{n},' ', '');
            nameID=strrep(nameID,'/', '');

            scatter(centerX+cos(2*pi/nID*find(regIdx==n)),sin(2*pi/nID*find(regIdx==n)),5,col.region.(nameID),'fill')

        end
        text(centerX,1.2,typeName{typeIdx},'horizontalAlign','center','verticalALign','baseline','fontsize',5)
        
        if typeIdx==3
            xlim([-1.1,1.1+centerX])
            ylim([-1.2,1.2])
            ax=fixAxis;
            text2(1,1.1,legChange,ax,'verticalALign','top','fontsize',5)
            text2(1,0,ensLeg,ax,'horizontalALign','right','verticalALign','top','fontsize',5)
        end
    end
    axis equal
    axis off
    
end

%%
function panel_03(x,y)

    width=15;
    totalHeigh=38;
    
    gapY=1;
    gapX=3.5;
    smSigma=20; %in ms
    cLim=0.01*[-1,1];
    nShowBin=21;
    ccg=poolVar('icaReacZNCCG.mat');
    ccgSig=poolVar('icaReacCCG_sig.mat');
    
    col=setCoactColor;
    
    ratList=fieldnames(ccg);

    reg={};
    peakVal=[];
    sig=[];
    ccgVal=[];
    for ratIdx=1:length(ratList)
        rat=ratList{ratIdx};
        reg=[reg;ccg.(rat)(2).region(ccg.(rat)(2).pairID)];
        peakVal=[peakVal;ccgSig.(rat)(2).nrem.peakValue(:,[2,3])];
        sig=[sig;ccgSig.(rat)(2).nrem.significance(:,[2,3])];
        ccgVal=cat(1,ccgVal,ccg.(rat)(2).nrem.real.ccg(:,:,2:3));
    end
    
    tBinSize=ccg.(ratList{1})(2).tBinSize*1e3;

    nSm=ceil(smSigma/tBinSize);
    xSm=(-nSm*4:nSm*4)*tBinSize;
    smCore=normpdf(xSm,0,smSigma);
    smCore=smCore/sum(smCore);
    
    cBin=(size(ccgVal,2)+1)/2;
     
    tBin=(-nShowBin:nShowBin)*tBinSize;
    
    for n=1:2
        ccgVal(:,:,n)=Filter0(smCore,ccgVal(:,:,n));
    end
    
    ccgVal=ccgVal(:,cBin+(-nShowBin:nShowBin),:);
    nPair=0;
   for n=1:3
        switch n
            case 1
                target={'BLA','PrL L5'};
            case 2
                target={'vCA1','PrL L5'};
            case 3
                target={'vCA1','BLA'};
            otherwise
        end
        nPair=nPair+sum(strcmp(reg(:,1),target{1})&strcmp(reg(:,2),target{2}));
   end
   
   eachHight=(totalHeigh-gapY*2)/nPair;
%%
    totalY=0;
   for n=1:3
        switch n
            case 1
                target={'BLA','PrL L5'};
            case 2
                target={'vCA1','PrL L5'};
            case 3
                target={'vCA1','BLA'};
            otherwise
                continue
        end

        idx=find(strcmp(reg(:,1),target{1})&strcmp(reg(:,2),target{2}));
        subSig=sig(idx,:);
        subPeak=peakVal(idx,:);

        [~,order]=sort(mean(ccgVal(idx,nShowBin+1+(-3:3),2),2),'descend');
        idx=idx(order);
        subSig=subSig(order,:);

        height=length(idx)*eachHight;% cm
        for m=0:1
            subplotInMM(x+(width+gapX)*m,y+totalY,width,height,true)
            imagesc(tBin,1:length(idx),ccgVal(idx,:,1+m))
            box off
            set(gca,'ytick',[])
            if n~=3
                set(gca,'xtick',[])
            else
                set(gca,'xtick',-400:400:400)
            end
            xlim(400*[-1,1])
            set(gca,'clim',cLim)
            colormap(gca,col.coact.map)

            if n==1
                if m==0
                    title({'Pre-conditioning' 'NREM'},'fontweight','normal','fontsize',5)
                else
                    title({'Post-conditioning' 'NREM'},'fontweight','normal','fontsize',5)
                end
            end
            if n==3
                xlabel('\Deltatime (ms)','fontsize',5)
            end
            if m==0
                ylabel(join(target, ' - '),'fontsize',5)
            end
        end
        subplotInMM(x+width+0.5,y+totalY,2.5,height)
        imagesc([subSig(:,1),-2*ones(size(subSig(:,1))),subSig(:,2)])
        set(gca,'clim',[-2,1])
        colormap(gca,[1,1,1;flipud(col.pVal)])
        box off
        axis off


        totalY=totalY+height+gapY;
   end

    subplotInMM(x+width*2+gapX+0.5,y,1,totalY-gapY)
    imagescXY([],cLim,linspace(cLim(1),cLim(2),512));
    set(gca,'clim',cLim)
    colormap(gca,col.coact.map)
    box off
    set(gca,'XTick',[])
    set(gca,'YAxisLocation','right')
    set(gca,'YTick',[cLim(1),cLim(1)/2,0,cLim(2)/2,cLim(2)])
    ax=fixAxis;
    text2(7,0.5,'Correlation',ax,'horizontalALign','center','Rotation',-90)
    
end

%%
function panel_04(x,y)
    coact=poolVar('icaReacZNCCG_sig.mat');
    tempIdx=2;
    beh='nrem';
    
    width=9;
    height=11;
    withinGap=0;
    acrossGap=4;
    
    sig=[];
    reg={};
    ratList=fieldnames(coact);
    for ratIdx=1:length(ratList)
        rat=ratList{ratIdx};
        
        sig=[sig;coact.(rat)(tempIdx).(beh).significance(:,2:3)];
        
        tempReg=relabel_region(coact.(rat)(tempIdx).region,'minCellNum',0);
        reg=[reg;tempReg(coact.(rat)(tempIdx).pairID)];
        
    end
    
    [pairList,~,pairIdx]=uniqueCellRows(reg);
    
    
    targetPair={'BLA','PrL L5'
                'vCA1','PrL L5'
                'vCA1','BLA'};
    targetPairIdx=[]
    for n=1:3
        targetPairIdx(n)=find(strcmp(pairList(:,1),targetPair{n,1})&strcmp(pairList(:,2),targetPair{n,2}));
    end
    temp=setCoactColor();
    col=flipud(temp.pVal);
    pFrac=[];
    frac=[];
    for n=1:length(targetPairIdx)
        subSig=sig(pairIdx==targetPairIdx(n),:);

        observed=[histcounts(subSig(:,1),-1.5:1.5);
            histcounts(subSig(:,2),-1.5:1.5)];
        frac(:,:,n)=observed;


        observed(:,sum(observed,1)==0)=[];

        if size(observed,2)<2
            pFrac(n)=1;
            continue
        end

        expect=sum(observed,2)*sum(observed,1)/sum(observed(:));

        chi2=sum((observed(:)-expect(:)).^2./expect(:));
        df=prod(size(observed)-1);
        pFrac(n)=chi2cdf(chi2,df,'upper');
    end    
    prePost={'Pre','Post'};
    for n=1:3
        subplotInMM(x+(2*width+withinGap+acrossGap)*(n-1),y,2*width+withinGap,height)
        
        sigFontSize=7;
        sigYshift=1.25;
        if pFrac(n)<0.001
            sigTxt='***';            
        elseif pFrac(n)<0.01
            sigTxt='**';
        elseif pFrac(n)<0.05
            sigTxt='*';
        else
            sigTxt='';
            sigFontSize=5;
            sigYshift=0;
        end
        if ~isempty(sigTxt)
            plot(width/2+[0,0,1,1]*(width+withinGap),height-2+[0,1,1,0]*0.5,'k-','LineWidth',0.5)
            text(width+withinGap/2,height-1.5-sigYshift,sigTxt,'HorizontalAlignment','center','VerticalAlignment','bottom','FontSize',sigFontSize)
        end
        
        title(join(pairList(targetPairIdx(n),:),' - '),'fontweight','normal','fontsize',5)
        xlim([0,width*2+withinGap])
        ylim([0,height])
        axis off
        for m=1:2
            subplotInMM(x+(width+withinGap)*(m-1)+(2*width+withinGap+acrossGap)*(n-1),y+2,width,height-2,[],true)
            h=pie(frac(m,:,n),{'','',''});
            for hIdx=1:length(h)
                if strcmpi(h(hIdx).Type,'patch')
                    h(hIdx).LineStyle='none'
                end
            end
            tempCol=col;
            tempCol(frac(1,:,n)<1,:)=[];
            colormap(gca,tempCol)
            ax=fixAxis;
            text2(0.5,0,prePost{m},ax,'verticalAlign','top','horizontalAlign','center')
        end
    end
    
    
    
end
%%
function panel_05(x,y)
    width=20;
    height=12;
    
    coact=poolVar('icaReacZNCCG_sig.mat');
    tempIdx=2;
    beh='nrem';
    
    peak=[];
    reg={};
    ratList=fieldnames(coact);
    for ratIdx=1:length(ratList)
        rat=ratList{ratIdx};
        
        peak=[peak;coact.(rat)(tempIdx).(beh).peakValue(:,2:3)];        
        
        tempReg=relabel_region(coact.(rat)(tempIdx).region,'minCellNum',0);
        reg=[reg;tempReg(coact.(rat)(tempIdx).pairID)];
        
    end
    
    [pairList,~,pairIdx]=uniqueCellRows(reg);
    
    
    targetPair={'BLA','PrL L5'
                'vCA1','PrL L5'
                'vCA1','BLA'};
    targetPairIdx=[];
    for n=1:3
        targetPairIdx(n)=find(strcmp(pairList(:,1),targetPair{n,1})&strcmp(pairList(:,2),targetPair{n,2}));
    end
    col=setCoactColor();
    peakDiffMean=[];
    peakDiffSTE=[];
    pDiff=[];
    for n=1:length(targetPairIdx)
        subPeak=peak(pairIdx==targetPairIdx(n),:);
        peakDiffMean(:,n)=mean(diff(subPeak,1,2));
        peakDiffSTE(:,n)=ste(diff(subPeak,1,2),[],1);
        if any(~isnan(diff(subPeak,1,2)))
            pDiff(n)=signrank(diff(subPeak,1,2));
        else
            pDiff(n)=1;
        end

    end    
    
    subplotInMM(x,y,width,height)
    hold on
    for n=1:3
        
        temp=strrep(strrep(targetPair(n,:),' ',''),'/','');
        pairName=[temp{:}];
        bar(n,peakDiffMean(n),0.8,'LineStyle','none','FaceColor',col.pair.(pairName))
        posErr=peakDiffMean(n)+(2*(peakDiffMean(n)>0)-1)*peakDiffSTE(n);
        plot(n+[0,0],[peakDiffMean(n),posErr],'-','color',col.pair.(pairName))
        
        
        sigFontSize=7;
        sigYshift=0.2e-3;
        if pDiff(n)<0.001
            sigTxt='***';            
        elseif pDiff(n)<0.01
            sigTxt='**';
        elseif pDiff(n)<0.05
            sigTxt='*';
        else
            sigTxt='';
            sigFontSize=5;
            sigYshift=1e-3;
        end
        if ~isempty(sigTxt)
            text(n,posErr+sigYshift,sigTxt,'FontSize',sigFontSize,'HorizontalAlignment','center')
        end
    end    
    set(gca,'xtick',1:3,'XTickLabel',join(targetPair, ' - '),'XTickLabelRotation',30)
    xlim([0,4])
    ylim([-0.004,0.008])
    set(gca,'YTick',-0.004:0.004:0.008)
    ylabel('\Deltapeak correlation')
end

function panel_06_old(x,y)
    width=20;
    height=12;
    subplotInMM(x,y,width,height)
    target=poolVar('icaCoactTimeHT.mat');
    ratList=fieldnames(target);

    tempBeh=2;
    targetHC=3;
    beh='nrem';
    
    reg.(beh)={};
    gap.(beh)=[];
    for ratIdx=1:length(ratList)
        rat=ratList{ratIdx};
        idx=find(cellfun(@(x,y) ~strcmpi(x,y),target.(rat)(tempBeh,targetHC).(beh).region(:,1),target.(rat)(tempBeh,targetHC).(beh).region(:,2)));
        
        
        gap.(beh)=[gap.(beh);target.(rat)(tempBeh,targetHC).(beh).tGap(idx)];
        reg.(beh)=[reg.(beh);target.(rat)(tempBeh,targetHC).(beh).region(idx,:)];
    end
    
    [pairList,~,pairIdx]=uniqueCellRows(reg.(beh));
    
    targetPair={'BLA','PrL L5'
                'vCA1','PrL L5'
                'vCA1','BLA'};
    targetPairIdx=[];
    for n=1:3
        targetPairIdx(n)=find(strcmp(pairList(:,1),targetPair{n,1})&strcmp(pairList(:,2),targetPair{n,2}));
    end
    col=setCoactColor;
    hold on
    plot([0,4],[0,0],'k-','LineWidth',0.5)
    legTxtTop{1}='';
    legTxtTop{2}='\color[rgb]{0,0,0}precedes';

    legTxtBottom{1}='';
    legTxtBottom{2}='\color[rgb]{0,0,0}precedes';
    for n=1:3
        subGap=gap.(beh)(pairIdx==targetPairIdx(n))*20;
        [yv,xv]=ksdensity(subGap);
        yv=yv/max(yv)*0.3;
        subMed=nanmedian(subGap);
        medHalfDen=interp1(xv,yv,subMed);

        temp=strrep(strrep(targetPair(n,:),' ',''),'/','');
        pairName=[temp{:}];

        fill(n+[yv,-fliplr(yv)],[xv,fliplr(xv)],col.pair.(pairName),'EdgeColor','none')
        plot(n+medHalfDen*[-1,1],subMed+[0,0],'-','Color',0.9*[1,1,1],'LineWidth',0.5)
        p=signrank(subGap);
        
        fprintf('%s - %s: mean %f +/- %f, median %f, p= %f\n',targetPair{n,:},nanmean(subGap),nanste(subGap),subMed,p)

        sigFontSize=7;
        sigYshift=10;
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
            text(n,max(xv)+sigYshift,sigTxt,'HorizontalAlignment','center','FontSize',sigFontSize)
        end        
        if ~isempty(legTxtTop{1})
            legTxtTop{1}=[legTxtTop{1} '\color[rgb]{0,0,0}/'];
        end
        legTxtTop{1}=[legTxtTop{1} sprintf('\\color[rgb]{%f %f %f}%s',col.pair.(pairName),targetPair{n,2})];

        if ~isempty(legTxtBottom{1})
            legTxtBottom{1}=[legTxtBottom{1} '\color[rgb]{0,0,0}/'];
        end
        legTxtBottom{1}=[legTxtBottom{1} sprintf('\\color[rgb]{%f %f %f}%s',col.pair.(pairName),targetPair{n,1})];        
    end
    xlim([0,4])    
    set(gca,'xtick',1:3,'XTickLabel',join(targetPair, ' - '),'XTickLabelRotation',30)
    xlim([0,4])
    ylim(150*[-1,1])
    ylabel('\Deltatime (ms)')    
    ax=fixAxis;
    for n=1:2
        text2(1.02,1+(2-n)*0.15,legTxtTop{n},ax);
        text2(1.02,0+(2-n)*0.15,legTxtBottom{n},ax);
    end
end


function panel_06(x,y)
    width=20;
    height=12;
    subplotInMM(x,y,width,height)
    target=poolVar('icaCoactTimeHT.mat');
    ratList=fieldnames(target);

    tempBeh=2;
    targetHC=3;
    beh='nrem';
    
    reg.(beh)={};
    gap.(beh)=[];
    for ratIdx=1:length(ratList)
        rat=ratList{ratIdx};
        idx=find(cellfun(@(x,y) ~strcmpi(x,y),target.(rat)(tempBeh,targetHC).(beh).region(:,1),target.(rat)(tempBeh,targetHC).(beh).region(:,2)));
        
        
        gap.(beh)=[gap.(beh);target.(rat)(tempBeh,targetHC).(beh).tGap(idx)];
        reg.(beh)=[reg.(beh);target.(rat)(tempBeh,targetHC).(beh).region(idx,:)];
    end
    
    [pairList,~,pairIdx]=uniqueCellRows(reg.(beh));
    
    targetPair={'BLA','PrL L5'
                'vCA1','PrL L5'
                'vCA1','BLA'};
    targetPairIdx=[];
    for n=1:3
        targetPairIdx(n)=find(strcmp(pairList(:,1),targetPair{n,1})&strcmp(pairList(:,2),targetPair{n,2}));
    end
    col=setCoactColor;
    hold on
    plot([0,4],[0,0],'k-','LineWidth',0.5)
    legTxtTop{4}='\color[rgb]{0,0,0}precedes';

    legTxtBottom{1}='';
    legTxtBottom{4}='\color[rgb]{0,0,0}precedes';
    for n=1:3
        subGap=gap.(beh)(pairIdx==targetPairIdx(n))*20;
        [yv,xv]=ksdensity(subGap);
        yv=yv/max(yv)*0.3;
        subMed=nanmedian(subGap);
        medHalfDen=interp1(xv,yv,subMed);

        temp=strrep(strrep(targetPair(n,:),' ',''),'/','');
        pairName=[temp{:}];

        fill(n+[yv,-fliplr(yv)],[xv,fliplr(xv)],col.pair.(pairName),'EdgeColor','none')
        plot(n+medHalfDen*[-1,1],subMed+[0,0],'-','Color',0.9*[1,1,1],'LineWidth',0.5)
        p=signrank(subGap);
        fprintf('%s-%s: mean = %f +/- %f ms,median=%f ms\n',targetPair{n,:},nanmean(subGap),nanste(subGap),nanmedian(subGap))
        sigFontSize=7;
        sigYshift=10;
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
            text(n,max(xv)+sigYshift,sigTxt,'HorizontalAlignment','center','FontSize',sigFontSize)
        end
        
        legTxtTop{n}=sprintf('\\color[rgb]{%f %f %f}%s',col.pair.(pairName),targetPair{n,2});


        legTxtBottom{n}=sprintf('\\color[rgb]{%f %f %f}%s',col.pair.(pairName),targetPair{n,1});        
    end
    xlim([0,4])    
    set(gca,'xtick',1:3,'XTickLabel',join(targetPair, ' - '),'XTickLabelRotation',30)
    xlim([0,4])
    ylim(150*[-1,1])
    ylabel('\Deltatime (ms)')    
    ax=fixAxis;
    for n=1:4
        text2(1.1,1+(2.5-n)*0.15,legTxtTop{n},ax);
        text2(1.1,0+(2.5-n)*0.15,legTxtBottom{n},ax);
    end
    text2(1.075,1,'\color[rgb]{0,0,0}\uparrow',ax,'horizontalALign','right','verticalAlign','middle','fontsize',7);
    text2(1.075,0,'\color[rgb]{0,0,0}\downarrow',ax,'horizontalALign','right','verticalAlign','middle','fontsize',7);
    
end



%%


