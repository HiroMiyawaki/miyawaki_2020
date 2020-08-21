function coactPaper_figS19()
labelWeight=true;
labelCase=true;
labelSize=9;
letGapX=9;
letGapY=3;

close all
fh=initFig('height',9);

x=16;y=4;
panel_01(x,y);
drawnow();


print(fh,figS19.pdf','-dpdf','-painters','-r300')

end
%%

%%
function panel_01(x,y)
yGapIntraBottom=14;

frzModwidth=30;
frzModheigth=15;
frzModGapXIntra=24;

coact=poolVar('coactCompCell.mat');
info=poolVar('okUnit.cellinfo.mat');
frz=poolVar('frzFRmod.mat');


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
frzMod=[];
for ratIdx=1:length(ratList)
    rat=ratList{ratIdx};
    frzMod=[frzMod;squeeze(frz.(rat).modIdx(:,[1,2,3,4,8]))];
        
end

%%
colTemp=setCoactColor();

piCol=[colTemp.cellType.inh;
       colTemp.cellType.nc;
       colTemp.cellType.ex];

targetReg={'PrL L5','BLA','vCA1'};
behList={'wake','nrem','rem'};
cellTypeList={'excitatory cells','inhibitory cells'};




sesName={'Baseline','Conditioning','Context','Cue ses. before first tone','Cue ses. after first tone'};

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
    
    eiFrzMod.mean=zeros(length(partnerList)+1,size(frzMod,2),2);
    eiFrzMod.ste=zeros(length(partnerList)+1,size(frzMod,2),2);
    eiFrzMod.raw=cell(length(partnerList)+1,size(frzMod,2),2);
    
    partnerName={};
    partnerNameCore={};
    pairCol=[];
    for n=1:length(partnerList)+1
        if n>length(partnerList)
            id=target(cellfun(@isempty,partner(target)));
            partnerName{n}='Others';
            partnerNameCore{n}='Others';
            pairCol(n,:)=colTemp.region.others;
            pairLeg{2*n-1}=sprintf('\\color[rgb]{%f %f %f}%s', pairCol(n,:),partnerNameCore{n});
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
            
            eiFrzMod.mean(n,:,eiIdx)=nanmean(frzMod(eiId{eiIdx},:));
            eiFrzMod.ste(n,:,eiIdx)=nanste(frzMod(eiId{eiIdx},:));
            for sesIdx=1:size(frzMod,2)
                eiFrzMod.raw{n,sesIdx,eiIdx}=frzMod(eiId{eiIdx},sesIdx)';
            end
        end
    end


    eiFrzMod.p.all=ones(size(frzMod,2),2);
    eiFrzMod.p.each=ones(size(frzMod,2),2,3,3);
    eiFrzMod.p.couple.all=ones(3,2);
    eiFrzMod.p.couple.each=ones(3,2,size(frzMod,2)*(size(frzMod,2)-1)/2,3);
    for sesIdx=1:size(frzMod,2)
        for eiIdx=1:2
            temp=cellfun(@(x) ones(size(x)),eiFrzMod.raw(:,sesIdx,eiIdx),'UniformOutput',false);
            grp=[];
            for regIdx=1:size(temp,1)
                grp=[grp,regIdx*temp{regIdx}]
            end
            val=cat(2,eiFrzMod.raw{:,sesIdx,eiIdx});
            [p,~,s]=kruskalwallis(val,grp,'off');
            eiFrzMod.p.all(sesIdx,eiIdx)=p;

            idx=0;
            nGrp=3;
            for n=1:nGrp-1
                for m=n+1:nGrp
                    idx=idx+1;
                    p=ranksum(eiFrzMod.raw{n,sesIdx,eiIdx},eiFrzMod.raw{m,sesIdx,eiIdx});
                    p=p*nGrp*(nGrp-1)/2;
                    eiFrzMod.p.each(sesIdx,eiIdx,idx,:)=[n,m,p];
                end
            end
        end
    end 
    for eiIdx=1:2
        for n=1:3
            temp=cat(1,eiFrzMod.raw{n,:,eiIdx})
            temp(:,any(isnan(temp),1))=[];
            eiFrzMod.p.couple.all(n,eiIdx)=friedman(temp',1,'off');

            idx=0;
            for sesIdx1=1:size(frzMod,2)-1
                for sesIdx2=sesIdx1+1:size(frzMod,2);
                    idx=idx+1;
                    p=signrank(temp(sesIdx1,:),temp(sesIdx2,:))*size(frzMod,2)*(size(frzMod,2)-1)/2;
                    eiFrzMod.p.couple.each(n,eiIdx,idx,:)=[sesIdx1,sesIdx2,p];
                end
            end
        end
    end
  
    % Freeze modulation
    xShift=0;
    yShift=(targetIdx-1)*(frzModheigth+yGapIntraBottom);

    for eiIdx=1:2
        subplotInMM(x+xShift+(frzModwidth+frzModGapXIntra)*(eiIdx-1),yTop+yShift,frzModwidth,frzModheigth)
        
        hold on
        posPool=[];
        poolAvg=[];

        for n=1:length(partnerName)
            errorbar((1:size(frzMod,2))+(2-n)*0.1,eiFrzMod.mean(n,:,eiIdx),eiFrzMod.ste(n,:,eiIdx),...
                'color',pairCol(n,:),'CapSize',0,'linewidth',0.5,'Marker','.','MarkerSize',4,'linestyle','none')

            posPool(n,:)=eiFrzMod.mean(n,:,eiIdx)+eiFrzMod.ste(n,:,eiIdx);
            poolAvg(n,:)=eiFrzMod.mean(n,:,eiIdx);
        end
     
        
        title([targetReg{targetIdx} ' ' cellTypeList{eiIdx}],'fontsize',5,'fontweight','normal')
        xlim([0.5,size(frzMod,2)+1])
        
        sigPosY=max(posPool(:));
        sigPosX=1:size(frzMod,2);
        curPos=0;
        ax=axis;
        stepSize=diff(ax(3:4))*0.1;
        sigTxtPool=cell(1,3);
        for n=length(partnerName):-1:1
            sigTxt=getSigTxt(eiFrzMod.p.couple.all(n,eiIdx));
            if ~isempty(sigTxt)
            sigTxtPool{n}=repmat('\dag',1,length(sigTxt));
                
                [sigPos,sigTxt]=findSigPos(squeeze(eiFrzMod.p.couple.each(n,eiIdx,:,:)));
                if ~isempty(sigPos)
                    for idx=1:size(sigPos,1)
                        plot(sigPosX(sigPos(idx,[1,1,2,2])), sigPosY+stepSize*((sigPos(idx,3)+curPos)+[0,0.5,0.5,0]),'-',...
                            'LineWidth',0.5,'color',pairCol(n,:))
                        text(mean(sigPosX(sigPos(idx,[1,2]))),sigPosY+stepSize*(sigPos(idx,3)+curPos+0.15),...
                            sigTxt{idx},...
                            'HorizontalAlignment','center','VerticalAlignment','baseline', 'Color',pairCol(n,:)...
                            )                        
                    end                
                    curPos=curPos+max(sigPos(:,3));
                end
            end
            
        end        
        
        ax=fixAxis;
        
        for sesIdx=1:size(frzMod,2)
            if eiFrzMod.p.all(sesIdx,eiIdx)<0.001
                sigTxt='###';
            elseif eiFrzMod.p.all(sesIdx,eiIdx)<0.01
                sigTxt='##';
            elseif eiFrzMod.p.all(sesIdx,eiIdx)<0.05
                sigTxt='#';
            else
                sigTxt='';
            end
            
            if ~isempty(sigTxt)
                text(sesIdx,max(posPool(:,sesIdx))+diff(ax(3:4))*0.05,sigTxt,'horizontalAlign','center','fontsize',5,'verticalAlign','middle')

                if any(eiFrzMod.p.each(sesIdx,eiIdx,:,3)<0.05)
                    [sigPos,sigTxt]=findSig(squeeze(eiFrzMod.p.each(sesIdx,eiIdx,:,:)),poolAvg(:,sesIdx));      
                    
                    for sigIdx=1:size(sigPos,1)
                        if sigPos(sigIdx,1)<0
                            sigX=[-0.18,-0.22,-0.22,-0.18];
                            sigTxtX=-0.24;
                            vAlign='middle';
                        else
                            sigX=[0.18,0.22,0.22,0.18];
                            sigTxtX=0.22;
                            vAlign='top';
                        end
                        sigY=sort(sigPos(sigIdx,3:4)).*[1.00,1.00];
                        
                        plot(sesIdx+sigX,sigY([1,1,2,2]),'k-','linewidth',0.5)
                        text(sesIdx+sigTxtX,mean(sigPos(sigIdx,3:4)),sigTxt{sigIdx},'fontsize',6,...
                            'Rotation',90,'VerticalAlignment',vAlign,'HorizontalAlignment','center')
                    end                   
                    
                end            
            end
                
        end 
        
        
            ylabel({'Freeze' 'modulation index'})
        for n=1:5
            text2(1,1-0.15*(n-1),pairLeg{n},ax,'verticalALign','top')
        end
        for n=1:3
            if ~isempty(sigTxtPool{n})
                text2(5/5.2,1-0.15*(2*n-2),sigTxtPool{n},ax,'verticalALign','top',...
                    'color',pairCol(n,:),'Interpreter','latex','horizontalALign','center')
            end
        end        
        set(gca,'xtick',1:size(frzMod,2),'XTickLabel',sesName,'XTickLabelRotation',25)
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



