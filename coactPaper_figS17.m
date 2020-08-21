function coactPaper_figS17()
labelWeight=true;
labelCase=true;
labelSize=9;
letGapX=9;
letGapY=3;

close all
fh=initFig('height',3);

x=11;y=5;
panel_01(x,y);
panelLetter2(x-letGapX,y-letGapY,alphabet(1,labelCase),'fontSize',labelSize,'isBold',labelWeight)
panelLetter2(x-letGapX+55,y-letGapY,alphabet(2,labelCase),'fontSize',labelSize,'isBold',labelWeight)
drawnow();

print(fh,'figS17.pdf','-dpdf','-painters','-r300')

end
%%

function panel_01(x,y)
%%
width=17.5;
height=15;
yGap=8;
xGap=5;
interGapX=9;
%%
coact=poolVar('coactCompCell.mat');

info=poolVar('okUnit.cellinfo.mat');
rip=poolVar('ripFrMod.mat');
spdl=poolVar('spdlFrMod.mat');
delta=poolVar('deltaFrMod.mat');
hfo=poolVar('hfoFrMod.mat');

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

%% get param in HC

hc.nSwrPart=[];
hc.nSwrGain=[];

for ratIdx=1:length(ratList)
    rat=ratList{ratIdx};
    
    if isfield(rip,rat)
        hc.nSwrPart=[hc.nSwrPart;cat(1,rip.(rat).hcNrem(2:3).participation)'*100];
        hc.nSwrGain=[hc.nSwrGain;cat(1,rip.(rat).hcNrem(2:3).gain)'*100];
        
    else
        temp=nan(size(info.(rat).channel,2),2);
        hc.nSwrPart=[hc.nSwrPart;temp];
        hc.nSwrGain=[hc.nSwrGain;temp];
    end
end
%%
hc.hfoPart=[];
hc.hfoGain=[];

for ratIdx=1:length(ratList)
    rat=ratList{ratIdx};
    
    hc.hfoPart=[hc.hfoPart;cat(1,hfo.(rat).hcNrem(2:3).participation)'*100];
    hc.hfoGain=[hc.hfoGain;cat(1,hfo.(rat).hcNrem(2:3).gain)'*100];
end
%%

hc.spdlPart=[];
hc.spdlGain=[];
for ratIdx=1:length(ratList)
    rat=ratList{ratIdx};
    
    hc.spdlPart=[hc.spdlPart;cat(1,spdl.(rat).pfc(2:3).participation)'*100];
    hc.spdlGain=[hc.spdlGain;cat(1,spdl.(rat).pfc(2:3).gain)'*100];
end

hc.preDelPart=[];
hc.preDelGain=[];
hc.postDelPart=[];
hc.postDelGain=[];
hc.delSpkPart=[];
for ratIdx=1:length(ratList)
    rat=ratList{ratIdx};
    
    hc.delSpkPart=[hc.delSpkPart;cat(1,delta.(rat).pfc.deltaSpike(2:3).participation)'*100];
    hc.preDelPart=[hc.preDelPart;cat(1,delta.(rat).pfc.preDelta(2:3).participation)'*100];
    hc.postDelPart=[hc.postDelPart;cat(1,delta.(rat).pfc.postDelta(2:3).participation)'*100];
    hc.preDelGain=[hc.preDelGain;cat(1,delta.(rat).pfc.preDelta(2:3).gain)'*100];
    hc.postDelGain=[hc.postDelGain;cat(1,delta.(rat).pfc.postDelta(2:3).gain)'*100];
end

%%
colTemp=setCoactColor();

col=[colTemp.cellType.inh
    colTemp.cellType.nc
    colTemp.cellType.ex];
%%
CTname={'Inhibitory','Not classified','Excitatory'};

mesList={'nSwrGain','hfoGain'};

mesNameList.nSwrPart={'SWR' 'participation (%)'};
mesNameList.nSwrGain={'FR modulation' 'by SWRs (%)'};
mesNameList.wSwrPart={'Wake SWR' 'participation (%)'};
mesNameList.wSwrGain={'Wake SWR' 'gain (%)'};
mesNameList.hfoPart={'HFO' 'participation (%)'};
mesNameList.hfoGain={'FR modulation' 'by HFOs (%)'};
mesNameList.spdlPart={'spindle' 'participation (%)'};
mesNameList.spdlGain={'spindle' 'gain (%)'};
mesNameList.delSpkPart={'Delta spike' 'participatio (%)'};
mesNameList.preDelPart={'Delta onset' 'participatio (%)'};
mesNameList.postDelPart={'Delta offset' 'participatio (%)'};
mesNameList.preDelGain={'Delta onset' 'gain (%)'};
mesNameList.postDelGain={'Delta offset' 'gain (%)'};

targetReg={'vCA1','BLA'};
partnerList={'PrL L5','BLA','vCA1'};

yRange.PrLL5.nSwrGain=[0,200];
yRange.PrLL5.hfoGain=[0,250];
yRange.PrLL5.delSpkPart=[0,4];

yRange.BLA.nSwrGain=[0,200];
yRange.BLA.hfoGain=[0,500];
yRange.BLA.delSpkPart=[0,4];

yRange.vCA1.nSwrGain=[0,700];
yRange.vCA1.hfoGain=[0,300];
yRange.vCA1.delSpkPart=[0,4];

yTick.PrLL5.nSwrGain=[];
yTick.PrLL5.hfoGain=[];
yTick.PrLL5.delSpkPart=[];

yTick.BLA.nSwrGain=[];
yTick.BLA.hfoGain=[0:200:400];
yTick.BLA.delSpkPart=[];

yTick.vCA1.nSwrGain=[0:300:600];
yTick.vCA1.hfoGain=[];
yTick.vCA1.delSpkPart=[];

cellLeg={};
for n=1:3
    cellLeg{n}=sprintf('\\color[rgb]{%f %f %f}%s',col(n,:),CTname{n})
end

for regIdx=1:length(targetReg)
    
    target=find(strcmp(reg,targetReg{regIdx}));
    
    patList=partnerList;
    patList(strcmp(patList,targetReg{regIdx}))=[];
    cList={};
    
    for n=1:length(patList)
        cList{n}=target(cellfun(@(x) any(strcmp(patList{n},x)),partner(target)));
    end
    
    cName=cellfun(@(x) ['' x], patList,'UniformOutput',false);
    
    
    for mesIdx=1:length(mesList)
        mesName=mesList{mesIdx};

        subplotInMM(x+(width+xGap)*(regIdx-1+2*(mesIdx-1))+interGapX*(mesIdx-1),y,width,height)
        hold on
        for CTidx=0:1
            
            avg=[];
            err=[];
            p=[];
            for n=1:length(cName)
                val=hc.(mesName)(cList{n},:);
                subCT=cellType(cList{n});
                subset=val(subCT==(1-2*CTidx),:);
                
                avg(n,:)=nanmean(subset,1);
                if size(subset,1)>1
                    err(n,:)=nanste(subset,[],1);
                    if any(all(~isnan(subset),2))
                        p(n)=signrank(subset(:,1),subset(:,2));
                    else
                        p(n)=nan;
                    end
                else
                    err(n,:)=[0,0];
                    p(n)=nan;
                end
            end
            hold on
            
            xVal{1}=(0:length(cName)-1)*(length(cName)+0.5)+1-0.2+CTidx;
            plot(xVal{1}+[0;0],(avg(:,1)+err(:,1).*[0,1])','Color',col(3-2*CTidx,:))
            bar(xVal{1},avg(:,1),0.1,'EdgeColor',col(3-2*CTidx,:),'facecolor',0.99*[1,1,1])
            
            xVal{2}=(0:length(cName)-1)*(length(cName)+0.5)+1+0.2+CTidx;
            plot(xVal{2}+[0;0],(avg(:,2)+err(:,2).*[0,1])','Color',col(3-2*CTidx,:))
            bar(xVal{2},avg(:,2),0.1,'EdgeColor',col(3-2*CTidx,:),'facecolor',col(3-2*CTidx,:))
            ax=fixAxis;
            for n=1:length(p)
                if p(n)<0.001
                    sig='***';
                    sigShift=0.02;
                elseif p(n)<0.01
                    sig='**';
                    sigShift=0;
                elseif p(n)<0.05
                    sig='*';
                    sigShift=0.02;
                else
                    continue
                end
                sigPosY=max((avg(n,:)+err(n,:)))+diff(yRange.(strrep(targetReg{regIdx},' ','')).(mesName))/20;
                plot((length(cName)+0.5)*(n-1)+CTidx+1+0.2*[-1,1],sigPosY+[0,0],'k-')
                text((length(cName)+0.5)*(n-1)+CTidx+1+sigShift,sigPosY,sig,'HorizontalAlignment','center','fontsize',7)
            end
            
        end
        
        if regIdx==1
            ylabel(mesNameList.(mesName),'fontsize',5,'fontweight','normal')
        end
        set(gca,'XTick',[])
        xlim([0.5,length(cName)*2+0.5+0.5])
        ylim(yRange.(strrep(targetReg{regIdx},' ','')).(mesName))
        if ~isempty(yTick.(strrep(targetReg{regIdx},' ','')).(mesName))
            set(gca,'YTick',yTick.(strrep(targetReg{regIdx},' ','')).(mesName))
        end
        ax=fixAxis;
        for cIdx=1:length(cName)
            text((length(cName)+0.5)*(cIdx-1)+1+0.5,ax(3)-diff(ax(3:4))*0.075,...
                cName{cIdx},'horizontalAlign','center','fontsize',5,'verticalAlign','top')
        end
        text(mean([0.5,length(cName)*2+0.5+0.5]),ax(3)-diff(ax(3:4))*0.3,'Coupled region',...
            'horizontalAlign','center','fontsize',5,'verticalAlign','top')

        title(['Cells in ' targetReg{regIdx}],'fontsize',5,'fontweight','normal')

        if regIdx==length(targetReg) && mesIdx==length(mesList)

            subplotInMM(x+(width+xGap)*(regIdx-1+2*(mesIdx-1))+interGapX*(mesIdx-1)+width+3,y,5,height)
            xlim([0,5])
            ylim([0,height])
            text(0,height,cellLeg{3},'verticalAlign','middle','fontsize',5)
            text(0,height-2,cellLeg{1},'verticalAlign','middle','fontsize',5)
            
            rectangle('Position',[0,height-5.25,2,1],'linestyle','-','EdgeColor','k')
            text(2.5,height-4,'Pre-','verticalAlign','middle','fontsize',5)
            text(2.5,height-5.5,'cond.','verticalAlign','middle','fontsize',5)
            
            rectangle('Position',[0,height-8.75,2,1],'linestyle','-','EdgeColor','k','FaceColor','k')
            text(2.5,height-7.5,'Post-','verticalAlign','middle','fontsize',5)
            text(2.5,height-9,'cond.','verticalAlign','middle','fontsize',5)
            
            axis off
        end
        
    end        
end

%%


end


