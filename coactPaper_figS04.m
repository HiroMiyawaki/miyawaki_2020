function coactPaper_figS04()
labelWeight=true;
labelCase=true;
labelSize=9;
letGapX=6;
letGapY=3;

close all
fh=initFig('height',6.5);

x=8.5;y=3;
panel_01_02(x,y);
panelLetter2(x-letGapX-2,y-letGapY++1,alphabet(1,labelCase),'fontSize',labelSize,'isBold',labelWeight)
panelLetter2(x-letGapX+52,y-letGapY+42,alphabet(2,labelCase),'fontSize',labelSize,'isBold',labelWeight)
drawnow();

print(fh,'figS04.pdf','-dpdf','-painters','-r300')

end
%%
function panel_01_02(x,y)
freeze=poolVar('freezeHMM.events.mat','subdir','','delimiter','.');
session=poolVar('sessions.events.mat','subdir','','delimiter','.');
cue=poolVar('cues.events.mat','subdir','','delimiter','.');
shock=poolVar('shocks.events.mat','subdir','','delimiter','.');

basicMetaData=poolVar('basicMetaData.mat','subdir','','delimiter','.');
%%
ratList=fieldnames(freeze);
sesName={basicMetaData.(ratList{1}).chamber.name};
col=setCoactColor;

shockCol=[1,0,0];
cueCol=[0.6,0.8,0.9];
lesionCol=[1,0.6,0];
hcCol=[0.8,0.8,1];
colChamA=[0,0.7,0];
colChamB=[0.95,0,0];

%%
plotW_max=115;
plotH=10;

headCap=@(x) [upper(x(2)),x(3:end)];

mmPerSec=(plotW_max-x)/diff(session.(ratList{1}).timestamps(4,:));
plotGapW=((plotW_max-x)-sum(diff(session.(ratList{1}).timestamps(1:3,:),1,2))*mmPerSec)/2;
plotGapH=10;
xPos=inf;
yPos=y-plotH-plotGapH;

for sIdx=1:length(sesName);
    
    frz=[];
    dur=[];
    for ratIdx=1:length(ratList);
        
        ratName=ratList{ratIdx};
        
        tRange=session.(ratName).timestamps(sIdx,:);
        
        subFrz=freeze.(ratName).timestamps(freeze.(ratName).timestamps(:,2)>tRange(1)&freeze.(ratName).timestamps(:,1)<tRange(2),:);
        subQ=cue.(ratName).timestamps.Tone(cue.(ratName).timestamps.Tone(:,2)>tRange(1)&cue.(ratName).timestamps.Tone(:,1)<tRange(2),:);
        %     subShock=shock.(ratName).timestamps.ShcokTrig(shock.(ratName).timestamps.ShcokTrig(:,1)>tRange(1)&shock.(ratName).timestamps.ShcokTrig(:,1)<tRange(2),1);
        subShock=shock.(ratName).timestamps.ShockTrig(shock.(ratName).timestamps.ShockTrig(:,1)>tRange(1)&shock.(ratName).timestamps.ShockTrig(:,1)<tRange(2),1);
        
        if ~isempty(subShock)
            subShock=subShock-tRange(1);
        end
        
        if subFrz(1)<tRange(1);subFrz(1)=tRange(1);end
        if subFrz(end)>tRange(2);subFrz(end)=tRange(2);end
        
        %     if sIdx==4
        %         tBorder=linspace(tRange(1),tRange(2),10);
        %     else
        %         tBorder=linspace(tRange(1),tRange(2),4);
        %     end
        tBorder=sort([tRange(1),subQ(:,1)',subQ(:,2)',tRange(2)]);
        
        intervals=[tBorder(1:2:end);tBorder(2:2:end)]'*[1,2;2,1]/3;
        tBorder=sort([tBorder,intervals(:)'])
        
        for tIdx=1:length(tBorder)-1;
            f=subFrz(subFrz(:,2)>tBorder(tIdx)&  subFrz(:,1)<tBorder(tIdx+1),:);
            
            if isempty(f)
                frz(ratIdx,tIdx)=0;
                continue
            end
            
            if f(1)<tBorder(tIdx);f(1)=tBorder(tIdx);end
            if f(end)>tBorder(tIdx+1);f(end)=tBorder(tIdx+1);end
            
            frz(ratIdx,tIdx)=sum(diff(f,1,2));
        end
        t=(tBorder(1:end-1)+tBorder(2:end))/2-tRange(1);
        tBorder=tBorder-tRange(1);
        dur(ratIdx,:)=diff(tBorder);
        subQ=subQ-tRange(1);
        
    end
    
    plotW=diff(tRange)*mmPerSec;
    if xPos+plotW>plotW_max
        xPos=x;
        yPos=yPos+plotGapH+plotH;
    else
        xPos=xPos+plotGapW;
    end
    subplotInMM(xPos,yPos,plotW,plotH,true)
    for tIdx=1:size(subQ,1);
        rectangle('position',[subQ(tIdx,1),0,diff(subQ(tIdx,:)),100],...
            'linestyle','none','facecolor',cueCol)
    end
    hold on
    if ~isempty(subShock)
        plot(subShock+[0,0],[0,100],'-','color',shockCol)
    end
    % for ratIdx=1:length(ratList)
    %     plot(t,frz(ratIdx,:)./dur(ratIdx,:)*100,'-','color',col(ratIdx,:),'linewidth',0.5)
    % end
    avg=mean(frz./dur*100);
    er=ste(frz./dur*100);
    fill([t,fliplr(t)],[avg+er,fliplr(avg-er)],0.8*[1,1,1],'LineStyle','none')
    plot(t,avg,'k-','linewidth',0.5,'markersize',6)
    xlim(tRange-tRange(1))
    ylim([0,100])
    % if  xPos==x+plotGapW;
    ylabel({'Freeze (%)'},'fontsize',5,'fontweight','normal')
    % end
    xlabel('Time (sec)')
    ax=fixAxis;
    set(gca,'TickLength',0.5/plotW*[1,1],'TickDir','out')

    text2(0,1.05,headCap(lower(regexprep(sesName{sIdx},'[A-Z]',' $0'))),ax,'horizontalAlign','left','verticalAlign','bottom','fontweight','normal','fontsize',5)
    xPos=xPos+plotW;
    
    if sIdx~=3
        legTxt=sprintf('\\color[rgb]{%f %f %f} Cue ',cueCol),
        
        if sIdx==2
            legTxt=[legTxt sprintf('\\color[rgb]{%f %f %f}Shock ',shockCol)]
        end
        text2(1,1.01,legTxt,ax,'verticalAlign','bottom','horizontalALign','right')
        
    end
    %     if sIdx==length(sesName)
    %     legTxt={};
    %     for ratIdx=1:length(ratList)
    %         legTxt{ratIdx}=['\color[rgb]{' num2str(col(ratIdx,:)) '}' headCap([' ' ratList{ratIdx}(1:end-6)])];
    %     end
    %     legTxt{end+1}='\color[rgb]{0 0 0}Mean';
    %     text2(1.02,1,legTxt,ax,'horizontalAlign','left','verticalAlign','top')
    %     end
    
end

xPos=xPos+plotGapW;
subplotInMM(xPos,yPos,plotW_max-xPos,plotH)
set(gca,'ytick',[])
xlabel('Time (min)')
conCol=[hcCol;colChamA;colChamB];
rectangle('Position',[basicMetaData.(ratName).detectionintervals.lfp(1)/60,0,diff(basicMetaData.(ratName).detectionintervals.lfp)/60,1],...
    'linestyle','none','facecolor',conCol(1,:))
cType=[2,3,3,2,2];
for sIdx=1:size(session.(ratName).timestamps,1)
    rectangle('position',[session.(ratName).timestamps(sIdx,1)/60,0,diff(session.(ratName).timestamps(sIdx,:)/60),1],...
        'linestyle','none','facecolor',conCol(cType(sIdx),:))
    if mod(sIdx,2)==0
        text(mean(session.(ratName).timestamps(sIdx,:)/60),1.05,headCap(lower(regexprep(sesName{sIdx},'[A-Z]',' $0'))),...
            'verticalALign','bottom','color',conCol(cType(sIdx),:),'fontsize',5,'horizontalALign','center')
    else
        text(mean(session.(ratName).timestamps(sIdx,:)/60),-0.1,headCap(lower(regexprep(sesName{sIdx},'[A-Z]',' $0'))),...
            'verticalALign','top','color',conCol(cType(sIdx),:),'fontsize',5,'horizontalALign','center')
        
    end
end
hold on
plot(basicMetaData.(ratName).detectionintervals.lfp(1)/60+[0,0],[0,1],'k-')
text(basicMetaData.(ratName).detectionintervals.lfp(1)/60,1.05,'Start recording','VerticalAlignment','bottom','HorizontalAlignment','left')
plot(basicMetaData.(ratName).detectionintervals.lfp(2)/60+[0,0],[0,1],'k-')
text(basicMetaData.(ratName).detectionintervals.lfp(2)/60,1.05,'End recording','VerticalAlignment','bottom','HorizontalAlignment','right')

% text(basicMetaData.(ratName).detectionintervals.lfp(2)/60+30,-0.1,...
%     'Lesioning','VerticalAlignment','top','HorizontalAlignment','right','color',lesionCol)
text(basicMetaData.(ratName).detectionintervals.lfp(2)/60+35,0.5,...
    'Lesioning','VerticalAlignment','bottom','HorizontalAlignment','center','color',lesionCol,'rotation',-90)


rectangle('Position',[basicMetaData.(ratName).detectionintervals.lfp(2)/60+20,0,10,1],...
    'linestyle','none','facecolor',lesionCol)

ylim([-1,2.5])
xlim(basicMetaData.(ratName).detectionintervals.lfp/60+[-10,35])
ax=fixAxis;
text2(0,1.1,'Behavior schedule',ax,'fontweight','normal','fontsize',5)
set(gca,'TickLength',0.5/(plotW_max-xPos)*[1,1],'TickDir','out')
box off

yExtend=2;
subplotInMM(xPos,yPos-yExtend,plotW_max-xPos,plotH+yExtend,true,true)

chamName={'Homecage','Chamber A','Chamber B'};
legTxt=[];
for n=1:3
    rectangle('position',[200+200*n+10,2,180,0.8],'linestyle','none','facecolor',conCol(n,:))
    text(200+200*(n+0.5),2.4,chamName{n},'horizontalALign','center','verticalAlign','middle','fontsize',5)
end
ylim([0,3.5*(plotH+yExtend)/plotH]-1)
xlim(basicMetaData.(ratName).detectionintervals.lfp/60+[-10,35])
axis off
end
%%
