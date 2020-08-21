function coactPaper_figS01()
labelWeight=true;
labelCase=true;
labelSize=9;
letGapX=5;
letGapY=4;
%
close all
fh=initFig('height',13.5);

x=4;y=7;
panel_01(x,y);
panelLetter2(x-letGapX+2,y-letGapY,alphabet(1,labelCase),'fontSize',labelSize,'isBold',labelWeight)
drawnow();

x=10;y=7+24;
panel_02(x,y);
panelLetter2(x-letGapX-4,y-letGapY+4,alphabet(2,labelCase),'fontSize',labelSize,'isBold',labelWeight)
drawnow();

x=10;y=7+24+35;
panel_03(x,y);
panelLetter2(x-letGapX-4,y-letGapY+5,alphabet(3,labelCase),'fontSize',labelSize,'isBold',labelWeight)
drawnow();

x=10;y=7+24+35+45;
panel_04(x,y);
panelLetter2(x-letGapX-4,y-letGapY+3,alphabet(4,labelCase),'fontSize',labelSize,'isBold',labelWeight)

print(fh,'figS01.pdf','-dpdf','-painters','-r300')

end

%%

function panel_01(x,y)

scaleFactor=0.66;

hpcFile='~/Dropbox/FearAnalyses/histo_example_vCA1.png';
amyFile='~/Dropbox/FearAnalyses/histo_example_BLA.png';
pfcFile='~/Dropbox/FearAnalyses/histo_example_PFC.png';

colLeg='\color[rgb]{1,0,0}Nissl \color[rgb]{0,0,1}DAPI';

subplotInMM(x,y,40*scaleFactor,31*scaleFactor)
image(imread(hpcFile))
axis equal
axis off
ax=fixAxis;
text2(0,-0.06/scaleFactor,'Probe in the ventral hippocampus',ax,'verticalALign','bottom','fontweight','normal','fontsize',5)
for n=0:2
    text2(13.5/40*n,-0.002/scaleFactor, sprintf('%0.2f mm',-5.00-0.05*n),ax,'verticalALign','bottom','fontweight','normal','fontsize',5)
end
text2(1,0,colLeg,ax,'verticalALign','bottom','fontweight','normal','fontsize',5,'rotation',-90)

subplotInMM(x+40*scaleFactor+3,y,40*scaleFactor,31*scaleFactor)
image(imread(amyFile))
axis equal
axis off
ax=fixAxis;
text2(0,-0.06/scaleFactor,'Probe in the amygdala',ax,'verticalALign','bottom','fontweight','normal','fontsize',5)
for n=0:2
    text2(13.5/40*n,-0.002/scaleFactor, sprintf('%0.2f mm',-2.55-0.05*n),ax,'verticalALign','bottom','fontweight','normal','fontsize',5)
end
text2(1,0,colLeg,ax,'verticalALign','bottom','fontweight','normal','fontsize',5,'rotation',-90)

subplotInMM(x+40*2*scaleFactor+3*2,y,80.5*scaleFactor,31*scaleFactor)
image(imread(pfcFile))
ax=fixAxis;
text2(0,-0.06/scaleFactor,'Probe in the prefrontal cortex',ax,'verticalALign','bottom','fontweight','normal','fontsize',5)
axis equal
axis off
for n=0:5
    text2(13.5/80.5*n,-0.002, sprintf('+%0.2f mm',4.20-0.2*n),ax,'verticalALign','bottom','fontweight','normal','fontsize',5)
end
text2(1,0,colLeg,ax,'verticalALign','bottom','fontweight','normal','fontsize',5,'rotation',-90)


end
function panel_02(x,y)
doUpdate=false;

offset=47.9;
dur=1.5;


tRange=offset+[0,dur];

chGap=50;
linewidth=0.25;



col=setCoactColor;
colList=[0.6,0.6,0.6; %PFC
    0.8,0.8,0.8; %BLA
    0.6,0.6,0.6; %HPC
    1.0,0.5,0.0; %EMG
    0.0,0.6,1.0; %OB
    1.0,0.0,0.7; %ECG
    0,0.8,0.6]; %Acc

txtCol=colList;
txtCol(1,:)=col.region.vCA1;
txtCol(2,:)=col.region.BLA;
txtCol(3,:)=col.region.PrLL5;


probe=[1,3,2];
chList={(1:64)+(probe(1)-1)*64,(1:64)+(probe(2)-1)*64,(1:64)+(probe(3)-1)*64,193,195,194,196};
tracName={'vCA1','BLA','PrL L5','EMG','EOG','ECG','Head acc'};
resGap=[0,2,2,8,8,8,8];
scale=[1,1,1,0.5,1/0.2,0.5,500];
scaleBarGap=[0,0,0,(-3:0)*8];
unit={'mV','mV','mV','mV','mV','mV','m/s^2'};
unitScale=[1000,1000,1000,1000,1000,1000,1];
barSize=1000;
width=163;
height=50;
if ~doUpdate && exist('~/Dropbox/FearAnalyses/png/example-trace.png','file') && ...
        exist('~/Dropbox/FearAnalyses/png/example-trace-info.mat','file')
    
else
    
    load('~/data/Fear/triple/hoegaarden181115/hoegaarden181115.basicMetaData.mat')
    
    if exist('~/Dropbox/FearData/example/example-trace.mat','file')
        load('~/Dropbox/FearData/example/example-trace.mat')
    else
        if ~exist(basicMetaData.dat,'file')
            subplotInMM(x,y,width,height)
            return
        end
        
        load('~/data/Fear/triple/hoegaarden181115/hoegaarden181115.acceleration.lfp.mat');
        nCh=basicMetaData.nCh;
        samplingRate=basicMetaData.SampleRates.dat;
        dur=60;
        tOffset=359*60+45-dur/2;
        
        fh=fopen(basicMetaData.dat);
        fseek(fh,tOffset*samplingRate*nCh*2,'bof');
        eeg=fread(fh,[nCh,dur*samplingRate+1],'int16');
        fclose(fh);
        
        emgScale=0.195;
        accScale=0.0012; %m/s^2 per bit
        
        param.setting.dat=basicMetaData.dat;
        param.setting.nCh=nCh;
        param.setting.dur=dur;
        param.setting.tOffset=tOffset;
        param.setting.emgScale=emgScale;
        param.setting.accScale=accScale;
        
        param.generatedate=today('datetime');
        param.generator=mfilename;
        
        param.unit.t='s';
        param.unit.lfp='/muV';
        param.unit.acc='m/s^2';
        
        
        
        basicMetaData.Ch.names(1:195)
        
        lfp=eeg(1:195,:)*emgScale;
        xyzAcc=eeg(196:198,:)*accScale;
        
        t=(0:dur*samplingRate)/samplingRate;
        
        acc=interp1(...
            accelerometer.timestamps(accelerometer.timestamps>=tOffset&accelerometer.timestamps<=tOffset+dur)-tOffset,...
            accelerometer.abs(accelerometer.timestamps>=tOffset&accelerometer.timestamps<=tOffset+dur),t);
        if ~exist('~/Dropbox/FearData/example','dir')
            mkdir('~/Dropbox/FearData/example')
        end
        save('~/Dropbox/FearData/example/example-trace.mat','lfp','acc','xyzAcc','t','param','-v7.3')
    end
    %%
    load('~/data/Fear/triple/hoegaarden181115/hoegaarden181115.okUnit.spikes.mat')
    lfp(196,:)=acc;
    
    %%
    
    fhTemp=figure();
    set(fhTemp, 'paperUnit','centimeters','Units','centimeters')
    set(fhTemp,'position',[0,20,width/10,height/10])
    set(fhTemp, 'Units','centimeters')
    set(fhTemp,'PaperSize',[width/10,height/10])
    set(fhTemp,'paperPosition',[0,0,width/10,height/10])
    
    set(fhTemp,'defaultAxesFontName','Helvetica')
    set(fhTemp,'defaultTextFontName','Helvetica')
    
    set(fhTemp,'defaultAxesXColor',[0,0,0]);
    set(fhTemp,'defaultAxesYColor',[0,0,0]);
    set(fhTemp,'defaultAxesZColor',[0,0,0]);
    
    set(fhTemp,'defaultAxesFontSize',5);
    set(fhTemp,'defaultTextFontSize',5);
    set(fhTemp,'defaultAxesLineWidth', 0.5);
    
    
    set(fhTemp,'DefaultAxesXGrid','off');
    set(fhTemp,'DefaultAxesYGrid','off');
    set(fhTemp,'DefaultAxesBox','off');
    set(fhTemp,'PaperPositionMode','auto');
    
    subplot('Position',[0,0,1,1]);
    hold on

    idx=(okUnit.spikeTime>param.setting.tOffset+offset & okUnit.spikeTime<param.setting.tOffset + offset +dur);
    res=okUnit.spikeTime(idx)-param.setting.tOffset-offset;
    origClu=okUnit.cluster(idx);
    [cluList,~,clu]=unique(origClu);
    sh=okUnit.cluInfo.shank(cluList);
    
    threshold=50;
    for idx=1:length(cluList)
        amp=max(abs(okUnit.waveform.wave(cluList(idx)).mean),[],2);
        if any(amp*0.192>threshold)
            onCh{idx}=find(amp*0.192>threshold);
        else
            [~,onCh{idx}]=max(amp);
        end
            
    end
    
    colSpk=zeros(length(cluList),3);
    rng(1);
    for n=1:3
        idx=(sh>7*(n-1)&sh<=7*n);
        temp=jet(sum(idx));
        temp=temp(randperm(size(temp,1)),:);
        colSpk(idx,:)=temp;
    end
    
    idx=(okUnit.spikeTime>param.setting.tOffset+offset & okUnit.spikeTime<param.setting.tOffset + offset +dur);
    res=okUnit.spikeTime(idx)-param.setting.tOffset-offset;
    origClu=okUnit.cluster(idx);
    [cluList,~,clu]=unique(origClu);
    sh=okUnit.cluInfo.shank(cluList);
    
    colSpk=zeros(length(cluList),3);
    rng(1);
    for n=1:3
        idx=(sh>7*(n-1)&sh<=7*n);
        temp=jet(sum(idx));
        temp=temp(randperm(size(temp,1)),:);
        colSpk(idx,:)=temp;
    end
    showIdx=t>tRange(1) & t<=tRange(2);
    
    
    nCnt=0;
    n=1;
    scaleLet={['\color[rgb]{0.6,0.6,0.6}' num2str(barSize/unitScale(n)/scale(n)) ' ' unit{n}]}
    for n=1:length(chList)
        subLFP=lfp(chList{n},showIdx)*scale(n);
        tSub=t(showIdx);
        plot(tSub,subLFP-((1:length(chList{n}))'+nCnt+sum(resGap(1:n)))*chGap,...
            '-','linewidth',linewidth,'color',colList(n,:))
        
        if n<4
            targetClu=find(sh>7*(probe(n)-1)&sh<=7*(probe(n)));
            for tarIdx=1:length(targetClu)
                loc=round(res(clu==targetClu(tarIdx))*20e3);
                chIdx=basicMetaData.chMap{sh(targetClu(tarIdx))};
                
                chIdx=chIdx(onCh{targetClu(tarIdx)})
                
                for idx=1:length(loc)
                    if loc(idx)<9
                        f(1)=1;
                    else
                        f(1)=loc(idx)-8;
                    end
                    
                    if loc(idx)>length(tSub)-8
                        f(2)=length(tSub);
                    else
                        f(2)=loc(idx)+8;
                    end
                    
                    plot(tSub(f(1):f(2)),subLFP(chIdx-64*(probe(n)-1),f(1):f(2))-(chIdx'-64*(probe(n)-1)+nCnt+sum(resGap(1:n)))*chGap,...
                        '-','linewidth',linewidth,'color',colSpk(targetClu(tarIdx),:))
                    
                end
            end
        end
        if n>3
            scaleLet{end+1}=['\color[rgb]{' num2str(colList(n,:)) '}' num2str(barSize/unitScale(n)/scale(n)) ' ' unit{n}]
        end
        nCnt=nCnt+length(chList{n});
    end
    n=1;
    axis tight
    axis off
    xRange=get(gca,'xlim')
    yRange=get(gca,'ylim')
    xRange(2)=tRange(2)+dur/100;
    yRange(1)=-(nCnt+sum(resGap)+8)*chGap;
    xlim(xRange)
    ylim(yRange)
    
    print(fhTemp,'~/Dropbox/FearAnalyses/png/example-trace.png','-dpng','-r600')
    save('~/Dropbox/FearAnalyses/png/example-trace-info.mat','xRange','yRange')
    close(fhTemp)
end
%%
img=imread('~/Dropbox/FearAnalyses/png/example-trace.png');
info=load('~/Dropbox/FearAnalyses/png/example-trace-info.mat','xRange','yRange')
scaleFactor=0.645;
subplotInMM(x,y,width*scaleFactor,height*scaleFactor)
image(info.xRange,info.yRange,img)
nCnt=0;
n=1
scaleLet={['\color[rgb]{0.6,0.6,0.6}' num2str(barSize/unitScale(n)/scale(n)) ' ' unit{n}]}
for n=1:length(chList)
    if n>3
        scaleLet{end+1}=['\color[rgb]{' num2str(colList(n,:)) '}' num2str(barSize/unitScale(n)/scale(n)) ' ' unit{n}]
    end
    text(tRange(1),info.yRange(1)+(nCnt+sum(resGap(1:n))+4)*chGap,tracName{n},'horizontalAlign','right','fontsize',5,'color',txtCol(n,:));
    nCnt=nCnt+length(chList{n});
end
text(tRange(2)-dur/20-0.1/2,info.yRange(1)+(nCnt+sum(resGap)+8)*chGap,'100 ms',...
    'verticalALign','bottom','horizontalAlign','center')
n=1;
text(tRange(2)+dur/100*1.5,info.yRange(1)+(scaleBarGap(n)+nCnt+sum(resGap(1:n))-12)*chGap-barSize/2+1000,...
    scaleLet)

hold on
plot(tRange(2)+dur/100+[0,0],info.yRange(1)-(-(scaleBarGap(n)+nCnt+sum(resGap(1:n))+8)*chGap+[0,-barSize]+1000),...
    '-','color','k','linewidth',1)
plot(tRange(2)-dur/20-[0,0.1],info.yRange(1)-(-(nCnt+sum(resGap)+8)*chGap+[0,0]),'k-','linewidth',1)

axis off
end

%%
function panel_03(x,y)
basename='~/data/Fear/triple/hoegaarden181115/hoegaarden181115';
load([basename '.basicMetaData.mat'])

load([basicMetaData.Basename '.SleepState.states.mat'])
load([basicMetaData.Basename '.sessions.events.mat'])
load([basicMetaData.Basename '.okUnit.spikes.mat'])

tBin=basicMetaData.detectionintervals.lfp(1):10:basicMetaData.detectionintervals.lfp(2);
cBin=unique(okUnit.cluster);
cBin=[cBin;max(cBin)+1];

cnt=histcounts2(okUnit.spikeTime,okUnit.cluster,tBin,cBin);
fr{3}=cnt(:,strcmpi(okUnit.cluInfo.region,'PrL L5'))/10;
fr{2}=cnt(:,strcmpi(okUnit.cluInfo.region,'BLA'))/10;
fr{1}=cnt(:,strcmpi(okUnit.cluInfo.region,'vCA1'))/10;

%%
beh=relabel_ma2sleep(SleepState.MECE.timestamps);
col=setCoactColor;

width=100;
hight=7;

yGap=0;
cellH=0.15

frTick=[0.01,0.1,1,10,100];
rName={'vCA1','BLA','PrL L5'};
for n=1:3
    nCell=size(fr{n},2)
    subplotInMM(x,y+yGap,width,nCell*cellH)
    [~,order]=sort(mean(fr{n}),'descend')
    imagesc(tBin/60,[],(log10(fr{n}(:,order)')))
    box off
    
    set(gca,'XTick',[])
    cLim=get(gca,'clim');
    set(gca,'ytick',[])
    ylabel(rName{n})
    xlim(tBin([1,end])/60)
    colormap(gca,col.fr)
    
    subplotInMM(x+width+0.5,y+yGap,1,nCell*cellH,true)
    imagesc([],cLim,linspace(cLim(1),cLim(2),256)')
    box off
    set(gca,'YAxisLocation','right','yDir','normal','xtick',[])
    set(gca,'ytick',log10(frTick),'YTickLabel',frTick)
    
    yGap=yGap+nCell*cellH+1.5;
    if n==2
        ylabel('Firing rate (Hz)','fontsize',5)
        set(get(gca,'ylabel'),'Rotation',-90,'Position',get(get(gca,'ylabel'),'Position')+[2,0,0])
    end
    colormap(gca,col.fr)
    
end

subplotInMM(x,y+yGap,width,hight,true)
for idx=1:size(beh,1)
    rectangle('Position',[beh(idx,1)/60,4-(beh(idx,3)+1)/2,diff(beh(idx,1:2)/60),1],...
        'linestyle','none','facecolor','k')
end
chType=[1,2,2,1,1];
chCol=[0,0.8,0;1,0,0];
for idx=1:5
    rectangle('Position',[sessions.timestamps(idx,1)/60,0,diff(sessions.timestamps(idx,1:2)/60),1],...
        'linestyle','none','facecolor',chCol(chType(idx),:))
end
xlim(tBin([1,end])/60)
set(gca,'YTick',0.5:3.5,'yticklabel',{'Chamber','REM','NREM','WAKE'})
xlabel('Time (min)')
set(gca,'TickDir', 'out','TickLength',[0.004,0.000])

end

%%
function panel_04(x,y)
xGap=3;
width=(100-xGap)/2;
height=20;

doUpdate=false;

col=setCoactColor;
behList=fieldnames(col.state);
for n=1:3
    beh=behList{n};
    legSlp{n}=sprintf('\\color[rgb]{%f %f %f}%s',col.state.(beh),upper(beh))
end

if ~doUpdate &&...
    exist('~/Dropbox/FearAnalyses/png/example-reac-pre.png','file') && ...
    exist('~/Dropbox/FearAnalyses/png/example-reac-pre-info.mat','file')&& ...
    exist('~/Dropbox/FearAnalyses/png/example-reac-post.png','file') && ...
    exist('~/Dropbox/FearAnalyses/png/example-reac-post-info.mat','file')
else
          
basename='~/data/Fear/triple/hoegaarden181115/hoegaarden181115';

load([basename '.basicMetaData.mat'])
load([basicMetaData.Basename '.sleepstate.states.mat'])
load([basicMetaData.Basename '.sessions.events.mat'])
load([basicMetaData.AnalysesName '-icaReac.mat'])
load([basicMetaData.AnalysesName '-icaCoactTimeCondHT.mat'])
load([basicMetaData.AnalysesName '-icaReacZNCCG_sig.mat'])

%%
slp=relabel_ma2sleep(SleepState.MECE.timestamps);
slp(:,3)=(slp(:,3)+1)/2;
slp(:,1:2)=slp(:,1:2)/60;
slpCol=[col.state.wake;
    col.state.nrem;
    col.state.rem];


tempSes=icaCoactTimeCond.param.templateIdx;
tBin=(1:size(icaReac(tempSes).strength,2))*20e-3/60;

withSlpState=true;


exID=[6,8,14,25,28,33];
regName=icaReac(tempSes).region(exID);
regList=unique(regName);
num=zeros(size(regList))
for regID=1:length(regName)
    idx=find(strcmp(regName{regID},regList));
    num(idx)=num(idx)+1;
    enName{regID}=[regName{regID} '_{En' num2str(num(idx)) '}'];
end

dur=12;
tRange=[287;471]+[0,dur];

reac=zscore(icaReac(tempSes).strength(exID,:),[],2);

yGapUit=ceil(diff(prctile(reac(:),[0.01,99.99]))/10)*10;

yGapStep=[0:-1:-length(exID)+1];

[ytick,order]=sort(yGapStep)
ytickLabel=enName(order);

tTxt={'Pre-conditioning homecage session','Post-conditioning homecage session'};

for prePost=1:2

    fhTemp=figure();
    set(fhTemp, 'paperUnit','centimeters','Units','centimeters')
    set(fhTemp,'position',[0,20,width/10,height/10])
    set(fhTemp, 'Units','centimeters')
    set(fhTemp,'PaperSize',[width/10,height/10])
    set(fhTemp,'paperPosition',[0,0,width/10,height/10])
    
    set(fhTemp,'defaultAxesFontName','Helvetica')
    set(fhTemp,'defaultTextFontName','Helvetica')
    
    set(fhTemp,'defaultAxesXColor',[0,0,0]); 
    set(fhTemp,'defaultAxesYColor',[0,0,0]);
    set(fhTemp,'defaultAxesZColor',[0,0,0]);
    
    set(fhTemp,'defaultAxesFontSize',5);
    set(fhTemp,'defaultTextFontSize',5);
    set(fhTemp,'defaultAxesLineWidth', 0.5);
    subplot('Position',[0,0,1,1]);
    
        subSlp=slp(slp(:,2)>tRange(prePost,1)&slp(:,1)<tRange(prePost,2),:);
        if subSlp(1,1)<tRange(prePost,1);subSlp(1,1)=tRange(prePost,1);end
        if subSlp(end,2)>tRange(prePost,2);subSlp(end,2)=tRange(prePost,2);end
        subSlp(:,1:2)=subSlp(:,1:2)-tRange(prePost,1);
        if withSlpState
            for sIdx=1:size(subSlp,1)
                rectangle('Position',[subSlp(sIdx,1),yGapUit*(min(yGapStep)-2),diff(subSlp(sIdx,1:2)),yGapUit*(-min(yGapStep)+4)],'LineStyle','none','FaceColor',slpCol(subSlp(sIdx,3),:))
            end
        end
        toShow=(tBin>=tRange(prePost,1)&tBin<=tRange(prePost,2));
        hold on
        plot(tBin(toShow)-tRange(prePost,1), reac(:,toShow)+yGapUit*yGapStep','k-','linewidth',0.25)
        ylim(yGapUit*[min(yGapStep)-2,2])

        if prePost==1
            for n=1:length(ytickLabel)
                text(0,yGapUit*ytick(n),ytickLabel{n},'horizontalAlign','right','fontsize',5)
            end
        end
        
        axis off

        xRange=get(gca,'XLim');
        yRange=get(gca,'YLim');
        tText=tTxt{prePost};
        if prePost==1
            fName='pre'
        else
            fName='post'
        end
        
        print(fhTemp,['~/Dropbox/FearAnalyses/png/example-reac-' fName '.png'],'-dpng','-r600')
        save(['~/Dropbox/FearAnalyses/png/example-reac-' fName '-info.mat'],...
            'xRange','yRange','tText','dur','ytickLabel','ytick','yGapUit')
        close(fhTemp)
        
    end
    
    
end
im{1}=imread('~/Dropbox/FearAnalyses/png/example-reac-pre.png');
info{1}=load('~/Dropbox/FearAnalyses/png/example-reac-pre-info.mat');
im{2}=imread('~/Dropbox/FearAnalyses/png/example-reac-post.png');
info{2}=load('~/Dropbox/FearAnalyses/png/example-reac-post-info.mat');

ytick=info{1}.ytick;
yGapUit=info{1}.yGapUit;
scalePos=1.1;
for prePost=1:2
    subplotInMM(x+(prePost-1)*(width+xGap),y,width,height)
    image(info{prePost}.xRange,info{prePost}.yRange,flipud(im{prePost}))
    set(gca,'YDir','normal')
    hold on
    plot(info{prePost}.dur*0.8+[0,1],yGapUit*(min(ytick)-scalePos)+[0,0],'k-','LineWidth',1)
    text(info{prePost}.dur*0.8+0.5,yGapUit*(min(ytick)-scalePos),'1 min','horizontalAlign','center','verticalAlign','top')
    plot(info{prePost}.dur*0.8+1.05+[0,0],yGapUit*(min(ytick)-scalePos)+[0,30],'k-','LineWidth',1)
    text(info{prePost}.dur*0.8+1.05,yGapUit*(min(ytick)-scalePos)+15,' 30 z','horizontalAlign','left','verticalAlign','middle')
    
    if prePost==1
        for n=1:length(info{prePost}.ytickLabel)
            text(info{prePost}.xRange(1)-diff(info{prePost}.xRange)*0.01,yGapUit*ytick(n),info{prePost}.ytickLabel{n},'horizontalAlign','right','fontsize',5)
        end
    end
        
        axis off
        title(info{prePost}.tText)
        ax=fixAxis;
        if prePost==2
            text2(1.01,1.0,legSlp,ax,'verticalAlign','top','fontsize',5)
        end
end

end

%%
function panel_04old(x,y);
xGap=2;
width=((162-10*2)/3-xGap*2-3)/2;
yGap=1.5;
totalHeight=60;
eachPETH=poolVar('swrTrigIcaReac.mat');

ratList=fieldnames(eachPETH);

col=setCoactColor;
%%
smSigma=20/1000; %in s
cLim=4*[-1,1];


sesIdx=2;
tempName=eachPETH.(ratList{1})(sesIdx).tempName;
tBinSize=eachPETH.(ratList{1})(sesIdx).tBinSize;
tWin=eachPETH.(ratList{1})(sesIdx).param.tWindow;
nWin=ceil(tWin/tBinSize);
tBin=(-nWin:nWin)*tBinSize*1e3;
peth=[];
reg={};
for ratIdx=1:length(ratList)
    rat=ratList{ratIdx};
    peth=cat(1,peth,eachPETH.(rat)(sesIdx).nrem.mean);
    reg=[reg,eachPETH.(rat)(sesIdx).region];
end
%%
[reg,regList]=relabel_region(reg);
regID=zeros(size(reg));
for regIdx=1:length(regList)
    regID(strcmp(reg,regList{regIdx}))=regIdx;
end

nPair=sum(ismember(reg,{'vCA1','BLA','PrL L5'}));
eachHeight=(totalHeight-yGap*2)/nPair;


xRange=[-420,420];

periodText={'Pre','Post'};
xTxt={'Time form' 'SWR peak (ms)'};

nSm=ceil(smSigma/tBinSize);
xSm=(-nSm*4:nSm*4)*tBinSize;
smCore=normpdf(xSm,0,smSigma);
smCore=smCore/sum(smCore);


totalY=0;
for n=1:3
    switch n
        case 1
            target='vCA1';
        case 2
            target='BLA';
        case 3
            target='PrL L5';
        otherwise
            continue
    end
    idx=find(strcmp(reg,target));
    
    temp=zscore([peth(idx,:,1)';peth(idx,:,2)']);
    subsetZ{1}=Filter0(smCore,temp(1:nWin*2+1,:));
    subsetZ{2}=Filter0(smCore,temp(nWin*2+1+1:end,:));
    
    peak=mean(subsetZ{2}(nWin+1+(-1:1),:),1);
    [~,order]=sort(peak);
    
    height=length(idx)*eachHeight;
    for prePost=1:2
        
        subplotInMM(x+(width+xGap)*(prePost-1),y+totalY,width,height)
        
        imagescXY(tBin,[],subsetZ{prePost}(:,order));
        box off
        xlim(xRange)
        set(gca,'YTick',[])
        set(gca,'clim',cLim)
        ax=fixAxis;
        if n==3
            xlabel(xTxt)
        else
            set(gca,'XTickLabel',[])
        end
        if n==1
            text2(0.5,1,{sprintf('%s-conditoning',periodText{prePost}),'NREM'},ax,...
                'fontsize',5,'horizontalAlign','center','verticalAlign','bottom')
        end
        if prePost==1
            ylabel(sprintf('%s ensemble',target),'fontsize',5,'fontweight','normal')
        end
        colormap(gca,col.react.map)
    end
    totalY=totalY+height+yGap;
end
subplotInMM(x+width*2+xGap*2,y,1,totalY-yGap)
imagescXY([0,1],cLim,linspace(cLim(1),cLim(2),256));
box off
set(gca,'XTick',[],'YTick',[cLim(1),cLim(1)/2,0,cLim(2)/2,cLim(2)],'YAxisLocation','right')
colormap(gca,col.react.map)
ax=fixAxis;
text2(5,0.5,'Reactivation strength (z)',ax,'fontsize',5,'horizontalALign','center','rotation',-90)
end

%%
function panel_05old(x,y);
xGap=2;
width=((162-10*2)/3-xGap*2-3)/2;
yGap=1.5;
eachHeight=0.32;
totalHeight=60;
eachPETH=poolVar('spdlTrigIcaReac.mat');

ratList=fieldnames(eachPETH);

col=setCoactColor;
%%
smSigma=100/1000; %in s
cLim=4*[-1,1];


sesIdx=2;
tempName=eachPETH.(ratList{1})(sesIdx).tempName;
tBinSize=eachPETH.(ratList{1})(sesIdx).tBinSize;
tWin=eachPETH.(ratList{1})(sesIdx).param.tWindow;
nWin=ceil(tWin/tBinSize);
tBin=(-nWin:nWin)*tBinSize;
peth=[];
reg={};
for ratIdx=1:length(ratList)
    rat=ratList{ratIdx};
    peth=cat(1,peth,eachPETH.(rat)(sesIdx).pfc.mean);
    reg=[reg,eachPETH.(rat)(sesIdx).region];
end
%%
[reg,regList]=relabel_region(reg);


nPair=sum(ismember(reg,{'vCA1','BLA','PrL L5'}));
eachHeight=(totalHeight-yGap*2)/nPair;

regID=zeros(size(reg));
for regIdx=1:length(regList)
    regID(strcmp(reg,regList{regIdx}))=regIdx;
end
xRange=[-2,2];

periodText={'Pre','Post'};
xTxt={'Time form' 'Spindle peak (s)'};

nSm=ceil(smSigma/tBinSize);
xSm=(-nSm*4:nSm*4)*tBinSize;
smCore=normpdf(xSm,0,smSigma);
smCore=smCore/sum(smCore);


totalY=0;
for n=1:3
    switch n
        case 1
            target='vCA1';
        case 2
            target='BLA';
        case 3
            target='PrL L5';
        otherwise
            continue
    end
    idx=find(strcmp(reg,target));
    
    temp=zscore([peth(idx,:,1)';peth(idx,:,2)']);
    subsetZ{1}=Filter0(smCore,temp(1:nWin*2+1,:));
    subsetZ{2}=Filter0(smCore,temp(nWin*2+1+1:end,:));
    
    peak=mean(subsetZ{2}(nWin+1+(-1:1),:),1);
    [~,order]=sort(peak);
    
    height=length(idx)*eachHeight;
    for prePost=1:2
        
        subplotInMM(x+(width+xGap)*(prePost-1),y+totalY,width,height)
        
        imagescXY(tBin,[],subsetZ{prePost}(:,order));
        box off
        xlim(xRange)
        set(gca,'YTick',[])
        set(gca,'clim',cLim)
        ax=fixAxis;
        if n==3
            xlabel(xTxt)
        else
            set(gca,'XTickLabel',[])
        end
        if n==1
            text2(0.5,1,{sprintf('%s-conditoning',periodText{prePost}),'NREM'},ax,...
                'fontsize',5,'horizontalAlign','center','verticalAlign','bottom')
        end
        if prePost==1
            ylabel(sprintf('%s ensemble',target),'fontsize',5,'fontweight','normal')
        end
        colormap(gca,col.react.map)
    end
    totalY=totalY+height+yGap;
end
subplotInMM(x+width*2+xGap*2,y,1,totalY-yGap)
imagescXY([0,1],cLim,linspace(cLim(1),cLim(2),256));
box off
set(gca,'XTick',[],'YTick',[cLim(1),cLim(1)/2,0,cLim(2)/2,cLim(2)],'YAxisLocation','right')
colormap(gca,col.react.map)
ax=fixAxis;
text2(5,0.5,'Reactivation strength (z)',ax,'fontsize',5,'horizontalALign','center','rotation',-90)
end

%%
function panel_06old(x,y);
xGap=2;
width=((162-10*2)/3-xGap*2-3)/2;
yGap=1.5;
totalHeight=60;
eachPETH=poolVar('hfoTrigIcaReac.mat');

ratList=fieldnames(eachPETH);

col=setCoactColor;
%%
smSigma=20/1000; %in s
cLim=4*[-1,1];


sesIdx=2;
tempName=eachPETH.(ratList{1})(sesIdx).tempName;
tBinSize=eachPETH.(ratList{1})(sesIdx).tBinSize;
tWin=eachPETH.(ratList{1})(sesIdx).param.tWindow;
nWin=ceil(tWin/tBinSize);
tBin=(-nWin:nWin)*tBinSize*1e3;
peth=[];
reg={};
for ratIdx=1:length(ratList)
    rat=ratList{ratIdx};
    peth=cat(1,peth,eachPETH.(rat)(sesIdx).nrem.mean);
    reg=[reg,eachPETH.(rat)(sesIdx).region];
end
%%
[reg,regList]=relabel_region(reg);
regID=zeros(size(reg));
for regIdx=1:length(regList)
    regID(strcmp(reg,regList{regIdx}))=regIdx;
end

nPair=sum(ismember(reg,{'vCA1','BLA','PrL L5'}));
eachHeight=(totalHeight-yGap*2)/nPair;


xRange=[-420,420];

periodText={'Pre','Post'};
xTxt={'Time form' 'HFO peak (ms)'};

nSm=ceil(smSigma/tBinSize);
xSm=(-nSm*4:nSm*4)*tBinSize;
smCore=normpdf(xSm,0,smSigma);
smCore=smCore/sum(smCore);


totalY=0;
for n=1:3
    switch n
        case 1
            target='vCA1';
        case 2
            target='BLA';
        case 3
            target='PrL L5';
        otherwise
            continue
    end
    idx=find(strcmp(reg,target));
    
    temp=zscore([peth(idx,:,1)';peth(idx,:,2)']);
    subsetZ{1}=Filter0(smCore,temp(1:nWin*2+1,:));
    subsetZ{2}=Filter0(smCore,temp(nWin*2+1+1:end,:));
    
    peak=mean(subsetZ{2}(nWin+1+(-1:1),:),1);
    [~,order]=sort(peak);
    
    height=length(idx)*eachHeight;
    for prePost=1:2
        
        subplotInMM(x+(width+xGap)*(prePost-1),y+totalY,width,height)
        
        imagescXY(tBin,[],subsetZ{prePost}(:,order));
        box off
        xlim(xRange)
        set(gca,'YTick',[])
        set(gca,'clim',cLim)
        ax=fixAxis;
        if n==3
            xlabel(xTxt)
        else
            set(gca,'XTickLabel',[])
        end
        if n==1
            text2(0.5,1,{sprintf('%s-conditoning',periodText{prePost}),'NREM'},ax,...
                'fontsize',5,'horizontalAlign','center','verticalAlign','bottom')
        end
        if prePost==1
            ylabel(sprintf('%s ensemble',target),'fontsize',5,'fontweight','normal')
        end
        colormap(gca,col.react.map)
    end
    totalY=totalY+height+yGap;
end
subplotInMM(x+width*2+xGap*2,y,1,totalY-yGap)
imagescXY([0,1],cLim,linspace(cLim(1),cLim(2),256));
box off
set(gca,'XTick',[],'YTick',[cLim(1),cLim(1)/2,0,cLim(2)/2,cLim(2)],'YAxisLocation','right')
colormap(gca,col.react.map)
ax=fixAxis;
text2(5,0.5,'Reactivation strength (z)',ax,'fontsize',5,'horizontalALign','center','rotation',-90)
end
%%
