function varargout = fear_FR_freezeModulation_KS2(basename,varargin)

%% load basic meta data
load([basename '.basicMetaData.mat']);
fprintf('%s start %s with data of %s\n',datestr(now),mfilename,basicMetaData.SessionName)
%% set default
param.minWakeDur=40; %in sec
param.minIsoDist=15;
param.maxIsiIdx=0.1;
param.minNcell=3;
param.minFR=1e-4;
param.varName='frzFR';
param.showPlot=true;
param.figFileName=[basicMetaData.AnalysesPdf '-frzFR.pdf'];
param.saveFileName=fullfile(basicMetaData.BaseDir,'analyses',[basicMetaData.SessionName '-frzFR.mat']);
%% set options
param=parseParameters(param,varargin);

%% load files

load([basicMetaData.Basename '.okUnit.spikes.mat'])
load([basicMetaData.Basename '.freezeHMM.events.mat'])
load([basicMetaData.Basename '.sessions.events.mat'])
load([basicMetaData.Basename '.cues.events.mat'])
load([basicMetaData.Basename '.sleepState.states.mat'])

%% get sleep periods (including MA)
slp=SleepState.MECE.timestamps(...
    SleepState.MECE.timestamps(:,3)>1 | ...
    (SleepState.MECE.timestamps(:,3)==1 & diff(SleepState.MECE.timestamps(:,1:2),1,2)<param.minWakeDur),...
    1:2);

slp=removeTransient(slp,eps,eps);
%%
temp=mergePeriod(freezeHMM.timestamps,slp,basicMetaData.detectionintervals.lfp(:));

if ~isempty(temp{2,2})
    warning('there is overlap between freeze and sleep: overlap ignored')
end

beh=sortrows([
    temp{2,1},1*ones(size(temp{2,1},1),1) %freeze
    temp{1,1},2*ones(size(temp{1,1},1),1) %non freeze wake
    temp{1,2},3*ones(size(temp{1,2},1),1) %sleep
    temp{2,2},nan*ones(size(temp{2,2},1),1) %error: freeze during sleep
    ]);

%% get period in the chambers

sesTime=sessions.timestamps;
% split Cue and Extinction in thirds
extIdx=find(strcmpi({basicMetaData.chamber.name},'CueAndExtinction'));
FirstPip=min(cues.timestamps.Pip(cues.timestamps.Pip(:,1)>=sesTime(extIdx,1),1));

temp=cues.timestamps.Tone(cues.timestamps.Tone(:,1)>=sesTime(extIdx,1),1);
extPipStart=min(cues.timestamps.Pip(cues.timestamps.Pip(:,1)>=temp(9)));

borders=[sesTime(extIdx,1),FirstPip,extPipStart,sesTime(extIdx,2)];

sesTime=[sesTime(1:extIdx-1,:)
    borders(1:end-1)',borders(2:end)'
    sesTime(extIdx+1:end,:)];

sesNameList={basicMetaData.chamber(1:extIdx-1).name};

sesNameList{end+1}=[basicMetaData.chamber(extIdx).name '-Base'];
sesNameList{end+1}=[basicMetaData.chamber(extIdx).name '-Cue'];
sesNameList{end+1}=[basicMetaData.chamber(extIdx).name '-Ext'];
sesNameList={sesNameList{:},basicMetaData.chamber(extIdx+1:end).name};
%% get cell list to be used
okCells=1:length(okUnit.cluInfo.channel);
%% get region list to be used
regList=unique(okUnit.cluInfo.region);
nCell=sum(strcmp(repmat(okUnit.cluInfo.region,length(regList),1),repmat(regList',1,length(okCells))),2)';

regList=regList(nCell>=param.minNcell);

%%
for sesIdx=1:length(sesNameList);
    sesName=sesNameList{sesIdx};
    tEdges=beh(beh(:,2)>sesTime(sesIdx,1)&beh(:,1)<sesTime(sesIdx,2),[1,3]);
    
    if any(diff(tEdges(:,1))<=0)
        error('wrong freeze subset');
    end
    
    if ~isempty(tEdges)
        if tEdges(1) >sesTime(sesIdx,1)
            error('first tEdges seems wrong in session %d of %s',sesIdx,basicMetaData.SessionName)
            tEdges=[sesTime(sesIdx,1);tEdges];
        else
            tEdges(1,1)=sesTime(sesIdx,1);
        end
        
        if tEdges(end) <sesTime(sesIdx,2)
            tEdges=[tEdges;sesTime(sesIdx,2),nan;];
        else
            error('last tEdges seems wrong in session %d of %s',sesIdx,basicMetaData.SessionName)
        end
    else
        error('no tEdges in session %d of %s',sesIdx,basicMetaData.SessionName)
    end
    
            
    subSpk=okUnit.spikeTime(okUnit.spikeTime>=sesTime(sesIdx,1) & okUnit.spikeTime<=sesTime(sesIdx,2));
    subClu=okUnit.cluster(okUnit.spikeTime>=sesTime(sesIdx,1) & okUnit.spikeTime<=sesTime(sesIdx,2));
    for regIdx=1:length(regList)        
%         regName=strrep(regList{regIdx},' ','_');
        regName=regList{regIdx};
        
        cList=sort(find(strcmpi(okUnit.cluInfo.region,regName)))';
        
        
        subT=subSpk(ismember(subClu,cList));
        subC=subClu(ismember(subClu,cList));
        
        cEdges=[cList-0.5;cList(end)+0.5];
        
        cnt=hist3([subT,subC],'edges',{tEdges(:,1),cEdges});
        
        cnt(:,end)=[];
        cnt(end,:)=[];
        
        dur=diff(tEdges(:,1));
        
        bType=tEdges(1:end-1,2);

        for bTypeIdx=1:3
            frzFR(regIdx).fr(bTypeIdx,:,sesIdx)=sum(cnt(bType==bTypeIdx,:),1)/sum(dur(bType==bTypeIdx));
            frzFR(regIdx).each(sesIdx).cnt=cnt;
            frzFR(regIdx).each(sesIdx).dur=dur;            
        end
        frzFR(regIdx).regName=regName;
        frzFR(regIdx).sesName=sesNameList;
        frzFR(regIdx).behType={'freeze','non-freeze','sleep'};
        frzFR(regIdx).generator=mfilename;
        frzFR(regIdx).generatedate=today('datetime');
        frzFR(regIdx).generateParam=param;
    end
end
%% plot results
if param.showPlot || (isempty(param.figFileName) && isempty(param.saveFileName))
    nRow=max(arrayfun(@(x) size(x.fr,3),frzFR))+1;

    nCol=length(frzFR);
%     close all
    fh=initFig4A4('fontsize',7);
    col=[0,0,1;1,0,0];
    for regIdx=1:length(frzFR)
        temp=frzFR(regIdx).fr(1:2,:,:);
        temp(temp<param.minFR)=param.minFR;
        frRange=log([min(temp(:)),    max(temp(:))]);        
        frRange= exp(frRange+diff(frRange)*0.05*[-1,1]);

        for sesIdx=1:size(frzFR(regIdx).fr,3)
            subplot2(nRow,nCol,sesIdx,regIdx)
            xlim(frRange)
            ylim(frRange);
            plotIdentityLine(gca,{'color',0.5*[1,1,1]})
            hold on
            scatter((temp(1,:,sesIdx)),(temp(2,:,sesIdx)),1,col(2,:))
            set(gca,'XScale','log','YScale','log')
            xlabel(['FR_{' frzFR(regIdx).behType{1} '} (Hz)'])
            ylabel(['FR_{' frzFR(regIdx).behType{2} '} (Hz)'])
            ax=fixAxis;
            if sesIdx==1
                title(frzFR(regIdx).regName)
            end
            if regIdx==1 
                text2(-0.125*nCol,0.5,frzFR(regIdx).sesName{sesIdx},ax,'horizontalAlign','center','rotation',90,'fontsize',7)
            end
        end
    end
    addScriptName(mfilename)
    subplot2(nRow,nCol,nRow,1)
    xlim([0,1])
    ylim([0,1])
    txt{1}='Only good cells were used';
    txt{2}=sprintf('FRs < %0.1e are plotted as %0.1e',param.minFR,param.minFR);
    txt{3}='Sleep periods are excluded from the analyses';
    ax=fixAxis;
    text2(0,0.5,txt,ax,'verticalAlign','top')
    axis off
    drawnow
    textInMM(15,10,basicMetaData.SessionName,'fontsize',18,'horizontalAlign','left')
    
    if ~isempty(param.figFileName)
        drawnow
        print(fh,param.figFileName,'-dpdf','-painters','-r300')
    end
end
%% save results
if ~isempty(param.varName) &&  ~isempty(param.saveFileName)
    if strcmp(param.varName,'frzFR')
        eval([param.varName '=frzFR;']);
    end
    save(param.saveFileName,param.varName,'-v7.3')
end
%% set output
if nargout>0
    varargout{1}=frzFR;
end

    
