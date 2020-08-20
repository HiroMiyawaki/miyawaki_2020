function fear_getfrzFRmod(basename,varargin)
%% load basic meta data
load([basename '.basicMetaData.mat']);
fprintf('%s start %s with data of %s\n',datestr(now),mfilename,basicMetaData.SessionName)

%% load files

load([basicMetaData.Basename '.okUnit.spikes.mat'])
load([basicMetaData.Basename '.freezeHMM.events.mat'])
load([basicMetaData.Basename '.sessions.events.mat'])
load([basicMetaData.Basename '.cues.events.mat'])
load([basicMetaData.Basename '.sleepState.states.mat'])

%% get sleep periods (including MA)
slp=relabel_ma2sleep(SleepState.MECE.timestamps);

slp=removeTransient(slp(slp(:,3)>2,1:2),eps,eps);

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
    sesTime(extIdx+1:end,:)
    borders([2,end])];

sesNameList={basicMetaData.chamber(1:extIdx-1).name};

sesNameList{end+1}=[basicMetaData.chamber(extIdx).name '-Base'];
sesNameList{end+1}=[basicMetaData.chamber(extIdx).name '-Cue'];
sesNameList{end+1}=[basicMetaData.chamber(extIdx).name '-Ext'];
sesNameList={sesNameList{:},basicMetaData.chamber(extIdx+1:end).name};

sesNameList{end+1}=[basicMetaData.chamber(extIdx).name '-CueWhole'];

%%
cList=unique(okUnit.cluster);
cEdges=[cList-0.5;cList(end)+0.5];

mi=@(x) diff(-x,1,1)./sum(x,1);
clear frzFR
frzFRmod.modIdx=zeros(length(cList),length(sesNameList));
frzFRmod.fr=zeros(4,length(cList),length(sesNameList));
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
        
    cnt=hist3([subSpk,subClu],'edges',{tEdges(:,1),cEdges});

    cnt(:,end)=[];
    cnt(end,:)=[];
        
    dur=diff(tEdges(:,1));
        
    bType=tEdges(1:end-1,2);

    for bTypeIdx=1:3
        frzFRmod.fr(bTypeIdx,:,sesIdx)=sum(cnt(bType==bTypeIdx,:),1)/sum(dur(bType==bTypeIdx));
    end
    frzFRmod.fr(4,:,sesIdx)=sum(cnt,1)/sum(dur);
    
    frzFRmod.modIdx(:,sesIdx)=mi(squeeze(frzFRmod.fr(1:2,:,sesIdx)));
end

frzFRmod.behType={'freeze','non-freeze','sleep','all'};
frzFRmod.session=sesNameList;
frzFRmod.generator=mfilename;
frzFRmod.generatedate=today('datetime');

%%

save([basicMetaData.AnalysesName '-frzFRmod.mat'],'frzFRmod','-v7.3')



