function fear_tonePETHeach(basename,varargin)
load([basename '.basicMetadata.mat'])
fprintf('%s start %s with data of %s\n',datestr(now),mfilename,basicMetaData.SessionName)

load([basicMetaData.Basename '.okUnit.spikes.mat'])

load([basicMetaData.Basename '.sessions.events.mat'])
load([basicMetaData.Basename '.cues.events.mat'])

load([basicMetaData.Basename '.sleepstate.states.mat'])
%%
param.smoothSigma=0.5;
param.tBinSize=0.1;
param.tRange=[-20,80];
%%
param=parseParameters(param,varargin);

%%
smSigma=param.smoothSigma;
tBinSize=param.tBinSize;
tRange=param.tRange;

%%
cluBin=unique(okUnit.cluster);
cluBorder=[cluBin-0.5;max(cluBin)+0.5];

nSmBin=ceil(smSigma*4/tBinSize);

tBorder=tRange(1)-tBinSize/2-tBinSize*nSmBin:tBinSize:tRange(2)+tBinSize/2+tBinSize*nSmBin;
tBin=(tBorder(1:end-1)+tBorder(2:end))/2;

fr=zeros(length(tBin)-2*nSmBin,length(cluBin),size(cues.timestamps.Tone,1));
smFR=fr;
%%
smX=0:tBinSize:smSigma*4;
smX=[-fliplr(smX),smX(2:end)];
smCore=normpdf(smX,0,smSigma);
smCore=smCore/sum(smCore);

%%
nTone=zeros(1,size(sessions.timestamps,1));
for sesIdx=1:size(sessions.timestamps,1)
    tones=cues.timestamps.Tone(cues.timestamps.Tone(:,2)<sessions.timestamps(sesIdx,2) & cues.timestamps.Tone(:,2)>sessions.timestamps(sesIdx,1),:);
    onset=[];
    for tIdx=1:size(tones,1)
        onset(tIdx)=cues.timestamps.Pip(find(cues.timestamps.Pip(:,1)>=tones(tIdx,1),1,'first'),1);
    end
    

    spk=okUnit.spikeTime(okUnit.spikeTime>sessions.timestamps(sesIdx,1) & okUnit.spikeTime<sessions.timestamps(sesIdx,2));
    clu=okUnit.cluster(okUnit.spikeTime>sessions.timestamps(sesIdx,1) & okUnit.spikeTime<sessions.timestamps(sesIdx,2));
    for tIdx=1:length(onset)    
        temp=histcounts2(spk,clu,tBorder+onset(tIdx),cluBorder)/tBinSize;;
        fr(:,:,sum(nTone)+tIdx)=temp(nSmBin+1:end-nSmBin,:);
        temp=Filter0(smCore,temp);
        smFR(:,:,sum(nTone)+tIdx)=temp(nSmBin+1:end-nSmBin,:);
    end
    nTone(sesIdx)=length(onset);
end
%%
tonePETHeach.fr=fr;
tonePETHeach.smoothed=smFR;
tonePETHeach.tBin=tBin(nSmBin+1:end-nSmBin);
tonePETHeach.nTonePerSession=nTone;
tonePETHeach.param=param;
tonePETHeach.generator=mfilename;
tonePETHeach.generatedate=datestr(now,'yyyy-mm-dd');
%%
save([basicMetaData.AnalysesName '-tonePETHeach.mat'],'tonePETHeach','-v7.3')
