function fear_evtTrigReact(basename,varargin)

load([basename '.basicMetaData.mat'])
fprintf('%s start %s with data of %s\n',datestr(now),mfilename,basicMetaData.SessionName)

%%
param.binSize=0.02; % in sec
param.halfWindow=5; % in sec
param.targetHC=3;
param.behavior='nrem'; %nrem,rem,wake or entire

param.varName='evtTrigReact';
param.saveFile=[basicMetaData.AnalysesName '-evtTrigReact.mat'];

param.reacFile=[basicMetaData.AnalysesName '-icaCoactTimeCondHT.mat'];
%%
param=parseParameters(param,varargin);

%%
if exist([basicMetaData.Basename '.ripples.events.mat'],'file')
    load([basicMetaData.Basename '.ripples.events.mat'])
else
    ripples.timestamps=zeros(0,2);
    ripples.peaks.timestamps=zeros(0,1);
end
load([basicMetaData.Basename '.pfcSlowwave.events.mat'])
load([basicMetaData.Basename '.pfcSpindle.events.mat'])
load([basicMetaData.Basename '.SleepState.states.mat'])
load([basicMetaData.Basename '.sessions.events.mat'])

temp=load(param.reacFile);
vName=fieldnames(temp);
coactTime=temp.(vName{1});

%%
trig{1}=ripples.peaks.timestamps;
trig{2}=ripples.timestamps(:,1);
trig{3}=ripples.timestamps(:,2);
trig{4}=pfcSpindle.peaktime';
trig{5}=pfcSpindle.timestamps(:,1);
trig{6}=pfcSpindle.timestamps(:,2);
trig{7}=pfcSlowWave.peak.timestamps;
trig{8}=pfcSlowWave.timestamps(:,1);
trig{9}=pfcSlowWave.timestamps(:,2);

trigName={'SWR peak','SWR onset','SWR offset','Spindle peak','Spindle onset','Spindle offset','DOWN peak','DOWN onset','DOWN offset'};

%%
slp=relabel_ma2sleep(SleepState.MECE.timestamps);

%%
hcIdx=param.targetHC;
tRange=sessions.homecage(hcIdx,:);

switch lower(param.behavior)
    case 'nrem'
        beh=slp(slp(:,3)==3,1:2);
    case 'rem'
        beh=slp(slp(:,3)==5,1:2);
    case 'wake'
        beh=slp(slp(:,3)==1,1:2);
    case 'entire'
        beh=inf*[-1,1];
    otherwise
        error('behavior must be nrem,rem,wake or entire')        
end
    
beh=beh(beh(:,2)>tRange(1) & beh(:,1)<tRange(2),:);
if beh(1)<tRange(1); beh(1)=tRange(1); end
if beh(end)>tRange(2); beh(end)=tRange(2); end

%%
for n=1:length(trig)
    temp=trig{n}(trig{n}>tRange(1) & trig{n}<tRange(2));
    trig{n}=temp(any(temp>beh(:,1)' & temp<beh(:,2)',2));
end

%%
target=find(coactTime.sigLevel>0);

evt={};
for n=1:length(target)
    temp=coactTime.timestamp{target(n)};
    temp=temp(temp>tRange(1) & temp<tRange(2))';
    evt{n}=temp(any(temp>beh(:,1)' & temp<beh(:,2)',2));
end
%%
t=cat(1,trig{:},evt{:});

nG=cumsum([cellfun(@length,trig),cellfun(@length,evt)]);
g=ones(1,nG(end));
for n=1:length(nG)-1
    g(nG(n)+1:nG(end))=g(nG(n)+1:nG(end))+1;
end
%%
nWin=ceil(param.halfWindow/param.binSize);
[cnt,time]=CCG(t,g,param.binSize,nWin,1);

rate=zeros(nWin*2+1,length(trig),length(evt));
for n=1:length(trig)
    rate(:,n,:)=cnt(:,n,length(trig)+1:end)/length(trig{n})/param.binSize;
end
%%

evtTrigReact.rate=rate;
evtTrigReact.triger.n=cellfun(@length,trig);
evtTrigReact.triger.name=trigName;
evtTrigReact.pairID=target;
evtTrigReact.reacID=coactTime.reacID(target,:);
evtTrigReact.region=coactTime.region(target,:);
evtTrigReact.sigLevel=coactTime.sigLevel(target);
evtTrigReact.tGap=coactTime.tGap(target);

evtTrigReact.param=param;
evtTrigReact.generator=mfilename;
evtTrigReact.generatedate=datestr(now,'yyyy-mm-dd');

if ~strcmp(param.varName,'evtTrigReact')
    eval(sprintf('%s=evtTrigReact;',param.varName))
end
save(param.saveFile,param.varName,'-v7.3');


