function fear_getSwrTrigHist_SubdivideHomecage(basename,varargin)

load([basename '.basicMetaData.mat'])
fprintf('%s start %s with data of %s\n',datestr(now),mfilename,basicMetaData.SessionName)

load([basicMetaData.Basename '.ripples.events.mat'])
load([basicMetaData.Basename '.sessions.events.mat'])
load([basicMetaData.Basename '.SleepState.states.mat'])
load([basicMetaData.Basename '.okUnit.spikes.mat'])

%%
param.nHalfBin=500;
param.tBinSize=1e-3;
param.tSmoothSigma=10e-3;

param.varName='swrTrigHistSub';
param.matFileName=[basicMetaData.AnalysesName '-swrTrigHistSub.mat'];
param.nDiv=5;
%%
param=parseParameters(param,varargin);

%%
beh=relabel_ma2sleep(SleepState.MECE.timestamps);
cluList=unique(okUnit.cluster);
nDiv=param.nDiv;

hc=sessions.homecage;

smT=0:param.tBinSize:param.tSmoothSigma*4;
smT=[-fliplr(smT),smT(2:end)];

smCore=normpdf(smT,0,param.tSmoothSigma);
smCore=smCore/sum(smCore);
nFil=length(smCore);
%%
cnt=zeros(size(cluList,1),2*(param.nHalfBin+nFil)+1,nDiv*size(hc,1));
tRange=[];
sesName={};
for hIdx=1:size(hc,1);
    tBorder=hc(hIdx,1)+diff(hc(hIdx,:))/nDiv*[0:nDiv];
    
    for tIdx=1:nDiv
        tRange(end+1,:)=tBorder(tIdx+(0:1));
        sesName{end+1}=sprintf('Homecage%d-%d',hIdx,tIdx);
        
        subBeh=beh(beh(:,1)<tBorder(tIdx+1) & beh(:,2)>tBorder(tIdx),:);
        subNrem=subBeh(subBeh(:,3)==3,1:2);
        
        subRes=okUnit.spikeTime(okUnit.spikeTime>tBorder(tIdx) & okUnit.spikeTime<tBorder(tIdx+1));
        subClu=okUnit.cluster(okUnit.spikeTime>tBorder(tIdx) & okUnit.spikeTime<tBorder(tIdx+1));
                
        ripT=ripples.peaks.negative(any(ripples.peaks.negative>subNrem(:,1)' & ripples.peaks.negative<subNrem(:,2)',2));
        n(tIdx+nDiv*(hIdx-1))=length(ripT);
        if isempty(ripT)
            cnt(:,:,tIdx+nDiv*(hIdx-1))=nan;
        else
            for cIdx=1:length(cluList)
                sIdx=find(subClu==cluList(cIdx));
                temp=CCG([subRes(sIdx);ripT],[ones(size(sIdx));2*ones(size(ripT))],param.tBinSize,param.nHalfBin+nFil,1);
                cnt(cIdx,:,tIdx+nDiv*(hIdx-1))=squeeze(temp(:,2,1));
            end
        end
    end
end
%%
zs=zscore(cnt,[],2);
smZ=zeros(size(zs));
for idx=1:size(zs,3)
smZ(:,:,idx)=Filter0(smCore,squeeze(zs(:,:,idx))')';
end

%%

swrTrigHistSub.cnt=cnt(:,nFil+1:end-nFil,:);
swrTrigHistSub.zscore=zs(:,nFil+1:end-nFil,:);
swrTrigHistSub.smoothed=smZ(:,nFil+1:end-nFil,:);
swrTrigHistSub.t=(-param.nHalfBin:param.nHalfBin)*param.tBinSize;
swrTrigHistSub.n=n;
swrTrigHistSub.clu=cluList;
swrTrigHistSub.param=param;
swrTrigHistSub.param.sesRange=tRange;
swrTrigHistSub.param.sesName=sesName;

swrTrigHistSub.generator=mfilename;
swrTrigHistSub.generatedate=datestr(now,'yyyy-mm-dd');
%%
if ~strcmpi(param.varName,'swrTrigHistSub')
    eval([param.varName '=swrTrigHistSub;']);
end

save(param.matFileName,param.varName,'-v7.3')
