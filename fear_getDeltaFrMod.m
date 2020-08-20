function fear_getDeltaFrMod(basename,varargin)

load([basename '.basicMetaData.mat'])
fprintf('%s start %s with data of %s\n',datestr(now),mfilename,basicMetaData.SessionName)

%%
load([basicMetaData.Basename '.sessions.events.mat'])
load([basicMetaData.Basename '.sleepstate.states.mat'])
load([basicMetaData.Basename '.okUnit.spikes.mat'])

%%
param.preWin=0.25;
param.postWin=0.25;
param.peakWin=0.03;
%%
param=parseParameters(param,varargin);
%%
cBorder=0.5:max(okUnit.cluster)+0.5;
beh=relabel_ma2sleep(SleepState.MECE.timestamps);

nrem=beh(beh(:,3)==3,1:2);
tBorder=[-inf,sort(nrem(:))',inf];
cnt=histcounts2(okUnit.spikeTime,okUnit.cluster,tBorder,cBorder);
cnt=cnt(2:2:end,:);
base=sum(cnt,1)/sum(diff(nrem,1,2));


%%
for regIdx=1:3 
    switch regIdx
        case 1
            reg='pfc';
        case 2
            reg='amy';
        case 3
            reg='hpc';
    end    
    if ~exist([basicMetaData.Basename '.' reg 'Slowwave.events.mat'],'file')
        warning('Slowwaves in %s were not detected on %s',reg,basicMetaData.SessionName)
        continue
    end    
    
    
    temp=load([basicMetaData.Basename '.' reg 'Slowwave.events.mat']);
    varName=fieldnames(temp);
    slowwave=temp.(varName{1});
    
    % delta spikes

    deltaPeak=slowwave.peak.timestamps+ param.peakWin/2*[-1,1];
    deltaPeak(~any(deltaPeak(:,1)>nrem(:,1)' & deltaPeak(:,2)<nrem(:,2)',2),:)=[];
    
    tBorder=[-inf,sort(deltaPeak(:))',inf];
    cnt=histcounts2(okUnit.spikeTime,okUnit.cluster,tBorder,cBorder);
    deltaSpk=cnt(2:2:end,:);
    
    pre=slowwave.timestamps(:,1)+[-param.preWin,0];    
    pre(~any(pre(:,1)>nrem(:,1)' & pre(:,2)<nrem(:,2)',2),:)=[];
    while any(pre(1:end-1,2)>pre(2:end,1))
        pre(find(pre(1:end-1,2)>pre(2:end,1))+1,:)=[];
    end
    tBorder=[-inf,sort(pre(:))',inf];
    cnt=histcounts2(okUnit.spikeTime,okUnit.cluster,tBorder,cBorder);
    preSpk=cnt(2:2:end,:);
    preGain=(preSpk/param.preWin)./base;
    
    post=slowwave.timestamps(:,2)+[0,param.postWin];
    post(~any(post(:,1)>nrem(:,1)' & post(:,2)<nrem(:,2)',2),:)=[];
    while any(post(1:end-1,2)>post(2:end,1))
        post(find(post(1:end-1,2)>post(2:end,1)),:)=[];
    end
    tBorder=[-inf,sort(post(:))',inf];
    cnt=histcounts2(okUnit.spikeTime,okUnit.cluster,tBorder,cBorder);
    postSpk=cnt(2:2:end,:);
    postGain=(postSpk/param.postWin)./base;
    
        

    
    %%
    for mesType=1:3
        switch mesType
            case 1
                mesName='deltaSpike';
                val=deltaPeak;
                gain=deltaSpk;
            case 2
                mesName='preDelta';
                val=pre;
                gain=preGain;
            case 3
                mesName='postDelta';
                val=post;
                gain=postGain;
        end
                
        for idx=1:5
            toUse=(val(:,1)>sessions.homecage(idx,1) & val(:,2) <sessions.homecage(idx,2));

            deltaFrMod.(reg).(mesName)(idx).participation=mean(gain(toUse,:)>0,1);
            if mesType==1
                deltaFrMod.(reg).(mesName)(idx).nSpk=mean(gain(toUse,:),1);
            else
                deltaFrMod.(reg).(mesName)(idx).gain=mean(gain(toUse,:),1);
            end
            deltaFrMod.(reg).(mesName)(idx).n=sum(toUse);

        end
    end
end
deltaFrMod.baseFR=base;
deltaFrMod.param=param;
deltaFrMod.generator=mfilename;
deltaFrMod.generatedate=datestr(now,'yyyy-mm-dd');

%%
save([basicMetaData.AnalysesName '-deltaFrMod.mat'],'deltaFrMod','-v7.3')