function fear_coactTrigWaveletHighThreshold(basename,varargin)
load([basename '.basicMetaData.mat'])
fprintf('%s start %s with data of %s\n',datestr(now),mfilename,basicMetaData.SessionName);

param.lfp=basicMetaData.lfp;
param.Ch=[basicMetaData.Ch.hpcTheta,basicMetaData.Ch.pfcDelta,basicMetaData.Ch.amyGamma];
param.templateSession=2;
param.targetHomecage=3;
param.tWin=1.1;
%%
param=parseParameters(param,varargin);
%%
load([basicMetaData.AnalysesName '-icaCoactTimeHT.mat'])
load([basicMetaData.AnalysesName '-instCoactTimeHT.mat'])
load([basicMetaData.Basename '.sleepState.states.mat'])
load([basicMetaData.Basename '.sessions.events.mat'])

%%
behList={};
noTarget=true;

tempSes=param.templateSession;
targetHC=param.targetHomecage;

temp=fieldnames(icaCoactTime(tempSes,targetHC));

for idx=1:length(temp)
    if isfield(icaCoactTime(tempSes,targetHC).(temp{idx}),'timestamp')||...
            isfield(instCoactTime(tempSes,targetHC).(temp{idx}),'timestamp')
        noTarget=false;
        behList{end+1}=temp{idx};
    elseif isfield(icaCoactTime(tempSes,targetHC).(temp{idx}),'pairID')||...
            isfield(instCoactTime(tempSes,targetHC).(temp{idx}),'pairID')
        behList{end+1}=temp{idx};
    end
end

%%

tMin=sessions.homecage(param.targetHomecage,1);
tMax=sessions.homecage(param.targetHomecage,2);

temp=fieldnames(instCoactTime(tempSes,targetHC));
for idx=1:length(temp)
    if isfield(instCoactTime(tempSes,targetHC).(temp{idx}),'timestamp')
        tTemp=[instCoactTime(tempSes,targetHC).(temp{idx}).timestamp{:}];
        tMin=min([tTemp,tMin]);
        tMax=max([tTemp,tMax]);
    end
end

temp=fieldnames(icaCoactTime(tempSes,targetHC));
for idx=1:length(temp)
    if isfield(icaCoactTime(tempSes,targetHC).(temp{idx}),'timestamp')
        tTemp=[icaCoactTime(tempSes,targetHC).(temp{idx}).timestamp{:}];
        tMin=min([tTemp,tMin]);
        tMax=max([tTemp,tMax]);
    end
end


tRange=[max(tMin-10,basicMetaData.detectionintervals.lfp(1)),min(tMax+10,basicMetaData.detectionintervals.lfp(2))];
fRange=[floor(tRange(1)*basicMetaData.SampleRates.lfp),ceil(tRange(2)*basicMetaData.SampleRates.lfp)];

%%

slp=relabel_ma2sleep(SleepState.MECE.timestamps);

%%
tBin=(fRange(1):fRange(2))/basicMetaData.SampleRates.lfp;
inNrem=false(size(tBin));
for idx=1:size(slp,1)
    if slp(idx,3)==3
        inNrem(tBin>slp(idx,1)&tBin<slp(idx,2))=true;
    end
end

%%

nFrame=diff(fRange)+1;
nBuff=100000;
nChank= floor(nFrame/nBuff);

fh=fopen(param.lfp);
lfp=zeros(length(param.Ch),nFrame);
fseek(fh,basicMetaData.nCh*(fRange(1)-1)*2,'bof');
for n=0:nChank-1
    if mod(n,10)==0
        fprintf('%s loading LFP of %s, %d/%d\n',datestr(now),basicMetaData.SessionName,n,nChank);
    end
    temp=fread(fh,[basicMetaData.nCh,nBuff],'int16');
    lfp(:,n*nBuff+(1:nBuff))=temp(param.Ch,:);
end

if mod(nFrame,nBuff)~=0
    temp=fread(fh,[basicMetaData.nCh, mod(nFrame,nBuff)],'int16');
    lfp(:,nChank*nBuff+1:end)=temp(param.Ch,:);
end

fclose(fh);
fprintf('%s finish loading LFP of %s\n',datestr(now),basicMetaData.SessionName);

%%
K0=6;
DJ=0.1;
FreqRange=[0.5,330];

fourier_factor=(4*pi)/(K0 + sqrt(2 + K0^2));
scaleMax=(1/FreqRange(1))/fourier_factor;
scaleMin=(1/FreqRange(2))/fourier_factor;

J1=ceil(log2(scaleMax/scaleMin)/DJ);

for idx=1:size(lfp,1);
    [temp,period] = wavelet(lfp(idx,:),1/basicMetaData.SampleRates.lfp,0,DJ,scaleMin,J1,'MORLET',K0);
    pow=abs(temp);
    zPow(:,:,idx)=(pow-mean(pow(:,inNrem),2))./std(pow(:,inNrem),[],2);
end
frequency=arrayfun(@(x) 1/x, period);
%%
tBin=0:ceil(param.tWin*basicMetaData.SampleRates.lfp);
tBin=[-fliplr(tBin),tBin(2:end)];

for reactType=1:2
    clear res
    for bIdx=1:length(behList)
        beh=behList{bIdx};
        fprintf('  %s processing %s (%d/%d) of %s\n',datestr(now),beh,bIdx,length(behList),basicMetaData.SessionName)
        if reactType==1
            coact=icaCoactTime(tempSes,targetHC).(beh);
            varName='icaCoactTrigWavelet';
        else
            coact=instCoactTime(tempSes,targetHC).(beh);
            varName='instCoactTrigWavelet';
        end
        
        across=find(cellfun(@(x,y) ~strcmp(x,y),coact.region(:,1),coact.region(:,2)));
        avg=zeros(size(zPow,1),length(tBin),length(across),size(zPow,3));
        
        for idx=1:length(across)
            targetID=across(idx);
            
            frame=round(coact.timestamp{targetID,1}*basicMetaData.SampleRates.lfp)-fRange(1)+1;
            
            frame(frame+tBin(1)<0)=[];
            frame(frame+tBin(end)>nFrame)=[];
            
            for n=1:length(param.Ch)
                temp=zPow(:,frame+tBin',n);
                temp=reshape(temp,size(zPow,1),length(tBin),length(frame));
                
                
                avg(:,:,idx,n)=squeeze(mean(temp,3));
            end
        end
        res.(beh).wavelet=avg;
        if isfield(coact,'timestamp')
            res.(beh).n=cellfun(@length,coact.timestamp(across,1));
        else
            res.(beh).n=[];
        end
        res.(beh).region=coact.region(across,:);
        res.(beh).pairID =coact.pairID(across);
        res.(beh).reacID =coact.reacID(across);
        res.(beh).tGap=coact.tGap(across);
        res.(beh).sigLevel=coact.sigLevel(across);
    end
    res.t=tBin;
    res.f=frequency;
    res.param=param;
    res.generator=mfilename;
    res.generatedate=datestr(now,'yyyy-mm-dd');

    eval(sprintf('%s=res;',varName));
    
    save([basicMetaData.AnalysesName '-' varName 'HT.mat'],varName,'-v7.3')
end
