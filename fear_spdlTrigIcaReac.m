function fear_spdlTrigIcaReac(basename,varargin)
load([basename '.basicMetaData.mat']);
fprintf('%s start %s with data of %s\n',datestr(now),mfilename,basicMetaData.SessionName)

%%
param.tWindow=2;
param.binSize=20; %20 or 100
param.varName='spdlTrigInstReac';
param.filename=[basicMetaData.AnalysesName '-spdlTrigIcaReac.mat'];

param=parseParameters(param,varargin);
%%
fprintf('start %s\n%s loading data for %s \n',mfilename, datestr(now), basicMetaData.SessionName)

load([basicMetaData.Basename '.Sleepstate.states.mat'])
load([basicMetaData.Basename '.sessions.events.mat'])
load([basicMetaData.Basename '.pfcSpindle.events.mat'])
load([basicMetaData.Basename '.amySpindle.events.mat'])
load([basicMetaData.Basename '.hpcSpindle.events.mat'])
if param.binSize==20
    load([basicMetaData.AnalysesName '-icaReac.mat'])
else
    wornign('bin size for ica reactivation should be 20 ms')
    return;
end

fprintf('%s data have been loaded\n',datestr(now))


%%
spdlPeak{1}=pfcSpindle.troughtime;
spdlPeak{2}=amySpindle.troughtime;
spdlPeak{3}=hpcSpindle.troughtime;
spdlName={'pfc','amy','hpc'};

tBinSize=icaReac(1).param.tBinSize;
tWindow=param.tWindow;
nWin=ceil(tWindow/tBinSize);
nBinTotal=size(icaReac(1).strength,2);
%%
preHC=[1,2,3,3,3,3,4];

for tempSes=1:size(icaReac,2);
    fprintf('%s processing %d/%d templates\n',datestr(now),tempSes,size(icaReac,2))
     for spdlTypeIdx=1:length(spdlName)
        spdlType=spdlName{spdlTypeIdx};
         spdlTrigIcaReac(tempSes).(spdlType)=struct();
     end
    for spdlTypeIdx=1:length(spdlName)
        temp=zeros(size(icaReac(tempSes).strength,1),2*nWin+1,2);
        spdlType=spdlName{spdlTypeIdx};
        numSpdl=[0,0];
        
        for prePost=1:2
            hcIdx=preHC(tempSes)+prePost-1;
            
            
            trig=spdlPeak{spdlTypeIdx}(spdlPeak{spdlTypeIdx}>sessions.homecage(hcIdx,1)&spdlPeak{spdlTypeIdx}<sessions.homecage(hcIdx,2));
            
            trigBin=ceil(trig/tBinSize);
            
            trigBin(trigBin<nWin)=[];
            trigBin(trigBin>nBinTotal-nWin)=[];
            
            
            each=reshape(icaReac(tempSes).strength(:,trigBin+(-nWin:nWin)'),[],2*nWin+1,length(trigBin));
            
            temp(:,:,prePost)=mean(each,3);
            numSpdl(prePost)=length(trigBin);
        end
        spdlTrigIcaReac(tempSes).(spdlType).mean=temp;
        spdlTrigIcaReac(tempSes).(spdlType).n=numSpdl;
%         spdlTrigIcaReac(tempSes).phi=icaReac(tempSes).phi;
        spdlTrigIcaReac(tempSes).region=icaReac(tempSes).region;
        spdlTrigIcaReac(tempSes).tempTime=icaReac(tempSes).tempTime;
        spdlTrigIcaReac(tempSes).tempName=icaReac(tempSes).tempName;
        spdlTrigIcaReac(tempSes).tBinSize=tBinSize;
        spdlTrigIcaReac(tempSes).param=param;
        spdlTrigIcaReac(tempSes).generator=mfilename;
        spdlTrigIcaReac(tempSes).generatedate=datestr(now,'yyyy-mm-dd');
    end
end
%%
fprintf('%s saving data\n',datestr(now,'yyyy-mm-dd'))

if ~strcmp(param.varName,'spdlTrigIcaReac');
    eval(sprintf('%s=spdlTrigIcaReac;',param.varName));
end

save(param.filename,param.varName,'-v7.3')











