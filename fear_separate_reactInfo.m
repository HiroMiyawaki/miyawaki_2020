function fear_separate_reactInfo(basename)
load([basename '.basicMetaData.mat'])
fprintf('  %s start %s with data of %s\n',datestr(now),mfilename,basicMetaData.SessionName)


load([basicMetaData.AnalysesName '-instReac20ms.mat'])

instReacInfo=rmfield(instReac20ms,'strength');
for n=1:length(instReacInfo)
    instReacInfo(n).generatedate=datestr(now,'yyyy-mm-dd');
    instReacInfo(n).generator=mfilename;
end
clear instReac20ms

load([basicMetaData.AnalysesName '-icaReac.mat'])
icaReacInfo=rmfield(icaReac,'strength');
for n=1:length(icaReacInfo)
    icaReacInfo(n).generatedate=datestr(now,'yyyy-mm-dd');
    icaReacInfo(n).generator=mfilename;
end

clear icaReac

save([basicMetaData.AnalysesName '-instReacInfo.mat'],'instReacInfo','-v7.3')
save([basicMetaData.AnalysesName '-icaReacInfo.mat'],'icaReacInfo','-v7.3')
