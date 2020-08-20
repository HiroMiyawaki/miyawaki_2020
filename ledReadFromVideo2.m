function ledReadFromVideo2(videofile,imageRange,tStart,tEnd,sessionName,saveFileName)
    close all
    nBaseFrame=250;%10 sec
    fprintf('%s start %s with data of %s\n',datestr(now),mfilename,videofile)

    if ~exist('tStart','var')|| isempty(tStart)
        tStart=0;
    end    
    
    vr=VideoReader(videofile);

    if ~exist('tEnd','var')|| isempty(tEnd)
        tEnd=vr.Duration;
    end
    
    if ~exist('sessionName','var') || isempty(sessionName)
        sessionName='_';
    else
        sessionName=['_' sessionName '_'];
    end    
    
    
    vr.CurrentTime=tStart; 
    if isscalar(imageRange)
        nRange=imageRange;
        imageRange=zeros(nRange,4);
        
        for idx=1:nBaseFrame
            if vr.hasFrame
                baseVideo(:,:,:,idx)=double(readFrame(vr));
            else
                baseVideo(:,:,:,idx:end)=[];
                break
            end
        end

        baseImg=mean(baseVideo,4)/255;        
        
        
        close all
        fh=figure('position',[200,1000,1280,960]);
        imagesc(baseImg);
        axis equal
        set(gca,'clim',[0,255]);
        colormap(gray);
        title('Select LED area')
        for n=1:nRange
            [y,x]=ginput(2);
            xRange=[floor(min(x)),ceil(max(x))];
            yRange=[floor(min(y)),ceil(max(y))];
            imageRange(n,:)=[xRange,yRange];
            rectangle('position',[yRange(1),xRange(1),diff(yRange),diff(xRange)],'EdgeColor','r')
            drawnow;
        end

    elseif size(imageRange,2)==4 && size(imageRange,1)>0
        nRange=size(imageRange,1);
    else
        displ('Give number of LEDs or list of range as n by 4 matrix')
        return
    end
%%


vr.CurrentTime=tStart; idx=0;
clear xRange yRange
for n=1:nRange
    xRange{n}=imageRange(n,1):imageRange(n,2);
    yRange{n}=imageRange(n,3):imageRange(n,4);
end

nFrame=round((tEnd-tStart)*vr.FrameRate);

progStep=0.05;
prog=progStep;

disp([datestr(now) ' start processing'])

ledVal=nan(nRange,nFrame);

try
    while vr.hasFrame
        idx=idx+1;
        
        if idx/nFrame>prog
            disp([datestr(now) ' ' num2str(prog*100) '% done'])
            prog=prog+progStep;
        end

        fr=vr.readFrame;
        
        if vr.CurrentTime>tEnd
            break;
        end        
        
        for n=1:nRange
            subFr=double(fr(xRange{n},yRange{n},:));
            brightness=0.299*subFr(:,:,1)+0.587*subFr(:,:,2)+0.114*subFr(:,:,3);
            ledVal(n,idx)=mean(brightness(:));
        end
    end
    disp(['whole video was processed (' num2str(idx) ' frames)'])
    
catch
    disp(['process ended at frame ' num2str(idx)])
end
lastFrame=idx;
       
%%



if ~exist('saveFileName','var')    
    ext=findstr(videofile,'.');
    saveFileCore=[videofile(1:ext(end)-1) sessionName 'ledBlink'];
    saveFileName=[saveFileCore '.mat'];
else
    if strcmpi(saveFileName(end-3:end),'.mat')
        saveFileCore=saveFileName(1:end-4);
    else
        saveFileCore=saveFileName;
        saveFileName=[saveFileName '.mat'];
    end
end


if exist(saveFileName,'file')
    
    backupName=[saveFileCore '-backup.mat'];
    idx=0;
    while exist(backupName,'file')
        idx=idx+1;
        backupName=[saveFileCore '-backup' num2str(idx) '.mat'];
    end
    disp(['Previous results was backed up as ' backupName])
    unix(['mv ' saveFileName ' ' backupName]);
end


param.videofile=videofile;
param.imageRange=imageRange;
param.frameRate=vr.FrameRate;

param.tStart=tStart;
param.tEnd=tEnd;

param.madeby=mfilename;

save(saveFileName,'ledVal','param','-v7.3');
disp([datestr(now) 'Results saved in ' saveFileName])
