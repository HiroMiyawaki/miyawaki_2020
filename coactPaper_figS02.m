function coactPaper_figS02()
labelWeight=true;
labelCase=true;
labelSize=9;
letGapX=6;
letGapY=3;

close all
fh=initFig('height',16);

x=3.5;y=3;
panel_01(x,y);
drawnow();

print(fh,'figS02.pdf','-dpdf','-painters','-r300')

end
%%
function panel_01(x,y)
pngDir='~/Dropbox/ratBrain';

if ~exist(pngDir,'dir')
    mkdir(pngDir)
    atlasFile='~/Box\ Sync/textbook/TheRatBrain.pdf';
    gsPath='/usr/local/bin/gs';
    
    pngName=fullfile(pngDir,'ratBrain%03d.png');
    pageList=34:2:354;
    pageList=strrep(num2str(pageList),' ',',');
    while contains(pageList,',,')
        pageList=strrep(num2str(pageList),',,',',');
    end
    
    com=sprintf('%s -dSAFER -dBATCH -dNOPAUSE -sDEVICE=png256 -o %s -r300 -sPageList=%s %s',gsPath,pngName,pageList,atlasFile);
    
    system(com);
end
%%
% for Paxinos & Watson, The Rat Brain, 6th edition
atlasAP={7.56:-0.48:3.72
    3.24:-0.24:2.52
    2.28:-0.12:1.44
    1.28
    1.20:-0.12:-1.56
    -1.72
    -1.80:-0.12:-2.76
    -2.92
    -3.00:-0.12:-4.20
    -4.36
    -4.44:-0.12:-5.04
    -5.20
    -5.28:-0.12:-13.20
    -13.36
    -13.44:-0.12:-13.68
    -13.76
    -13.92:-0.12:-14.64
    -14.76:-0.24:-15.96};
atlasAP=[atlasAP{:}];

getML=@(x,z) 2250+236*x;
getDV=@(y,z) 70+236*(y-(z>=(5.64+5.16)/2)+(z<(2.76+3.00)/2)-(z<(-0.12-0.24)/2)-(z<(-8.28-8.40)/2)+(z<(-13.68-13.76)/2));

fprintf('%s loading png data\n',datestr(now))
for n=1:length(atlasAP)
    [img(:,:,n),map(:,:,n)]=imread(fullfile(pngDir,sprintf('ratBrain%03d.png',n)));
end
fprintf('%s png data have loaded\n',datestr(now))

%%
prPos=[];
sessionList={
    'achel180320'
    'booyah180430'
    'chimay180612'
    'duvel190505'
    'estrella180808'
    'feuillien180830'
    'guiness181002'
    'hoegaarden181115'
    'innis190601'
    'jever190814'
    'karmeliet190901'
    'leffe200124'
    'maredsous200224'
    'nostrum200304'
    'oberon200325'
    };


%achel180320 pfc
ratIdx=1;
prb=2;
tempAP=4-(0:5)*0.2;
tempML=-0.2*ones(1,6);
tempDV=3.2*ones(1,6);

prPos=[prPos;[tempAP;tempML;tempDV;ratIdx*ones(size(tempAP));prb*ones(size(tempAP))]'];

%achel180320 bla
ratIdx=1;
prb=3;
tempAP=-2.7*ones(1,6);
tempML=3.4+0.18*(0:5);
tempDV=8.85*ones(1,6);

prPos=[prPos;[tempAP;tempML;tempDV;ratIdx*ones(size(tempAP));prb*ones(size(tempAP))]'];

%achel180320 hpc
ratIdx=1;
prb=1;
tempAP=-4.55*ones(1,6);
tempML=4.15+0.15*(0:5)*cos(-14/180*pi);
tempDV=7.45+0.15*(0:5)*sin(-14/180*pi);

prPos=[prPos;[tempAP;tempML;tempDV;ratIdx*ones(size(tempAP));prb*ones(size(tempAP))]'];

%booyah180430 pfc
ratIdx=2;
prb=2;
tempAP=3.7+(0:5)*0.2;
tempML=0.5*ones(1,6);
tempDV=3.5*ones(1,6);

prPos=[prPos;[tempAP;tempML;tempDV;ratIdx*ones(size(tempAP));prb*ones(size(tempAP))]'];

%booyah180430 bla
ratIdx=2;
prb=3;
tempAP=-2.15*ones(1,6);
tempML=4.6+0.19*(0:5);
tempDV=8.8*ones(1,6);

prPos=[prPos;[tempAP;tempML;tempDV;ratIdx*ones(size(tempAP));prb*ones(size(tempAP))]'];

%booyah180430 hpc
ratIdx=2;
prb=1;
tempAP=-4.55*ones(1,6);
tempML=4.15+0.15*(0:5)*cos(-14/180*pi);
tempDV=9.3+0.15*(0:5)*sin(-14/180*pi);

prPos=[prPos;[tempAP;tempML;tempDV;ratIdx*ones(size(tempAP));prb*ones(size(tempAP))]'];

%chimay180612 pfc
ratIdx=3;
prb=2;
tempAP=3.2+(0:5)*0.2;
tempML=0.5*ones(1,6);
tempDV=3*ones(1,6);

prPos=[prPos;[tempAP;tempML;tempDV;ratIdx*ones(size(tempAP));prb*ones(size(tempAP))]'];

%chimay180612 bla
ratIdx=3;
prb=3;
tempAP=-3.2*ones(1,6);
tempML=4.7+0.2*(0:5);
tempDV=8.2*ones(1,6);

prPos=[prPos;[tempAP;tempML;tempDV;ratIdx*ones(size(tempAP));prb*ones(size(tempAP))]'];

%chimay180612 hpc
ratIdx=3;
prb=1;
tempAP=-4.45*ones(1,6);
tempML=4.15+0.15*(0:5)*cos(-14/180*pi);
tempDV=8.45+0.15*(0:5)*sin(-14/180*pi);

prPos=[prPos;[tempAP;tempML;tempDV;ratIdx*ones(size(tempAP));prb*ones(size(tempAP))]'];

%duvel190505 pfc
ratIdx=4;
prb=2;
tempAP=3+(0:5)*0.2;
tempML=0.4*ones(1,6);
tempDV=3.5*ones(1,6);

prPos=[prPos;[tempAP;tempML;tempDV;ratIdx*ones(size(tempAP));prb*ones(size(tempAP))]'];


%duvel190505 bla
ratIdx=4;
prb=3;
tempAP=-2.6*ones(1,6);
tempML=3.8+0.15*(0:5);
tempDV=8.65*ones(1,6);

prPos=[prPos;[tempAP;tempML;tempDV;ratIdx*ones(size(tempAP));prb*ones(size(tempAP))]'];

%duvel190505 hpc
ratIdx=4;
prb=1;
tempAP=-4.75*ones(1,6);
tempML=4.3+0.15*(0:5)*cos(-14/180*pi);
tempDV=8.7+0.15*(0:5)*sin(-14/180*pi);

prPos=[prPos;[tempAP;tempML;tempDV;ratIdx*ones(size(tempAP));prb*ones(size(tempAP))]'];

%estrella180808 pfc
ratIdx=5;
prb=2;
tempAP=3.6+(0:5)*0.2;
tempML=0.75*ones(1,6);
tempDV=3.5*ones(1,6);

prPos=[prPos;[tempAP;tempML;tempDV;ratIdx*ones(size(tempAP));prb*ones(size(tempAP))]'];

%estrella180808 bla
ratIdx=5;
prb=3;
tempAP=-2.1*ones(1,6);
tempML=4.4+0.2*(0:5);
tempDV=8.9*ones(1,6);

prPos=[prPos;[tempAP;tempML;tempDV;ratIdx*ones(size(tempAP));prb*ones(size(tempAP))]'];

%estrella180808 hpc
ratIdx=5;
prb=1;
tempAP=-5.5*ones(1,6);
tempML=4.0+0.15*(0:5)*cos(-14/180*pi);
tempDV=8.75+0.15*(0:5)*sin(-14/180*pi);

prPos=[prPos;[tempAP;tempML;tempDV;ratIdx*ones(size(tempAP));prb*ones(size(tempAP))]'];

%feuillien180830 pfc
ratIdx=6;
prb=2;
tempAP=3.2+(0:5)*0.2;
tempML=0.2*ones(1,6);
tempDV=3.9*ones(1,6);

prPos=[prPos;[tempAP;tempML;tempDV;ratIdx*ones(size(tempAP));prb*ones(size(tempAP))]'];

%feuillien180830 bla
ratIdx=6;
prb=3;
tempAP=-3.0*ones(1,6);
tempML=4.3+0.2*(0:5);
tempDV=8.3*ones(1,6);

prPos=[prPos;[tempAP;tempML;tempDV;ratIdx*ones(size(tempAP));prb*ones(size(tempAP))]'];

%feuillien180830 hpc
ratIdx=6;
prb=1;
tempAP=-5.3*ones(1,6);
tempML=3.75+0.15*(0:5)*cos(-14/180*pi);
tempDV=8.9+0.15*(0:5)*sin(-14/180*pi);

prPos=[prPos;[tempAP;tempML;tempDV;ratIdx*ones(size(tempAP));prb*ones(size(tempAP))]'];

%guiness181002 pfc
ratIdx=7;
prb=2;
tempAP=2.8+(0:5)*0.2;
tempML=0.2*ones(1,6);
tempDV=2.9*ones(1,6);

prPos=[prPos;[tempAP;tempML;tempDV;ratIdx*ones(size(tempAP));prb*ones(size(tempAP))]'];

%guiness181002 bla
ratIdx=7;
prb=3;
tempAP=-2.8*ones(1,6);
tempML=4.15+0.16*(0:5);
tempDV=8.375*ones(1,6);

prPos=[prPos;[tempAP;tempML;tempDV;ratIdx*ones(size(tempAP));prb*ones(size(tempAP))]'];

%guiness181002 hpc
ratIdx=7;
prb=1;
tempAP=-5.4*ones(1,6);
tempML=3.55+0.15*(0:5)*cos(-14/180*pi);
tempDV=8.95+0.15*(0:5)*sin(-14/180*pi);

prPos=[prPos;[tempAP;tempML;tempDV;ratIdx*ones(size(tempAP));prb*ones(size(tempAP))]'];


%hoegaarden181115 pfc
ratIdx=8;
prb=2;
tempAP=3.2+(0:5)*0.2;
tempML=0.8*ones(1,6);
tempDV=2.95*ones(1,6);

prPos=[prPos;[tempAP;tempML;tempDV;ratIdx*ones(size(tempAP));prb*ones(size(tempAP))]'];

%hoegaarden181115 bla
ratIdx=8;
prb=3;
tempAP=-2.6*ones(1,6);
tempML=4.45+0.18*(0:5);
tempDV=8.75*ones(1,6);

prPos=[prPos;[tempAP;tempML;tempDV;ratIdx*ones(size(tempAP));prb*ones(size(tempAP))]'];

%hoegaarden181115 hpc
ratIdx=8;
prb=1;
tempAP=-5.05*ones(1,6);
tempML=4.2+0.15*(0:5)*cos(-14/180*pi);
tempDV=8.95+0.15*(0:5)*sin(-14/180*pi);

prPos=[prPos;[tempAP;tempML;tempDV;ratIdx*ones(size(tempAP));prb*ones(size(tempAP))]'];

%innis190601 pfc
ratIdx=9;
prb=2;
tempAP=3+(0:5)*0.2;
tempML=0.7*ones(1,6);
tempDV=3.8*ones(1,6);

prPos=[prPos;[tempAP;tempML;tempDV;ratIdx*ones(size(tempAP));prb*ones(size(tempAP))]'];

%innis190601 bla
ratIdx=9;
prb=3;
tempAP=-2.8*ones(1,6);
tempML=4.2+0.18*(0:5);
tempDV=8.3*ones(1,6);

prPos=[prPos;[tempAP;tempML;tempDV;ratIdx*ones(size(tempAP));prb*ones(size(tempAP))]'];

%innis190601 hpc
ratIdx=9;
prb=1;
tempAP=-4.95*ones(1,6);
tempML=4.4+0.15*(0:5)*cos(-14/180*pi);
tempDV=8.8+0.15*(0:5)*sin(-14/180*pi);

prPos=[prPos;[tempAP;tempML;tempDV;ratIdx*ones(size(tempAP));prb*ones(size(tempAP))]'];

%jever190814 pfc
ratIdx=10;
prb=2;
tempAP=3.1+(0:5)*0.2;
tempML=0.6*ones(1,6);
tempDV=3.7*ones(1,6);

prPos=[prPos;[tempAP;tempML;tempDV;ratIdx*ones(size(tempAP));prb*ones(size(tempAP))]'];

%jever190814 bla
ratIdx=10;
prb=3;
tempAP=-2.8*ones(1,6);
tempML=3.75+0.18*(0:5);
tempDV=8.8*ones(1,6);

prPos=[prPos;[tempAP;tempML;tempDV;ratIdx*ones(size(tempAP));prb*ones(size(tempAP))]'];

%jever190814 hpc
ratIdx=10;
prb=1;
tempAP=-5.8*ones(1,6);
tempML=4.4+0.18*(0:5)*cos(-14/180*pi);
tempDV=8.55+0.18*(0:5)*sin(-14/180*pi);

prPos=[prPos;[tempAP;tempML;tempDV;ratIdx*ones(size(tempAP));prb*ones(size(tempAP))]'];

%karmeliet190901 bla

ratIdx=11;
prb=3;
tempAP=-2.5*ones(1,6);
tempML=4.1+0.18*(0:5);
tempDV=8.15*ones(1,6);

prPos=[prPos;[tempAP;tempML;tempDV;ratIdx*ones(size(tempAP));prb*ones(size(tempAP))]'];

%karmeliet190901 hpc
ratIdx=11;
prb=1;
tempAP=-4.75*ones(1,6);
tempML=4.35+0.17*(0:5)*cos(-14/180*pi);
tempDV=8.75+0.17*(0:5)*sin(-14/180*pi);

prPos=[prPos;[tempAP;tempML;tempDV;ratIdx*ones(size(tempAP));prb*ones(size(tempAP))]'];

%karmeliet190901 pfc

ratIdx=11;
prb=2;
tempAP=3.1+(0:5)*0.2;
tempML=0.7*ones(1,6);
tempDV=4*ones(1,6);

prPos=[prPos;[tempAP;tempML;tempDV;ratIdx*ones(size(tempAP));prb*ones(size(tempAP))]'];


%leffe200124 hpc

ratIdx=12;
prb=1;
tempAP=-4.9*ones(1,6);
tempML=4.35+0.17*(0:5)*cos(-12/180*pi);
tempDV=8.73+0.17*(0:5)*sin(-12/180*pi);

prPos=[prPos;[tempAP;tempML;tempDV;ratIdx*ones(size(tempAP));prb*ones(size(tempAP))]'];

%leffe200124 pfc

ratIdx=12;
prb=2;
tempAP=2.5+(0:5)*0.2;
tempML=0.7*ones(1,6);
tempDV=4.05*ones(1,6);

prPos=[prPos;[tempAP;tempML;tempDV;ratIdx*ones(size(tempAP));prb*ones(size(tempAP))]'];

%leffe200124 bla

ratIdx=12;
prb=3;
tempAP=-3.2+0.08*(0:5);
tempML=4.7+0.11*(0:5);
tempDV=8.8*ones(1,6);
prPos=[prPos;[tempAP;tempML;tempDV;ratIdx*ones(size(tempAP));prb*ones(size(tempAP))]'];

%maredsous200224 hpc

ratIdx=13;
prb=1;
tempAP=-5.30*ones(1,6);
tempML=4.2+0.18*(0:5)*cos(-14/180*pi);
tempDV=8.75+0.18*(0:5)*sin(-14/180*pi);

prPos=[prPos;[tempAP;tempML;tempDV;ratIdx*ones(size(tempAP));prb*ones(size(tempAP))]'];

%maredsous200224 pfc

ratIdx=13;
prb=2;
tempAP=2.8+(0:5)*0.2;
tempML=0.6*ones(1,6);
tempDV=3.5*ones(1,6);

prPos=[prPos;[tempAP;tempML;tempDV;ratIdx*ones(size(tempAP));prb*ones(size(tempAP))]'];

%maredsous200224 bla

ratIdx=13;
prb=3;
tempAP=-3.1+0*(0:5);
tempML=4.4+0.2*(0:5);
tempDV=8.55*ones(1,6);

prPos=[prPos;[tempAP;tempML;tempDV;ratIdx*ones(size(tempAP));prb*ones(size(tempAP))]'];


%nostrum200304 hpc

ratIdx=14;
prb=1;
tempAP=-4.60*ones(1,6);
tempML=4.6+0.2*(0:5)*cos(-14/180*pi);
tempDV=8.4+0.2*(0:5)*sin(-14/180*pi);

prPos=[prPos;[tempAP;tempML;tempDV;ratIdx*ones(size(tempAP));prb*ones(size(tempAP))]'];

%nostrum200304 pfc

ratIdx=14;
prb=2;
tempAP=4.2-(0:5)*0.2;
tempML=0.7*ones(1,6);
tempDV=3.8*ones(1,6);

prPos=[prPos;[tempAP;tempML;tempDV;ratIdx*ones(size(tempAP));prb*ones(size(tempAP))]'];

%nostrum200304 bla

ratIdx=14;
prb=3;
tempAP=-2.35+0*(0:5);
tempML=4.35+0.22*(0:5);
tempDV=8.7*ones(1,6);
prPos=[prPos;[tempAP;tempML;tempDV;ratIdx*ones(size(tempAP));prb*ones(size(tempAP))]'];

%oberon200325 hpc

ratIdx=15;
prb=1;
tempAP=-5.2*ones(1,6);
tempML=4.6+0.2*(0:5)*cos(-14/180*pi);
tempDV=9+0.2*(0:5)*sin(-14/180*pi);

prPos=[prPos;[tempAP;tempML;tempDV;ratIdx*ones(size(tempAP));prb*ones(size(tempAP))]'];


%oberon200325 pfc

ratIdx=15;
prb=2;
tempAP=4.0-(0:5)*0.2;
tempML=0.8*ones(1,6);
tempDV=3.8*ones(1,6);
prPos=[prPos;[tempAP;tempML;tempDV;ratIdx*ones(size(tempAP));prb*ones(size(tempAP))]'];

%oberon200325 bla

ratIdx=15;
prb=3;
tempAP=-2.7+0*(0:5);
tempML=3.7+0.2*(0:5);
tempDV=8.6*ones(1,6);
prPos=[prPos;[tempAP;tempML;tempDV;ratIdx*ones(size(tempAP));prb*ones(size(tempAP))]'];
%%

col=flatColorMap(ceil(length(unique(prPos(:,4)))/3),[],[],[],false);
sym=repmat('ov^',1,size(col,1));
col=reshape(repmat(col,1,3)',3,[])';
gSize=1;

headCapital=@(x) [upper(x(1)),x(2:end)];


targetNames={'vental hippocampus','prefrontal cortex','amygdala',};

prOrder=[1,3,2];
height=[0];
wGap=3;
width=(110-wGap*4)/5;
hGap=4;
for m=1:3
    prIdx=prOrder(m);
    
    tempAP=prPos(prPos(:,5)==prIdx,1);
    tempML=prPos(prPos(:,5)==prIdx,2);
    tempDV=prPos(prPos(:,5)==prIdx,3);
    tempRat=prPos(prPos(:,5)==prIdx,4);
    
    [~,sliceID]=min(abs(atlasAP-tempAP),[],2);
    sliceList=unique(sliceID);
    switch prIdx
        case 1
            MLRange=[2.0,6.5];
            DVRange=[5.5,10.5];
        case 3
            MLRange=[2.0,6.5];
            DVRange=[5.5,10.5];
        case 2
            MLRange=[-0.5,4.0];
            DVRange=[0.3,5.3];
    end
    
    height(m+1)=width/diff(MLRange)*diff(DVRange);
    for n=1:length(sliceList)
        subplotInMM(x+(width+wGap)*mod(n-1,5),y+2*sum(height(1:m))+2*hGap*(m-1)+(n>5)*(height(m+1)+hGap),width,height(m+1))
        
        
        xRange=getML(MLRange,atlasAP(sliceList(n)));
        yRange=getDV(DVRange,atlasAP(sliceList(n)));
        xMin=xRange(1)-1;
        yMin=yRange(1)-1;
        image(ind2rgb(img(yRange(1):yRange(2),xRange(1):xRange(2),sliceList(n)),rgb2gray(map(:,:,sliceList(n)))))
        hold on
        
        grp=tempRat(sliceID==sliceList(n));
        gCol=col(unique(grp),:);
        gSym=sym(unique(grp));
        gscatter(getML(tempML(sliceID==sliceList(n)),tempAP(sliceID==sliceList(n)))-xMin,...
            getDV(tempDV(sliceID==sliceList(n)),tempAP(sliceID==sliceList(n)))-yMin,...
            grp,gCol,gSym,gSize,'doleg','off')
        
        xlim(xRange-xMin)
        ylim(yRange-yMin)
        ax=fixAxis;
        text2(0,-0.002, sprintf('%+0.2f mm',atlasAP(sliceList(n))),ax,'verticalALign','bottom','fontweight','normal','fontsize',5)
        
        axis off
        ax=fixAxis;
    end
    
    drawnow()
end
subplotInMM(x+(width+wGap)*(mod(n-1,5)+1),y+2*sum(height(1:m))+2*hGap*(m-1)+(n>5)*(height(m+1)+hGap),120-(width+wGap)*(mod(n-1,5)+1),height(m+1))

nCol=3;
legX=repmat(((1:nCol)-0.8)/nCol,1,ceil(length(sessionList)/nCol));
legY=repmat((ceil(length(sessionList)/nCol):-1:1)-0.5,nCol,1)/ceil(length(sessionList)/nCol);
legY=legY(:)';

legX=legX(1:length(sessionList));
legY=legY(1:length(sessionList));

gscatter(...
    legX,legY,...
    1:length(sessionList),...
    col,sym,gSize,'doleg','off')

for nn=1:length(sessionList)
    text(legX(nn)+0.02,legY(nn),headCapital(sessionList{nn}(1:end-6)),'fontsize',5)
end
xlim([0,1]);ylim([0,1])
axis off
ax=fixAxis;
text2(0,1,'Rat names',ax,'verticalAlign','bottom','horizontalALign','left','fontsize',5)
end


