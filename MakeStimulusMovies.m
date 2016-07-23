function stim=MakeStimulus(display, stim)



% needs to be able to contain the max amount of eyewobble.
[X,Y] = meshgrid(linspace(-1,1,display.resolution(1)*stim.expFac));
T=atan2(Y,X)./pi;
XDeg=X.*display.screenAngle*stim.expFac/2;
YDeg=Y.*display.screenAngle*stim.expFac/2;
RDeg=sqrt(XDeg.^2+YDeg.^2);

stim.rect=CenterRect([0 0 display.resolution(1) display.resolution(2)],[ 0 0 size(X)]);

%% create the background
stim.background.nsectors=18; % controls how fine the texture is, should be different from the mf stimulus
stim.background.nrings=40;
stim.background.n=3; % how many jitter positions for the texture
bNum=1;
stim.background.ringpos=exp(linspace(.0001, log(max(XDeg(:))), stim.background.nrings+1));
stim.background.img=zeros([size(XDeg), stim.background.n.^2]);
disp('Generating background ...');
writerObj = VideoWriter([filename, '_Background.avi'], 'Grayscale AVI');
open(writerObj);
for rh=1:stim.background.n % for each radial background jitter
    rings=zeros(size(XDeg));
    for i=1:length(stim.background.ringpos)-1; % fill in each ring
        offset=rh*[diff(stim.background.ringpos(i:i+1))./stim.background.n];
        rings(RDeg>stim.background.ringpos(i)+offset & RDeg<=stim.background.ringpos(i+1)+offset)=i;
    end
    for th=1:stim.background.n % for each angular jitter position of the background
        segs=zeros(size(XDeg));
        segs(mod(T*stim.background.nsectors+[th*1/stim.background.n], 1)>0.5)=1;
        frame.cdata=mod(rings+segs, 2)==0;
        frame.colormap=[];
        writeVideo(writerObj,frame);
        bNum=bNum+1;
    end
end
break

%% make the scotoma

if stim.scotoma.on>0
    stim.scotoma.img=uint8(NaN(size(XDeg, 1), size(XDeg, 2)));
    disp('Generating scotoma ...');
    scot=zeros(size(XDeg));
    scot(sqrt((XDeg-stim.scotoma.center(2)).^2+(YDeg-stim.scotoma.center(1)).^2)<stim.scotoma.rad)=0;
    scot(sqrt((XDeg-stim.scotoma.center(2)).^2+(YDeg-stim.scotoma.center(1)).^2)>=stim.scotoma.rad)=1;
    stim.scotoma.img=scot;
else
    stim.scotoma.img=uint8(ones(size(X)));
end

%% create the bar stimulus
% will create a mask for the bar stimulus. whatever the Hz rate of updating
% the stimulus is ~8Hz, will have the necessary number of bars.

if strcmp(stim.type, 'bar')
    stim.img=uint8(NaN(size(XDeg, 1), size(XDeg, 2),stim.npos));
    disp('Generating bar stimulus ...');
    stim.bar.speed = [2,3];  %speed range degrees/sec
    stim.bar.width=6;
    direction = rand(1)*360;
    speed = rand(1)*(stim.bar.speed(2)-stim.bar.speed(2))+stim.bar.speed(1);
    center = -(max(XDeg(:))+stim.bar.width/2)*[cos(direction*pi/180),sin(direction*pi/180)];
    dx = cos(direction*pi/180)*speed*[1/stim.hz];
    dy = sin(direction*pi/180)*speed*[1/stim.hz];
    
    stim.center = zeros(stim.npos,2);
    stim.ang = zeros(stim.npos,1);
    for flipNum = 1:stim.npos
        center = center + [dx,dy];
        if norm(center)>max(XDeg(:))+stim.bar.width/2;
            direction = rand(1)*360;
            center = -(max(XDeg(:))+stim.bar.width/2)*[cos(direction*pi/180),sin(direction*pi/180)];
            dx = cos(direction*pi/180)*speed*[1/stim.hz];
            dy = sin(direction*pi/180)*speed*[1/stim.hz];
            speed = rand(1)*(stim.bar.speed(2)-stim.bar.speed(2))+stim.bar.speed(1);
        end
        stim.center(flipNum,:) = center;
        stim.ang(flipNum) = direction+90;
    end
    %    mask = false([size(XDeg),stim.npos]);
    
    for flipNum = 1:stim.npos
        %       a = cos(stim.ang(flipNum)*pi/180)*(XDeg-stim.center(flipNum,1)) + sin(stim.ang(flipNum)*pi/180)*(YDeg-stim.center(flipNum,2));
        b = -sin(stim.ang(flipNum)*pi/180)*(XDeg-stim.center(flipNum,1)) +cos(stim.ang(flipNum)*pi/180)*(YDeg-stim.center(flipNum,2));
        stim.img(:,:,flipNum) = uint8((abs(b)<stim.bar.width/2));
    end
elseif strcmp(stim.type, 'mf')
    stim.img=uint8(NaN(size(XDeg, 1), size(XDeg, 2),stim.npos));
    disp('Generating multifocal stimulus ...');
    
    %% create the MF stimulus
    stim.mf.nrings=8; % total number of rings
    stim.mf.inner=.5; % radius of the inner ring
    stim.mf.outer=max(XDeg(:)); %radius outer ring in degrees
    stim.mf.sectors=12; % number of angular sectors
    stim.mf.nrot=1; % number of times and direction it rotates back and forth the angular width of a segment in the scan
    stim.mf.nEC=3; % number of times it expands and contracts
    stim.mf.noneighbors=0; %
    stim.stimSeq = generateMultifocalSequence(stim.mf.sectors,stim.mf.nrings,stim.npos,  stim.mf.noneighbors);
    
    stim.mf.ringpos=exp(linspace(log(stim.mf.inner), log(stim.mf.outer), stim.mf.nrings+1));% radii position of each ring
    
    % describe how the stimulus rotates across the scan
    a=linspace(0, (360/stim.mf.sectors), ceil(stim.npos/(stim.mf.nrot*2))+1);
    a=[a fliplr(a(2:end-1))]; % add the reversal
    stim.mf.rotpos=repmat(a, 1, ceil(stim.npos/length(a)));
    stim.mf.rotpos=stim.mf.rotpos(1:stim.npos);
    
    % describe how the stimulus expands/contracts across the scan
    maxscFac=(max(XDeg(:))-diff(stim.mf.ringpos(1:2)))./max(XDeg(:));
    a=linspace(1, maxscFac, ceil(stim.npos/(stim.mf.nEC*2))+1);
    a=[a fliplr(a(2:end-1))]; % add the reversal
    stim.mf.ECscale=repmat(a, 1, ceil(stim.npos/length(a)));
    stim.mf.ECscale=stim.mf.ECscale(1:stim.npos);
    
    rad=zeros(size(XDeg));
    for i=1:length(stim.mf.ringpos)-1 % fill in the rings
        rad(RDeg>stim.mf.ringpos(i) & RDeg<=stim.mf.ringpos(i+1))=i.*stim.mf.sectors;
    end
    % calculate the sectors
    theta=ceil(scaleif(T,0,stim.mf.sectors)).*(rad>0);% window by the rings
    finalimg=((theta+rad)-stim.mf.sectors).*(rad>0);
    
    for flipNum = 1:stim.npos % for each version of the MF sequence
        img=zeros(size(finalimg));
        onlist=find(stim.stimSeq(:, flipNum));
        for o=1:length(onlist)
            img(finalimg==onlist(o))=1;
        end
        % add black lines to avoid aliasing
        img(mod(T, 1)<0.05)=0;
        
        % rotate the image
        img = imrotate(img,stim.mf.rotpos(flipNum),'nearest', 'crop');
        % scale the image
        img = imresize(img, stim.mf.ECscale(flipNum), 'nearest');
        stim.img(:,:,flipNum) = uint8(insertImg(zeros(size(X)), img));
    end
end
end
