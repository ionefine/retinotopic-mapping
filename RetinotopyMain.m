function [stim, eye, display]=RetinotopyMain(subjID, typeID, scanID)

% subjID='JMD1';
% typeID='bar'
% scanID='1';
close all hidden

place='Ione';
if strcmp(place, 'Ione')
    homedir =[ 'C:' filesep 'Users' filesep 'Ione Fine' filesep 'Documents' filesep ...
        'Work' filesep 'Science' filesep 'Projects' filesep 'Ione Fine' filesep '\retinotopic-mapping'];
elseif strcmp(place, 'You');
    homedir='what you want here';
end

cd(homedir);
rng('default'); rng(sum(100*clock)); stim.seed=rng; % saves the random seed in case you want to recreate the scan
stim.filename=[subjID, '_',typeID, '_', num2str(scanID)];

%% screen variables and trigger key
display.dist = 36.5;     %distance from screen (cm)
display.width = 27.5;  %width of screen (cm)
display.screenNum = 0;
display.skipChecks = 1;
display.bkColor = [0 0 0];
display.gamma = 0:255;% This is if the monitor is linear!!
% open and close the display in order to set some more generic display parameters
display = OpenWindow(display);
display.rect = Screen('Rect', display.windowPtr);
Screen('CloseAll');
display.screenAngle = pix2angle( display, display.resolution(1));
display.blankImg = 128*ones(display.rect(4), display.rect(3));

%% stimulus parameters, timing
% To alter other aspects of the stimuli go into MakeStimulus
stim.triggerKey = '5'; %key code for starting the scanner
stim.dur = 30;  %duration of scan (seconds)
if stim.dur<30 && strcmp(typeID, 'mf')
    error('Multifocal will produce weird segment orders if not adequate duration');
end

stim.type=typeID; %'mf' or 'bar'    

stim.logscaled=1; % determines whether eccentricity rings are log scaled to match cortical expansion
stim.expFac=1.2; %scales the virtual screen to to be larger than real screen. Can be 1 for real patients. 
% If simulating eye-movements needs to be proportional to the size of the
% eye-movements
if strcmp(stim.type, 'mf')
    stim.hz=1/3; % the rate at which mf stimulus updates
    % number of positions the bar/mf stimulus takes
elseif strcmp(stim.type, 'bar')
    stim.hz=8;% the rate for the updating of the bar, convenient to set it to match the background flicker rate
end
stim.npos = ceil(stim.dur/(1./stim.hz)); 
stim.switchtime=0:(1/stim.hz):stim.dur;

%% calculate some generic parameters relating the stimulus to the display screen
[stim.X,stim.Y] = meshgrid(linspace(-1,1,display.resolution(1)*stim.expFac));
stim.XDeg=stim.X.*display.screenAngle*stim.expFac/2;
stim.YDeg=stim.Y.*display.screenAngle*stim.expFac/2;
stim.RDeg=sqrt(stim.XDeg.^2+stim.YDeg.^2);
stim.T=atan2(stim.Y,stim.X)./pi; % angular rotation
stim.rect=CenterRect([1 1 display.resolution(1) display.resolution(2)],[ 1 1 size(stim.X)]);

%%  eye-movements
% modeling two types of eye-movements, slow drift and jitter
% for a real patient you want everything fixed (0)
eye.jitter.move=0; % this jitters the stimulus while keeping fix fixed, good for mimicking small rapid movements
% 0 - no jitter, 1- instability based on Bethlehem, Dumoulin, Dalmaijer, et al. Decreased Fixation Stability of the Preferred Retinal Location in Juvenile Macular Degeneration. Baker CI, ed. PLoS ONE. 2014;9(6)
eye.drift.move=[0 0; 0 0];% [4, 5; 6 3]; % this moves the fixation spot, good for big slow eye-movments.
% First row describes x, y drift distance (in degrees). Second row
% descripes the oscillation rate in degrees/seconds
eye.frameRate=display.frameRate; 
if eye.jitter.move+sum(eye.drift.move(1,:))>0 && stim.expFac==1; 
    error('stim.expFac must be greater than 1 to allow for eye-movements'); 
end
[eye, stim]=SimulateEyeMovements(eye, stim);

%% fixation parameters
fix.offset=[0 0]; %location of the fixation spot in the display, x y degrees for - +ve x y values the fix spot moves to the left or to the upper part of the screen
fix.scFac=8; % size of the fix spot, 1 is the default

%%  background
stim.background.nsectors=18; % controls how fine the background texture is, should be finer scale than the mf stimulus
stim.background.nrings=40;
stim.background.n=3; % how many jitter positions for the texture
stim.background.hz=8; % rate at which the background flickers
stim.background.switchtime=0:(1/stim.background.hz):stim.dur;
if stim.background.hz<stim.hz
    error('Background must update at the same or faster rate than stimulus')
end
stim=MakeBackground(stim);

%% scotoma 
stim.scotoma.rad = 1;
stim.scotoma.center = [0 0]; % x, y in degrees
% for -  positive x y values the fix spot moves to the left or to the upper part of the screen

if stim.scotoma.rad>0
    stim=MakeScotoma(stim);
else
    stim.scotoma.img=ones(size(stim.X));
end
if strcmp(stim.type, 'bar')
    stim.bar.speed = [2,3];  %speed range degrees/sec
    stim.bar.width=6;
    stim=MakeBar(stim);
elseif strcmp(stim.type, 'mf')
    stim.mf.nrings=8; % total number of rings
    stim.mf.inner=.25; % radius of the inner ring
    stim.mf.sectors=12; % number of angular sectors
    stim.mf.nrot=1; % number of times and direction it rotates back and forth the angular width of a segment in the scan
    stim.mf.nEC=3; % number of times it expands and contracts
    stim.mf.noneighbors=0; %neigborFlag flag to determine if neighbors cannot be on at the same time
    stim=MakeMF(stim);
end

%try
    display = OpenWindow(display); HideCursor;
    
    %% initialize Simon task
    task.apertures=[1,3,3.3,.1]; % inner edge of the segments, outer edge of segments, outer edge, inner hole
    task.action = 'init';   task.fixoffset=round(angle2pix(display, fix.offset)/2);
    task.rotAng=45; %  can be 0 (vert), -45 , 45 (diag, pointing towards fix) or 90 (horiz)
    task.ISI = 1/3;  %seconds
    task.dur = .25;   %seconds
    task.pauseDur = .25; %seconds
    task.errDur = 1;  %seconds
    task.errFreq = 4;  %Hz
    task.keys = {'1' '2','3' '4'};
    task.goodKeys = [1:4];

    task= doSimonScotoma2color(display,task);

    
    task.rectInit=task.rect; task.cInit=task.c; % this is the 'center' eye-position
    stimflipNum=1; % keeps track of stimulus flips
    backflipNum=1; on=1;% keeps track of background flips and sign
    
    t = stim.background.switchtime; % keeps track of when to switch the stimulus
     startTime = GetSecs;  %read the clock
    for fN=1:length(t)    
        curTime = GetSecs-startTime;
        if curTime>stim.switchtime(1) % this decides whether to update the stimulus
            stim.switchtime=stim.switchtime(2:end);
            stimflipNum=mod(stimflipNum,size(stim.img, 3)-1)+1;
        end
        if on
            backimg=stim.background.img(:,:,mod(backflipNum, size(stim.background.img, 3))+1).*255; on=0;
        else
            backimg=(~stim.background.img(:,:,mod(backflipNum, size(stim.background.img, 3))+1)).*255; on=1;
        end
        % decide what part of the stimulus/scotoma to sample from, given the eye-position of the subject,
        % includes jitter but not drift, this is for the stimulus and background
        rect(1,:)=stim.rect-angle2pix(display, [eye.jitter.pos(fN,1) eye.jitter.pos(fN,2) eye.jitter.pos(fN,1) eye.jitter.pos(fN,2)]);
        % includes drift but not jitter, this is for the scotoma
        rect(2,:)=stim.rect-angle2pix(display, [eye.drift.pos(fN,1) eye.drift.pos(fN,2) eye.drift.pos(fN,1) eye.drift.pos(fN,2)]);
        
        img=squeeze(double(stim.img(rect(1, 2):rect(1, 4), rect(1, 1):rect(1, 3), stimflipNum))) .* ...
            squeeze(backimg(rect(1, 2):rect(1, 4), rect(1, 1):rect(1, 3))) .* ...
            squeeze(stim.scotoma.img(rect(2, 2):rect(2, 4), rect(2, 1):rect(2, 3)));
        
        img = img(1:end-1,1:end-1);
        % do the simon task 
        task.fixoffset=angle2pix(display, fix.offset)/2+angle2pix(display, eye.drift.pos(fN,:));       
        
        
        textureIndex=Screen('MakeTexture', display.windowPtr, img);
        Screen('DrawTexture', display.windowPtr, textureIndex, display.rect);
        task=doSimonScotoma2color(display,task,curTime);
        Screen('Flip', display.windowPtr);
      %  s.fixoffset=fix.offsetP+eye.drift.posP(fN,:);
    end
% catch ME
%     Screen('CloseAll');
%     ListenChar(0);
%     ShowCursor;
% end 
  Screen('CloseAll');
    ListenChar(0);
    ShowCursor;
    
end

function stim=MakeBackground(stim)
if stim.logscaled==1
    stim.background.ringpos=exp(linspace(.0001, log(max(stim.XDeg(:))*1.2),stim.background.nrings+1));
else
    stim.background.ringpos=linspace(.0001, max(stim.XDeg(:))*1.2, stim.background.nrings+1);
end
stim.background.img=zeros([size(stim.XDeg), stim.background.n.^2]);
disp('Generating background ...');
bNum=1;
for rh=1:stim.background.n % for each radial background jitter
    rings=zeros(size(stim.XDeg));
    for i=1:length(stim.background.ringpos)-1; % fill in each ring
        offset=rh*(diff(stim.background.ringpos(i:i+1))./stim.background.n);
        rings(stim.RDeg>stim.background.ringpos(i)+offset & stim.RDeg<=stim.background.ringpos(i+1)+offset)=i;
    end
    for th=1:stim.background.n % for each angular jitter position of the background
        segs=zeros(size(stim.XDeg));
        segs(mod(stim.T*stim.background.nsectors+(th*1/stim.background.n), 1)>0.5)=1;
        stim.background.img(:, :, bNum)=mod(rings+segs, 2)==0;
        bNum=bNum+1;
    end
end
end

function [eye, stim]=SimulateEyeMovements(eye, stim)
eye.nFrames=ceil((stim.dur*eye.frameRate)+50); % guess at how many frames will be required

if eye.jitter.move==0;
    eye.jitter.pos=zeros(eye.nFrames, 2);
elseif eye.jitter.move==1 % estimated instability based on
    % Bethlehem, Dumoulin, Dalmaijer, et al. Decreased Fixation Stability of the Preferred Retinal Location in Juvenile Macular Degeneration. Baker CI, ed. PLoS ONE. 2014;9(6)
    eye.jitter.pos(:, 1)=linspace(0, .5,eye.nFrames)+.25.*(randn(1, eye.nFrames)); % in degrees
    eye.jitter.pos(:, 2)=linspace(0, -.5,eye.nFrames)+.25.*(randn(1, eye.nFrames));
end
if sum(eye.drift.move(1, :))==0
    eye.drift.pos=zeros(eye.nFrames, 2);
else 
    eye.drift.pos(:, 1)=eye.drift.move(1,1)*sin(linspace(0, eye.drift.move(2, 1)*pi,eye.nFrames)+.25); % in degrees, simulated drift % IFCHANGE
    eye.drift.pos(:, 2)=eye.drift.move(1,2)*sin(linspace(0, eye.drift.move(2,2)*pi,eye.nFrames)+.25); % in degrees
end
end

function stim=MakeScotoma(stim)
disp('Generating scotoma ...');
scot=zeros(size(stim.XDeg));
scot(sqrt((stim.XDeg-stim.scotoma.center(2)).^2+(stim.YDeg-stim.scotoma.center(1)).^2)<stim.scotoma.rad)=0;
scot(sqrt((stim.XDeg-stim.scotoma.center(2)).^2+(stim.YDeg-stim.scotoma.center(1)).^2)>=stim.scotoma.rad)=1;
stim.scotoma.img=scot;
end

function stim=MakeBar(stim)
%% create the bar stimulus
% will create a mask for the bar stimulus. whatever the Hz rate of updating
% the stimulus is ~8Hz, will have the necessary number of bars.
h = waitbar(0,'Generating bar stimulus');
stim.bar.speed = [2,3];  %speed range degrees/sec
stim.bar.width=6;
direction = rand(1)*360;
speed = rand(1)*(stim.bar.speed(2)-stim.bar.speed(2))+stim.bar.speed(1);
center = -(max(stim.XDeg(:))+stim.bar.width/2)*[cos(direction*pi/180),sin(direction*pi/180)];
dx = cos(direction*pi/180)*speed*(1/stim.background.hz);
dy = sin(direction*pi/180)*speed*(1/stim.background.hz);

stim.center = zeros(stim.npos,2);
stim.ang = zeros(stim.npos,1);
stim.img = zeros([size(stim.XDeg),stim.npos],'uint8');
for flipNum = 1:stim.npos
    center = center + [dx,dy];
    if norm(center)>max(stim.XDeg(:))+stim.bar.width/2;
        direction = rand(1)*360;
        center = -(max(stim.XDeg(:))+stim.bar.width/2)*[cos(direction*pi/180),sin(direction*pi/180)];
        dx = cos(direction*pi/180)*speed*(1/stim.background.hz);
        dy = sin(direction*pi/180)*speed*(1/stim.background.hz);
        speed = rand(1)*(stim.bar.speed(2)-stim.bar.speed(2))+stim.bar.speed(1);
    end
    stim.center(flipNum,:) = center;
    stim.ang(flipNum) = direction+90;
    %    mask = false([size(XDeg),stim.npos]);
    %       a = cos(stim.ang(flipNum)*pi/180)*(XDeg-stim.center(flipNum,1)) + sin(stim.ang(flipNum)*pi/180)*(YDeg-stim.center(flipNum,2));
    b = -sin(stim.ang(flipNum)*pi/180)*(stim.XDeg-stim.center(flipNum,1)) +cos(stim.ang(flipNum)*pi/180)*(stim.YDeg-stim.center(flipNum,2));
    stim.img(:,:,flipNum) = uint8((abs(b)<stim.bar.width/2));
    waitbar(flipNum/stim.npos,h);
end
delete(h);
end

function stim=MakeMF(stim)

h = waitbar(0,'Generating  multifocal stimulus');
stim.mf.outer=max(stim.XDeg(:)); %radius outer ring in degrees
stim.stimSeq = generateMultifocalSequence(stim.mf.sectors,stim.mf.nrings,stim.npos, stim.mf.noneighbors);
if stim.logscaled==1
    stim.mf.ringpos=exp(linspace(log(stim.mf.inner), log(stim.mf.outer), stim.mf.nrings+1));% radii position of each ring
else
    stim.mf.ringpos=linspace(stim.mf.inner, stim.mf.outer, stim.mf.nrings+1);
end
% describe how the stimulus rotates across the scan, this gradual rotation
% allows pRFs to be estimated more accurately visual space is sampled more
% finely
a=linspace(0, (360/stim.mf.sectors), ceil(stim.npos/(stim.mf.nrot*2))+1);
a=[a fliplr(a(2:end-1))]; % add the reversal
stim.mf.rotpos=repmat(a, 1, ceil(stim.npos/length(a)));
stim.mf.rotpos=stim.mf.rotpos(1:stim.npos);

% describe how the stimulus expands/contracts across the scan, ditto
maxscFac=(max(stim.XDeg(:))-diff(stim.mf.ringpos(1:2)))./max(stim.XDeg(:));
a=linspace(1, maxscFac, ceil(stim.npos/(stim.mf.nEC*2))+1);
a=[a fliplr(a(2:end-1))]; % add the reversal
stim.mf.ECscale=repmat(a, 1, ceil(stim.npos/length(a)));
stim.mf.ECscale=stim.mf.ECscale(1:stim.npos);

rad=zeros(size(stim.XDeg));
for i=1:length(stim.mf.ringpos)-1 % fill in the rings
    rad(stim.RDeg>stim.mf.ringpos(i) & stim.RDeg<=stim.mf.ringpos(i+1))=i.*stim.mf.sectors;
end
% calculate the sectors
theta=ceil(scaleif(stim.T,0,stim.mf.sectors)).*(rad>0);% window by the rings
finalimg=((theta+rad)-stim.mf.sectors).*(rad>0);

for flipNum = 1:stim.npos % for each version of the MF sequence
    img=zeros(size(finalimg));
    onlist=find(stim.stimSeq(:, flipNum));
    for o=1:length(onlist)
        img(finalimg==onlist(o))=1;
    end
    % add black lines to avoid aliasing
    %  img(mod(T, 1)<0.05)=0;
    % rotate the image
    img = imrotate(img,stim.mf.rotpos(flipNum),'nearest', 'crop');
    % scale the image
    img = imresize(img, stim.mf.ECscale(flipNum), 'nearest');
    stim.img(:,:,flipNum) = insertImg(zeros(size(stim.X)), img);
    waitbar(flipNum/stim.npos,h);
end
delete(h);

end

