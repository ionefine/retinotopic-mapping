%YorkRetinotopy


clear all
clear MEX

place='Ione';
if strcmp(place, 'Ione')
homedir =[ 'C:' filesep 'Users' filesep 'Ione Fine' filesep 'Documents' filesep ...
    'Work' filesep 'Science' filesep 'Projects' filesep 'Ione Fine' filesep 'York_retinotopy' filesep ...
    'Programs_stimulus' filesep 'FineCode'];
elseif strcmp(place, 'Geoff');
    homedir='C:\Users\Geoff Boynton\Documents\MATLAB\YorkRetinotopy';
end
addpath(genpath(homedir));
cd(homedir); 

rng('default')
rng(sum(100*clock));
stim.seed=rng; % saves the random seed in case you want to recreate the scan

subjID='JMD1';
typeID='bar'
scanID='1';

c=clock;
stim.filename=[subjID, '_',typeID, '_', scanID];

%% screen variables and trigger key
display.dist = 36.5;     %distance from screen (cm)
display.width = 27.5;  %width of screen (cm)
display.screenNum = 0;
display.skipChecks = 1;
display.bkColor = [0 0 0];
display.gamma = 0:255;% York monitor is linear
triggerKey = '5'; %key code for starting the scanner

%% open and close the display and collect some parameters
display = OpenWindow(display);
display.rect = Screen('Rect', display.windowPtr);
display.screenAngle = pix2angle( display, display.resolution(1));
display.blankImg = 128*ones(display.rect(4), display.rect(3));
Screen('CloseAll');

%% stimulus parameters
% To alter other aspects of the stimuli go into MakeStimulus

stim.dur = 330;  %duration of scan (seconds)
stim.type='bar'; %'mf' or 'bar'
if strcmp(stim.type, 'bar')
    stim.hz=8; % the rate at which bar stimulus updates
else
    stim.hz=1/3; % the rate at which mf stimulus updates
end
stim.nFrames=stim.dur*ceil(display.frameRate)+50; % guess at how many frames will be required
stim.switchtime=0:1/stim.hz:stim.dur;
stim.npos = stim.dur/(1./stim.hz); % number of positions the bar/mf stimulus takes

%% background parameters
stim.background.hz=8; % the rate at which the background updates;
stim.background.switchtime=0:1/[stim.background.hz*2]:stim.dur; %on and off presented in the hz

%%  eye-movement parameters
% modeling two types of eye-movements, slow drift and jitter
% % for a patient you want everything fixed
% eye.jitter.move=0;
% eye.drift.move=0; s

eye.jitter.move=0; % this jitters the stimulus while keeping fix fixed, good for mimicking small rapid movements
eye.drift.move=0; % this moves the fixation spot, good for big slow movements
if [eye.jitter.move+eye.drift.move]==0
    stim.expFac=1.2; % expansion factor creating a supersized screen on which the stimulus is made
else
     stim.expFac=1;
end

if eye.jitter.move==0;
    eye.jitter.pos=zeros(stim.nFrames, 2);
elseif eye.jitter.move==1 % estimated instability based on
    % Bethlehem, Dumoulin, Dalmaijer, et al.
    % Decreased Fixation Stability of the Preferred Retinal Location
    % in Juvenile Macular Degeneration. Baker CI, ed. PLoS ONE. 2014;9(6)
    eye.jitter.pos(:, 1)=linspace(0, .5,stim.nFrames)+.25.*(randn(1, stim.nFrames)); % in degrees
    eye.jitter.pos(:, 2)=linspace(0, -.5,stim.nFrames)+.25.*(randn(1, stim.nFrames));
end
if eye.drift.move==0
    eye.drift.pos=zeros(stim.nFrames, 2);
elseif eye.drift.move==1
    eye.drift.pos(:, 1)=4*sin(linspace(0, 6*pi,stim.nFrames)+.25); % in degrees
    eye.drift.pos(:, 2)=5*sin(linspace(0, 6*pi,stim.nFrames)+.25); % in degrees
end

eye.jitter.posP=angle2pix(display, eye.jitter.pos);
eye.drift.posP=angle2pix(display, eye.drift.pos);

%% scotoma parameters
stim.scotoma.on =0; % 0: no scotoma, 1: scotoma,  2 is York JMD1
if stim.scotoma.on>0
    % scotoma is registered to the center of the stimulus
    if stim.scotoma.on==1 % default scotoma
        stim.scotoma.center = [0 0]; % x, y in degrees 
        % for -  x y values the fix spot moves to the left or to the upper part of the screen
        stim.scotoma.rad = [1];
    elseif stim.scotoma.on==2 % personalized scotoma
        stim.scotoma.center=[0 0];
        stim.scotoma.rad=[7.5];
    end
end
stim=MakeStimulus(display, stim); % details of the scotoma, background, mf and bar can be changed inside this
