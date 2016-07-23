% PlayYorkRetinotopy
% This is f

clear all
clear MEX

subjID='JMD1'; % subject ID
stimID='mf'; % or can be 'bar'
scanID='1'; 

% avi file to read inq
filename=[subjID, '_', stimID, '_', scanID]; 

% create filename to save to
c=clock;

filename_out=[filename, '_', ... 
    num2str(c(5)), '_', num2str(c(4)),'_', num2str(c(3)), '_',num2str(c(2)), '_',num2str(c(1))];
triggerKey = '5'; %key code for starting the scanner

place='Ione';
if strcmp(place, 'Ione')
    homedir =[ 'C:' filesep 'Users' filesep 'Ione Fine' filesep 'Documents' filesep ...
        'Work' filesep 'Science' filesep 'Projects' filesep 'Ione Fine' filesep 'York_retinotopy' filesep ...
        'Programs_stimulus' filesep 'FineCode'];
elseif strcmp(place, 'Geoff');
    homedir='C:\Users\Geoff Boynton\Documents\MATLAB\YorkRetinotopy';
end
addpath(genpath(homedir));
cd(homedir); % cd('FineCode');

load(filename);

%% open and close the display and collect some parameters
display = OpenWindow(display);
display.rect = Screen('Rect', display.windowPtr);
display.screenAngle = pix2angle( display, display.resolution(1));
display.blankImg = 128*ones(display.rect(4), display.rect(3));
Screen('CloseAll');

%% start running it
try
    
    xyloObj = VideoReader([filename,'.avi']);
    display = OpenWindow(display);
   % commandwindow; % puts command window in front, good for debugging but puts the menu bar up
   
    %initialize Simon
    % update simon position & simon
    s.action = 'init';
    s.fixoffset=fix.offsetP/2;
    s.rotAng=45; % IFCHANGE, can be 0 (vert), -45 , 45 (diag, pointing towards fix) or 90 (horiz)
    s= doSimonScotoma2color(display,s);
    s.ISI = 1/3;  %seconds
    s.dur = .25;   %seconds
    s.pauseDur = .25; %seconds
    s.errDur = 1;  %seconds
    s.errFreq = 4;  %Hz
    s.keys = {'1' '2','3' '4'};

    s.rectInit=s.rect; s.cInit=s.c; % this i5s the 'center' eye-position
    Screen('Flip',display.windowPtr);
    
    wait4T(triggerKey);  %wait for 't' from scanner.
    startTime = GetSecs;  %read the clock

    fN=0; % keeps track of frames
    escPressed=0;
    count = 0;
    
    while GetSecs-startTime< dur && ~escPressed  %loop until 'esc' pressed or time runs out
        
        %calculate the current time elapsed, in seconds
        curTime = GetSecs-startTime;
        count = count+1;
        fN = find(t>=curTime,1,'first');
        img = read(xyloObj,fN);
        
        % check for escape key, to quit early
        %         [ keyIsDown, timeSecs, keyCode ] = KbCheck;
        %         if keyIsDown
        %             keyPressed= KbName(keyCode);
%             if strcmp(keyPressed, 'esc')
%                 escPressed=1;
%             end
%         end
%          
        % put up the stimulus
        textureIndex=Screen('MakeTexture', display.windowPtr, img);
        Screen('DrawTexture', display.windowPtr, textureIndex, display.rect);
       
        % do the simon task 
        s.fixoffset=fix.offsetP/2+eye.drift.posP(fN,:);       
        s=doSimonScotoma2color(display,s,curTime);
        Screen('Flip',display.windowPtr,t(count)+startTime);
        
    end
    %clean up the Simon process
    s.action = 'done';
    s= doSimonScotoma2color(display,s,curTime);
    save(filename_out);
    Screen('Flip',display.windowPtr);
    ShowCursor;  
    Screen('CloseAll');
    ListenChar(0);
    plotSimon2(s);
    
catch ME
    save(filename_out);
    Screen('CloseAll');
    ListenChar(0);
    rethrow(ME);
    ShowCursor;
end










