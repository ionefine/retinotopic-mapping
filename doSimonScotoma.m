
function  s=doSimonScotoma(display,s,curTime)
%s=doSimon(display,s,curTime)
escPressed=0;

switch s.action
    case 'init'  %set up simon variables
        % 'simon' parameters
        if ~isfield(s, 'scFac')
            s.scFac=1;
        end
        s.size = [6,20,24,3]*s.scFac;  %pixels
        s.gap = 2; %pixels
        s.c = ceil(display.resolution/2)+s.fixoffset;
        for i=1:length(s.size)
            s.rect{i} = [-s.size(i)/2+s.c(1),-s.size(i)/2+s.c(2),+s.size(i)/2+s.c(1),+s.size(i)/2+s.c(2)];
        end
        
        s.ISI = 1;  %seconds
        s.dur = .75;   %seconds
        s.pauseDur = .5; %seconds
        s.errDur = 1;  %seconds
        s.errFreq = 4;  %Hz
        
        s.color{1} = [0,0,255];
        s.color{2} = [255,255,0];
        s.color{3} = [0,255,0];
        s.color{4} = [255,0,0];
        
        s.goodKeys = [1,2,3,4];
        s.keys = {'b','y','g','r'};
        
        %s.keys = {'w','s','a','q'};  %for 1,2,3,4
        
        s.seq = [];
        s.maxLen = [];
        s.i = 0;
        
        s.state = 'off';
        
        s.seq = s.goodKeys(ceil(rand(1)*length(s.goodKeys)));
        s.pauseTill = 1;
        s.action = 'play';
        
        s.event.type = 'switch to play';
        s.event.time = 0;
        s.event.num = 0;
        s.numEvents = 1;
        
        showSimon(display,s);
        
    case 'play'  %sequence playback mode
        if s.pauseTill < curTime
            t = mod(curTime,s.ISI);
            
            if strcmp(s.state,'off') && t<=s.dur
                
                s.i = s.i+1;
                
                if s.i<=length(s.seq);
                    s.show = s.seq(s.i);
                    s.state = 'on';
                    
                    s.numEvents = s.numEvents+1;
                    s.event(s.numEvents).type= 'play: on';
                    s.event(s.numEvents).time = curTime;
                    s.event(s.numEvents).num = s.show;
                    
                else
                    s.numEvents = s.numEvents+1;
                    s.event(s.numEvents).type= 'switch to recall';
                    s.event(s.numEvents).time = curTime;
                    s.event(s.numEvents).num = NaN;
                    
                    s.action = 'recall';
                    s.state = 'off';
                    s.i = 0;
                end
                
            elseif strcmp(s.state,'on') && t> s.dur
                s.numEvents = s.numEvents+1;
                s.event(s.numEvents).type= 'play: off';
                s.event(s.numEvents).time = curTime;
                s.event(s.numEvents).num = s.show;
                s.state = 'off';
                
            end
        end
        showSimon(display,s);
        
    case 'recall'
        if s.pauseTill<curTime
            
            %see if a key is pressed
            [ keyIsDown, timeSecs, keyCode ] = KbCheck;
            
            if keyIsDown && strcmp(s.state,'off') %a key is down: record the key and time pressed
                
                keyPressed= KbName(keyCode);
                
                keyNum = find(strcmp(keyPressed(end),s.keys));
                
                if ~isempty(keyNum)
                    s.show = keyNum;
                    s.state = 'on';
                    s.numEvents = s.numEvents+1;
                    s.event(s.numEvents).type= 'recall: on';
                    s.event(s.numEvents).time = curTime;
                    s.event(s.numEvents).num = s.show;
                end
                
            elseif ~keyIsDown && strcmp(s.state,'on') %key just lifted
                s.state = 'off';
                s.numEvents = s.numEvents+1;
                s.event(s.numEvents).type= 'recall: off';
                s.event(s.numEvents).time = curTime;
                s.event(s.numEvents).num = NaN;
                
                s.i = s.i+1;
                if s.show ~= s.seq(s.i)
                    %incorrect
                    s.maxLen = [s.maxLen,length(s.seq)-1];
                    s.pauseTill = ceil(curTime+ s.errDur);
                    s.action = 'error';
                    
                    s.numEvents = s.numEvents+1;
                    s.event(s.numEvents).type= 'error: start';
                    s.event(s.numEvents).time = curTime;
                    s.event(s.numEvents).num = s.seq(end);
                    
                elseif s.i == length(s.seq)  %got the last one right
                    s.seq = [s.seq,s.goodKeys(ceil(rand(1)*length(s.goodKeys)))];
                    s.state = 'off';
                    s.action = 'play';
                    s.i = 0;
                    s.pauseTill = ceil(curTime+ s.pauseDur) ;
                    
                    s.numEvents = s.numEvents+1;
                    s.event(s.numEvents).type= 'switch to play';
                    s.event(s.numEvents).time = curTime;
                    s.event(s.numEvents).num = NaN;
                end
            end
        end
        showSimon(display,s);
    case 'error'
        if curTime < s.pauseTill
            ab = mod(curTime*s.errFreq,1);
            s.show = s.seq(end);
            if ab<.5
                s.state = 'on';
            else
                s.state = 'off';
            end
            showSimon(display,s);
        else
            s.i = 0;
            s.seq = s.goodKeys(ceil(rand(1)*length(s.goodKeys)));
            s.action = 'play';
            s.state = 'off';
            s.pauseTill = ceil(curTime+ s.pauseDur);
            s.numEvents = s.numEvents+1;
            s.event(s.numEvents).type= 'error: stop';
            s.event(s.numEvents).time = curTime;
            s.event(s.numEvents).num = NaN;
        end
    case 'done'
        s.maxLen = [s.maxLen,length(s.seq)-1];
        s.numEvents = s.numEvents+1;
        s.event(s.numEvents).type= 'done';
        s.event(s.numEvents).time = curTime;
        s.event(s.numEvents).num = NaN;
        s.seq = [];
end
end

function showSimon(display,s)
s.c = ceil(display.resolution/2)+s.fixoffset;
for i=1:length(s.size)
    s.rect{i} = [-s.size(i)/2+s.c(1),-s.size(i)/2+s.c(2),+s.size(i)/2+s.c(1),+s.size(i)/2+s.c(2)];
end
switch s.state
    case 'on';
        Screen('FillOval',display.windowPtr,[100,100,100],s.rect{3});
        Screen('FillArc',display.windowPtr,s.color{s.show},s.rect{2},90*(s.show-1),90);
        Screen('FillRect',display.windowPtr,[128,128,128],[s.c(1)-s.size(2)/2,s.c(2)-s.gap/2,s.c(1)+s.size(2)/2,s.c(2)+s.gap/2]);
        Screen('FillRect',display.windowPtr,[128,128,128],[s.c(1)-s.gap/2,s.c(2)-s.size(2)/2,s.c(1)+s.gap/2,s.c(2)+s.size(2)/2]);
        Screen('FillOval',display.windowPtr,[128,128,128],s.rect{1});
        Screen('FillOval',display.windowPtr,[255,255,255],s.rect{4});
    case 'off'
        Screen('FillOval',display.windowPtr,[100,100,100],s.rect{3});
        Screen('FillRect',display.windowPtr,[128,128,128],[s.c(1)-s.size(2)/2,s.c(2)-s.gap/2,s.c(1)+s.size(2)/2,s.c(2)+s.gap/2]);
        Screen('FillRect',display.windowPtr,[128,128,128],[s.c(1)-s.gap/2,s.c(2)-s.size(2)/2,s.c(1)+s.gap/2,s.c(2)+s.size(2)/2]);
        Screen('FillOval',display.windowPtr,[128,128,128],s.rect{1});
        Screen('FillOval',display.windowPtr,[255,255,255],s.rect{4});
    otherwise
        disp(sprintf('state %s not recognized',s.state));
end
end




