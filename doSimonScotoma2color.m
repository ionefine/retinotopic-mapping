
function  task=doSimonScotoma2color(display,task,curTime)
%s=doSimon(display,s,curTime)
escPressed=0;

switch task.action
    case 'init'  %set up simon variables
        % 'simon' parameters
        if ~isfield(task, 'scFac')
            task.scFac=1;
        end
        
        task.size = angle2pix(display, task.apertures);  %pixels
        task.gap = 2; %pixels
        task.c = ceil(display.resolution/2)+task.fixoffset;
        for i=1:length(task.size)
            task.rect{i} = [-task.size(i)/2+task.c(1),-task.size(i)/2+task.c(2),+task.size(i)/2+task.c(1),+task.size(i)/2+task.c(2)];
        end
                
        task.color{1} = [0,0,255];
        task.color{2} = [255,255,0];
        task.color{3} = [0,255,0];
        task.color{4} = [255,0,0];

        task.seq = [];
        task.maxLen = [];
        task.i = 0;
        
        task.state = 'off';
        
        task.seq = task.goodKeys(ceil(rand(1)*length(task.goodKeys)));
        task.pauseTill = 1;
        task.action = 'play';
        
        task.event.type = 'switch to play';
        task.event.time = 0;
        task.event.num = 0;
        task.numEvents = 1;
        
        showSimon(display,task);
        
    case 'play'  %sequence playback mode
        if task.pauseTill < curTime
            t = mod(curTime,task.ISI);
            
            if strcmp(task.state,'off') && t<=task.dur
                
                task.i = task.i+1;
                
                if task.i<=length(task.seq);
                    task.show = task.seq(task.i);
                    task.state = 'on';
                    
                    task.numEvents = task.numEvents+1;
                    task.event(task.numEvents).type= 'play: on';
                    task.event(task.numEvents).time = curTime;
                    task.event(task.numEvents).num = task.show;
                    
                else
                    task.numEvents = task.numEvents+1;
                    task.event(task.numEvents).type= 'switch to recall';
                    task.event(task.numEvents).time = curTime;
                    task.event(task.numEvents).num = NaN;
                    
                    task.action = 'recall';
                    task.state = 'off';
                    task.i = 0;
                end
                
            elseif strcmp(task.state,'on') && t> task.dur
                task.numEvents = task.numEvents+1;
                task.event(task.numEvents).type= 'play: off';
                task.event(task.numEvents).time = curTime;
                task.event(task.numEvents).num = task.show;
                task.state = 'off';
                
            end
        end
        showSimon(display,task);
        
    case 'recall'
        if task.pauseTill<curTime
            
            %see if a key is pressed
            [ keyIsDown, timeSecs, keyCode ] = KbCheck;
            
            if keyIsDown && strcmp(task.state,'off') %a key is down: record the key and time pressed
                
                keyPressed= KbName(keyCode);
                
                keyNum = find(strcmp(keyPressed(1),task.keys));
                
                if ~isempty(keyNum)
                    task.show = keyNum;
                    task.state = 'on';
                    task.numEvents = task.numEvents+1;
                    task.event(task.numEvents).type= 'recall: on';
                    task.event(task.numEvents).time = curTime;
                    task.event(task.numEvents).num = task.show;
                end
                
            elseif ~keyIsDown && strcmp(task.state,'on') %key just lifted
                task.state = 'off';
                task.numEvents = task.numEvents+1;
                task.event(task.numEvents).type= 'recall: off';
                task.event(task.numEvents).time = curTime;
                task.event(task.numEvents).num = NaN;
                
                task.i = task.i+1;
                if task.show ~= task.seq(task.i)
                    %incorrect
                    task.maxLen = [task.maxLen,length(task.seq)-1];
                    task.pauseTill = ceil(curTime+ task.errDur);
                    task.action = 'error';
                    
                    task.numEvents = task.numEvents+1;
                    task.event(task.numEvents).type= 'error: start';
                    task.event(task.numEvents).time = curTime;
                    task.event(task.numEvents).num = task.seq(end);
                    
                elseif task.i == length(task.seq)  %got the last one right
                    task.seq = [task.seq,task.goodKeys(ceil(rand(1)*length(task.goodKeys)))];
                    task.state = 'off';
                    task.action = 'play';
                    task.i = 0;
                    task.pauseTill = ceil(curTime+ task.pauseDur) ;
                    
                    task.numEvents = task.numEvents+1;
                    task.event(task.numEvents).type= 'switch to play';
                    task.event(task.numEvents).time = curTime;
                    task.event(task.numEvents).num = NaN;
                end
            end
        end
        showSimon(display,task);
    case 'error'
        if curTime < task.pauseTill
            ab = mod(curTime*task.errFreq,1);
            task.show = task.seq(end);
            if ab<.5
                task.state = 'on';
            else
                task.state = 'off';
            end
            showSimon(display,task);
        else
            task.i = 0;
            task.seq = task.goodKeys(ceil(rand(1)*length(task.goodKeys)));
            task.action = 'play';
            task.state = 'off';
            task.pauseTill = ceil(curTime+ task.pauseDur);
            task.numEvents = task.numEvents+1;
            task.event(task.numEvents).type= 'error: stop';
            task.event(task.numEvents).time = curTime;
            task.event(task.numEvents).num = NaN;
        end
    case 'done'
        task.maxLen = [task.maxLen,length(task.seq)-1];
        task.numEvents = task.numEvents+1;
        task.event(task.numEvents).type= 'done';
        task.event(task.numEvents).time = curTime;
        task.event(task.numEvents).num = NaN;
        task.seq = [];
end
end

function showSimon(display,task)
task.c = ceil(display.resolution/2)+task.fixoffset;
for i=1:length(task.size)
    task.rect{i} = [-task.size(i)/2+task.c(1),-task.size(i)/2+task.c(2),+task.size(i)/2+task.c(1),+task.size(i)/2+task.c(2)];
end
switch task.state
    case 'on';
        Screen('FillOval',display.windowPtr,[100,100,100],task.rect{3});
        if strcmp(task.action, 'error')
            Screen('FillArc',display.windowPtr,[ 0 0 0 ],task.rect{2},180*(task.show-1)+task.rotAng,180);
        else
            Screen('FillArc',display.windowPtr,task.color{task.show},task.rect{2},180*(task.show-1)+task.rotAng,180);
        end
        
        Screen('FillOval',display.windowPtr,[128,128,128],task.rect{1});
        Screen('FillOval',display.windowPtr,[255,255,255],task.rect{4});
        
        if task.rotAng==0
            Screen('FillRect',display.windowPtr,[128,128,128],[task.c(1)-task.gap/2,task.c(2)-task.size(2)/2,task.c(1)+task.gap/2,task.c(2)+task.size(2)/2]);
        elseif task.rotAng==-45
            Screen('DrawLine',display.windowPtr, [128,128,128], ...
                task.c(1)-(.75*task.size(2))/2, task.c(2)-(.75*task.size(2))/2, ...
                task.c(1)+(.75*task.size(2))/2 , task.c(2)+(.75*task.size(2))/2, min([task.scFac 4]));
        elseif task.rotAng==45
            Screen('DrawLine',display.windowPtr, [128,128,128], ...
                task.c(1)-(.75*task.size(2))/2, task.c(2)+(.75*task.size(2))/2 , ...
                task.c(1)+(.75*task.size(2))/2 , task.c(2)-(.75*task.size(2))/2, min([task.scFac 4]));          
        elseif task.rotAng==90
            Screen('FillRect',display.windowPtr,[128,128,128],[task.c(1)-task.size(2)/2,task.c(2)-task.gap/2,task.c(1)+task.size(2)/2,task.c(2)+task.gap/2]);
        end
    case 'off'
        Screen('FillOval',display.windowPtr,[100,100,100],task.rect{3});
        Screen('FillOval',display.windowPtr,[128,128,128],task.rect{1});
        Screen('FillOval',display.windowPtr,[255,255,255],task.rect{4});
        
         if task.rotAng==0
            Screen('FillRect',display.windowPtr,[128,128,128],[task.c(1)-task.gap/2,task.c(2)-task.size(2)/2,task.c(1)+task.gap/2,task.c(2)+task.size(2)/2]);
        elseif task.rotAng==-45
            Screen('DrawLine',display.windowPtr, [128,128,128], ...
                task.c(1)-(.75*task.size(2))/2, task.c(2)-(.75*task.size(2))/2, ...
                task.c(1)+(.75*task.size(2))/2 , task.c(2)+(.75*task.size(2))/2, min([task.scFac 4]));
        elseif task.rotAng==45
            Screen('DrawLine',display.windowPtr, [128,128,128], ...
                task.c(1)-(.75*task.size(2))/2, task.c(2)+(.75*task.size(2))/2 , ...
                task.c(1)+(.75*task.size(2))/2 , task.c(2)-(.75*task.size(2))/2, min([task.scFac 4]));
            
        elseif task.rotAng==90
            Screen('FillRect',display.windowPtr,[128,128,128],[task.c(1)-task.size(2)/2,task.c(2)-task.gap/2,task.c(1)+task.size(2)/2,task.c(2)+task.gap/2]);
        end
        
    otherwise
        disp(sprintf('state %s not recognized',task.state));
end
end




