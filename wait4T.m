function wait4T(triggerKey)

 KbName('UnifyKeyNames'); 
 
 
if ~exist('triggerKey','var')
    triggerKey = 't';
end


ch = '';
while ~strcmp(ch,triggerKey);
    [ keyIsDown, timeSecs, keyCode ] = KbCheck;
    keyPressed= KbName(keyCode);
    if ~isempty(keyPressed)
        ch = keyPressed(1);
    end
end

