function mask = getBarMask(display,stim)

n = display.resolution(2);

[x,y] = meshgrid(linspace(-stim.maxRad,stim.maxRad,n));
 

r = sqrt(x.^2+y.^2);
nMasks = stim.dur/(1./stim.hz);
mask = false([size(x),nMasks]);


for frameNum = 1:nMasks
    a = cos(stim.ang(frameNum)*pi/180)*(x-stim.center(frameNum,1)) + sin(stim.ang(frameNum)*pi/180)*(y-stim.center(frameNum,2));
    b = -sin(stim.ang(frameNum)*pi/180)*(x-stim.center(frameNum,1)) +cos(stim.ang(frameNum)*pi/180)*(y-stim.center(frameNum,2));  
    mask(:,:,frameNum) = (abs(b)<stim.width/2 & r<stim.maxRad);
end




