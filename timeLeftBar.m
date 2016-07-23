function timeLeftBar(h,i,n,str)
% timeLeftBar(h,i,n,[str])

if ~exist('str','var')
    str = 'Time left:';
end

timeLeft = toc*(n-i)/i;  % seconds

 
if timeLeft<60
    timeLeftStr = sprintf('%d sec',floor(timeLeft));
elseif timeLeft <60*60
    timeLeftStr = sprintf('%d min %d sec',floor(timeLeft/60),floor(mod(timeLeft,60)));
else
    timeLeftStr = sprintf('%d hr %d min',floor(timeLeft/3600),...
        floor(mod(timeLeft,3600)));
end

waitbar(i/n,h,sprintf('%s %s',str,timeLeftStr));
if i==n
    delete(h)
end
