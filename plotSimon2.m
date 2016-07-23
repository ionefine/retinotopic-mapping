function plotSimon2(s)

tmp=[s.event(:).time];
tmp=tmp(tmp~=1000);
dur = max(tmp);
tsample=1000; 
figure(1)
clf
hold on
mat=zeros(1, round(dur*tsample));

ct=1;
status=s.event(1).type;
til=s.event(2).time;
s.event(end+1).time=1000; 
for t=1:size(mat, 2); 
    if strcmp(status, 'play: on');
    mat(1,t)=1;
    elseif strcmp(status, 'recall: on')
        mat(1,t)=2;
     elseif strcmp(status, 'error: start')
        mat(1,t)=3;
    end
    if t/tsample>til % event change
        ct=ct+1;
        status=s.event(ct).type;
        til=s.event(ct+1).time;
    end
    
end

% image it
colormap([0 0 0 ; 0 1 0 ; 0 0 1; 1 0 0 ; .5 .5 .5])
imagesc(mat); 
axis tight
%axis off
