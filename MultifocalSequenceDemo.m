% MultifocalSequenceDemo.m


nAng = 24;  %number of segments in each ring
nEcc = 8;   %number of eccecentric rings of segments
nT = 200;  %number of 'blocks' (time steps) needs to be even

stimSeq = generateMultifocalSequence(nAng,nEcc,nT);


%%
% Visualize the sequence matrix

figure(1)
clf
imagesc(1-stimSeq);
colormap(gray)
axis equal;
axis tight
xlabel('Time');
ylabel('Segment number');

%%
% Correlations of time-courses between segments
%
% The constraint that neighboring segments are never on at the same time
% induces negative correlations between neighboring segments and positive
% correlations between segments within 'even' and 'odd' subsets.  This is
% not ideal.  This issue exists in the sequence that Omar and Geoff Aguirre
% sent us.

 c = corrcoef(stimSeq);
 
 % set the diagonals to NaN;
 c(find(eye(size(c)))) =NaN;
 
 figure(2)
 clf
 imagesc(c);
 colorbar
 colormap(gray);
 
  
 
 
 %%
 % autocorrelations

clear r

s = stimSeq(1,:);
shiftS = stimSeq(1,:);  % stimSeq(1,:) is the autocorrelation.  Other numbers give correlations between segments  
for i=1:length(s)
    tmp = corrcoef(s,shiftS);
    r(i) = tmp(1,2);
    shiftS = [shiftS(2:end),shiftS(1)];
end

figure(3)
clf
plot(r,'.-');
xlabel('shift')

%%
% Visualize the segments

[ecc,ang] = meshgrid(1:nEcc,1:nAng);
nSegments = nAng*nEcc;

segA = mod(ecc+mod(ang,2),2);  %checkerboard
segA = logical(segA(:));

[x,y] = meshgrid(linspace(-1,1,501));

ang = atan2(y,x);
ecc = sqrt(x.^2+y.^2);

angList = linspace(-pi,pi,nAng+1);
eccList = linspace(.1,1,nEcc+1);


%h= image(zeros(size(x)));
img = zeros(size(x));
count= 1;
clear wedgeId
for i = 1:nEcc
    for j=1:nAng
        count = count+1;
     
        img(ecc>eccList(i) &ecc<=eccList(i+1) & ang>angList(j) & ang<=angList(j+1)) = count;
    end
end

figure(4)
clf
h= image(x(1,:),y(:,1),img);
axis equal
axis tight
axis off

cmap=  hsv(nSegments+1);
cmap(1,:) = [0,0,0];
colormap(cmap);

colorbar


%%
% Show a movie of the sequence

frameRate = 4;  %Hz

colorbar off
tic
for i=1:nT
    img = zeros(size(x));
    cmap = zeros(nSegments+2,3);
    cmap(1,:) = [.5,.5,.5];
    for j=1:(nAng*nEcc)
        if stimSeq(j,i)
        %if segA(j);
            cmap(j+1,:) = [1,1,1];
        end
    end
    colormap(cmap);
    title(num2str(i));
    drawnow
    while(toc<i/frameRate)
    end
    
end




