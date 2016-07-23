function stimSeq = generateMultifocalSequence(nAng,nEcc,nT,neighborFlag)
% stimSeq = generateMultifocalSequence(nAng,nEcc,nT)
%
% Generates a matrix that determines which segments are 'on' during a
% multifocal sequence in which neighboring segments are never on at the
% same time. Segments are defined to be arc-shaped regions and are numbered
% clockwise starting from the innermost eccentricity.
%
% Inputs:
%   nAng   number of segments for each eccentricity
%   nEcc   number of rings of segments
%   nT     number of temporal frames in the sequence
%   neigborFlag flag to determine if neighbors cannot be on at the same time
%   (default: true, neighbors are never simultaneously on)
%
% Outputs:
%  stimSeq a [nSegments x nT] size binary matrix, where nSegments = nAng*nEcc
%
% Segments are divided into 'even' and 'odd' complementary checkerboard
% patterns.  For any time point, only the even or the odd segments will be
% set to 'on', which prevents neighboring segments to be on at the same
% time.
%
% Two m-squences are used to determine the time-series of the 'even' and
% 'odd' segments.  At any TR, a third m-sequence determines whether the
% even or odd sequence is on.
%
% Written by G.M. Boynton, summer of 2005

if ~exist('neighborFlag','var')
    neighborFlag = true;
end

primeNum = 11;  %can be any prime number
    nSegments = nAng*nEcc;


%%
% Define segA, a logical vector that is an index into the 'odd' segments.
% ~segA indexes the even segments.


if neighborFlag
    
    [ecc,ang] = meshgrid(1:nEcc,1:nAng);

    segA = mod(ecc+mod(ang,2),2);   % checkerboard
    segA = logical(segA(:));        % unwrap into a vector
    
    % Define the m-sequences for odd (seqA) and even (seqB) segments
    
    pow = ceil(log2((nT+2)/2)); % power value for the shortest possible m-sequence
    seqA = mseq(2,pow,1,1)';    % 'odd' segments
    seqB = mseq(2,pow,1,2)';    % 'even' segments
    % Segments within even or odd groups are circularly shifted versions of the
    % same m-sequence.  The shift for each successive segment is a prime
    % number, modulo the number of segments within the group (nSegments/2). See
    % Vanni et al. NeuroImage (2005) p.98.
    
    shift = mod([0:(nSegments/2-1)]*primeNum,nSegments/2);
    
    % Define the sequence matrix for 'even' and 'odd' segments
    
    seqAmat = zeros(nSegments/2,length(seqA));
    seqBmat = zeros(nSegments/2,length(seqB));
    for i=1:(nSegments/2)
        segAmat(i,:) = circshift(seqA,shift(i),2);
        segBmat(i,:) = circshift(seqB,shift(i),2);
    end
    
    % Get the m-sequence that determines which segments are applied for each
    % time-frame
    segSwitch = mseq(2,pow+1,1,1); % segSwitch defines whether even or odd is applied for each time-point
    segSwitch = logical(segSwitch(1:nT));  %truncate to length nT
    nA = sum(segSwitch==1); %number of time-frames that point to 'odd' segments
    
    % Stick the sequences for the even and odd segments into the overall
    % sequence
    
    stimSeq = zeros(nSegments,nT);
    stimSeq(segA,segSwitch) = segAmat(:,1:nA); % stick in the odd segments
    stimSeq(~segA,~segSwitch) = segBmat(:,1:(nT-nA));  % even segments
    
else
    pow = ceil(log2((nT+1))); % power value for the shortest possible m-sequence
    seq1 = mseq(2,pow,1,1)';    % 'odd' segments
        seq2 = mseq(2,pow,1,2)';    % 'odd' segments

    seq = seq1*2+seq2;
    seq = (seq(1:nT) == 0);
        
    shift = mod([0:(nSegments-1)]*primeNum,nSegments);
    stimSeq = zeros(nSegments,nT);
    
    for i=1:nSegments
        stimSeq(i,:) = circshift(seq,shift(i),2);
    end
end









