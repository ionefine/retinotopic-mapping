# retinotopic-mapping
simple flexible retinotopy code

This retinotopic code was written to be extremely flexible, work on a number of scanners, and serve a number of purposes, these include
1. Including code for multifocal (Vanni Neuroimage 2005) as well as drifting bars
2. Code is written has having a background image windowed by an aperture. Allows you to have anything from a 
flickering checkerboard to a pixar movie in the background
3. Generates a binary mask, useful as input into prf mapping
4. Includes a very attentionally demanding Simon task
5. Easy to change the size/position of fixation spot
6. Can also create a 'fake patient' by a) jittering the fixation to match fixation instability 
(you can even read in a real patients eye-trace) and adds a scotoma. 

Written by Ione Fine and Geoff Boynton 2015-2016

****************************************************
MainRetinotopy.m is the main experiment


INSTRUCTIONS FOR THE SIMON TASK

Basically press the corresponding button (currently ‘1’ ’2’ ’3’ ’4’) to match the pattern of the colored segments. 
The sequence will get longer each time you get it right. If you get it wrong it will flash and then start again with a sequence of 1.
(If you have no idea what a Simon task is, look it up on the web before you run this, and it will make a lot more sense)

THINGS TO CHECK BEFORE A SESSION

Make sure the random seed isn’t fixed!!
Make sure the homedir is set correctly (line 9 in YorkRetintopy)
If you are running a real session, make sure the following are correct
display.dist = 36.5;     %distance from screen (cm)
display.width = 27.5; %width of screen (cm)
display.screenNum = 0;
stim.dur = 30; %duration of scan (seconds)

Within Wait4T make sure you are waiting for the right keycode from the scanner button box (default is ‘t’)
s.keys = {'w','s','a','q'}; Make sure the right keys in the button box are placed here for the Simon task

Other things to check before a scan
subjID='JMD1';
scanID='bar1'; 
stim.type='mf'; %also can be 'bar'

STIMULUS 

Checkerboard Background:  (you can replace this background with anything you want)
  can change flicker rate (stim.background.hz), 
  check sizes (in inside MakeSimulus.m: 
  stim.background.nsectors (number of segments in the multifocal)
  stim.background.nrings (number of rings in the multifocal)
  the number of phase-jitters for the checkboard edges aren’t always in the same place (altered in MakeStimulus: stim.background.n)

The retinotopic stimulus can be a drifting bar or multifocal
Drifting bars: most parameters (e.g. speed or width) are changed inside MakeStimulus
Multifocal: most parameters are changed inside MakeStimulus. 
One important flag is: stim.mf.noneighbors=0 which controls whether neighboring checks can be on, see note below.

SCOTOMA
stim.scotoma.on = 2; % 1: scotoma, 0: no scotoma, 2 is a scotoma matched to a real juvenile macular degeneration patient
Currently the scotoma is always modeled as a circle, with parameters center, (x, y in degrees) and radius (in degrees).

FIXATION
fix.offset, sets the fixation offset, x, y in degrees
s.scFac=2;% controls fixation Simon size

EYE-MOVEMENT
This is the weird bit. There are two parts of the eye-movement that are modelled separately.

eye.jitter.move=1; 
This jitters the stimulus (and scotoma) while keeping the fixation spot fixed, 
good for mimicking small rapid movements of the eyes due to fixation instability in patients. 

eye.drift.move=1; 
This moves the fixation spot (and scotoma) while keeping the stimulus fixed, good for mimicking larger drifts.

They can both be on simultaneously. For a real patient they should both be turned off.

NOTE FROM OMAR BUTT ABOUT USING NEIGHBORING CHECKS

So Butt (and I think Vanni) using multifocal had the constraint that neighboring checks could never be on. But, as seen by his email 
below, that was because they were trying to minimize lateral inhibition. If you aren't too worried about that, then
you may as well allow neightboring checks to be on.

"So to answer your question the short answer was it had to do with optimizing the detection power while 
(attempting to) minimize the effect of lateral inhibition. In our experiments, MF stimuli required 
stimuli presented much longer to provide equivalent maps in area V1. We thought this may be ameliorated 
by introducing constraints where adjacent regions could never be active at the same time. Even then, we 
still needed twice as much MF stimuli as bar to get the same map...and the maps themselves were limited to area V1 
(as also seen by Vanni). 

Looking back, while I think part of had to do with lateral inhibition, 
I think part also had to do with the unpredictable nature of the stimuli itself. 
With a traveling bar a subject can somewhat predict where a stimulus will be at a given time, 
even if your focusing a fixation task in the center. I very briefly experimented with this 
by 'dropping' one bar in my traveling bar stimulus while pretending it was still there in my analysis...
the results were almost identical (unpublished, very pilot data, there of course could be a huge number of confounds here). 
I think the proper way to tease this out will be to use a traveling bar stimulus which has the bars randomly 
appear on the screen vs using a fixed Left to right, then diagonals or whatever pattern, and see if there is a 
qualitative difference in the produced maps or quantitative difference in the average beta values. 

I have code for the stimuli, but it was all based on Vanni's toolbox
And you are right, I could have used a different sequence...at the time that was all I had (early/mid 2010)

This was (tangentially) discussed on my 2011 Poster: https://cfn.upenn.edu/aguirre/wiki/_media/public:presentations:vss2011_barvsmf.pdf"

