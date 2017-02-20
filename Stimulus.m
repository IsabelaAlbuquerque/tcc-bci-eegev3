% Possible frequencies (from 6 to 30 Hz)
possFreq = zeros(9,3);
possFreq(:,1) = [5 5 4 4 3 3 2 2 1]';
possFreq(:,2) = [5 4 4 3 3 2 2 1 1]';
possFreq(:,3) = 60./(possFreq(:,1)+possFreq(:,2));
possFreq(2,3) = 6.7; % Needs finite number for comparison with the detection algorithm
possFreq(4,3) = 8.6;% For comparison with the detection algorithm

% Create the screen
whichScreen = 0;
window = Screen(whichScreen, 'OpenWindow');

% Get black/white colors
white = WhiteIndex(window);
black = BlackIndex(window);

% Create the rectangles
[winW, winH] = Screen('WindowSize',window);
rect = zeros(4,2);
rect(:,1) = [winW/15; winH/3; winW*(1/15+1/10); winH*2/3];
rect(:,2) = [winW-winW*(1/10+1/15); winH/3; winW-winW/15; winH*2/3];
%rectRightArrow = [winW/2+50,winH/2,winW/2+100,winH/2+50];
%rectLeftArrow = [winW/2-100,winH/2,winW/2-50,winH/2+50];

% Initialize the screen (completely black)
Screen(window, 'FillRect', black); 
Screen(window,'Flip');

% For rectangle flashing
fCount = zeros(2,1);

% Task parameters
refreshRate = 60; % In Hertz
freqInd = [2 4]; 
stimTime = 10000;     

% For rectangle flashing
fCount = zeros(2,1);
for i=1:stimTime*refreshRate
    
            for k = 1:length(freqInd) % For every frequency to be displayed
                if fCount(k) >= possFreq(freqInd(k),1) && fCount(k) < possFreq(freqInd(k),1)+possFreq(freqInd(k),2)
                    Screen(window, 'FillRect', black, rect(:,k));
                elseif fCount(k) >= 0 && fCount(k) < possFreq(freqInd(k),1)
                    Screen(window, 'FillRect', white, rect(:,k));
                end
                fCount(k) = fCount(k) + 1;
                if fCount(k) == possFreq(freqInd(k),1)+possFreq(freqInd(k),2)
                    fCount(k) = 0;
                end
            end          
            
            Screen(window,'Flip');
           

end

KbWait;
Screen('CloseAll');