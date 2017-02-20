function varargout = calibtask2(varargin)
% CALIBTASK2 MATLAB code for calibtask2.fig
%      CALIBTASK2, by itself, creates a new CALIBTASK2 or raises the existing
%      singleton*.
%
%      H = CALIBTASK2 returns the handle to a new CALIBTASK2 or the handle to
%      the existing singleton*.
%
%      CALIBTASK2('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in CALIBTASK2.M with the given input arguments.
%
%      CALIBTASK2('Property','Value',...) creates a new CALIBTASK2 or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before calibtask2_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to calibtask2_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help calibtask2

% Last Modified by GUIDE v2.5 03-Jul-2014 13:49:19

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                  'gui_Singleton',  gui_Singleton, ...
                  'gui_OpeningFcn', @calibtask2_OpeningFcn, ...
                  'gui_OutputFcn',  @calibtask2_OutputFcn, ...
                  'gui_LayoutFcn',  [] , ...
                  'gui_Callback',   []);
if nargin && ischar(varargin{1})
   gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
   [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
   gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before calibtask2 is made visible.
function calibtask2_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to calibtask2 (see VARARGIN)

% Choose default command line output for calibtask2
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes calibtask2 wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = calibtask2_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


function edit7_Callback(hObject, eventdata, handles)    % Get the IP address. 
% hObject    handle to edit7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit7 as text
%        str2double(get(hObject,'String')) returns contents of edit7 as a double

handles.ipadress = get(hObject,'String');   % Get the IP address and save on handles structure
guidata(hObject, handles);                  % Uptade handles structure


% --- Executes during object creation, after setting all properties.
function edit7_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
   set(hObject,'BackgroundColor','white');
end

% --- Executes on button press in pushbutton1.
function pushbutton1_Callback(hObject, eventdata, handles)  % Button responsible for the connection with Enobio
% hObject    handle to pushbutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

host = handles.ipadress; % ip address
%Creating the connection
message = sprintf('Press YES on NICs dialog box');
nicdialog = msgbox(message,' ','warn');
[ret, status, socket] = MatNICConnect(host);    % Creates the connection
%Start EEG
ret = MatNICStartEEG(socket);                   % Start streaming EEG data
handles.socket = socket;                        % Save the socket number on handles
guidata(hObject, handles);                      % Uptade handles structure  
set(handles.pushbutton1,'Enable','off');        % Disable the button


% --- Executes on button press in pushbutton2.
function pushbutton2_Callback(hObject, eventdata, handles)  % Button responsible for the connection with the robot
% hObject    handle to pushbutton2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


%robot = Brick('ioType','usb'); % Usb connection. Just necessary when the bluetooth connection doesn't work
robot = Brick('ioType','instrbt','btDevice','EV3','btChannel',1);   % Creates a "Brick" object connected with Matlab by bluetooth
handles.robot = robot;
%b=Brick('ioType','bt','serPort','/dev/rfcomm0'); %bluetooth connection
robot.beep();       % One beep to show that the connection is working properly

% setting motor power
robot.outputPower(0,Device.MotorB,10);
robot.outputPower(0,Device.MotorD,10);
robot.outputStart(0,Device.MotorB);         % The robot moves forward for a while
robot.outputStart(0,Device.MotorD);
pause(1)
robot.outputStopAll()
set(handles.pushbutton2,'Enable','off');    % Disable the button
guidata(hObject, handles);                  % Uptade handles structure   


function edit1_Callback(hObject, eventdata, handles)    % Get the number of seconds to look to each rectangle
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit1 as text
%        str2double(get(hObject,'String')) returns contents of edit1 as a double

secr = str2double(get(hObject,'String'));   % Get the number of seconds and save on secr
handles.secr = secr;                        % Save secr on handles
guidata(hObject, handles);                  % Uptade handles structure


% --- Executes during object creation, after setting all properties.
function edit1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
   set(hObject,'BackgroundColor','white');
end



function edit2_Callback(hObject, eventdata, handles)    % Get the number of seconds to look to baseline
% hObject    handle to edit2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit2 as text
%        str2double(get(hObject,'String')) returns contents of edit2 as a double

secb = str2double(get(hObject,'String'));    % Get the number of seconds and save on secb
handles.secb = secb;                         % Save secb on handles   
guidata(hObject, handles);                   % Uptade handles structure



% --- Executes during object creation, after setting all properties.
function edit2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
   set(hObject,'BackgroundColor','white');
end



function edit3_Callback(hObject, eventdata, handles)    % Get the number of calibrations
% hObject    handle to edit3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit3 as text
%        str2double(get(hObject,'String')) returns contents of edit3 as a double

nc = str2double(get(hObject,'String'));       % Get the number of calibrations and save on nc
handles.nc = nc;                              % Save nc on handles   
guidata(hObject, handles);                    % Uptade handles structure  


% --- Executes during object creation, after setting all properties.
function edit3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
   set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton3.
function pushbutton3_Callback(hObject, eventdata, handles)  % Calibration task
% hObject    handle to pushbutton3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


message = sprintf('When you listen: \n One long beep: look to the rectangle on the left \n Two beeps: look to the dark screen \n Three beeps: look to the rectangle on the right \n One short beep: recording has started', 'ATTENTION!','warn');
attention = msgbox(message,'ATTENTION!','warn');
pause(10);

secr = handles.secr;
secb = handles.secb;
nc = handles.nc;
host = handles.ipadress;

   dataf1 = zeros(500*secr*nc,8);   % Sampling rate = 500; 8 channels
   baseline = zeros(500*secb*nc,8);
   dataf2 = zeros(500*secr*nc,8);

   %Beep creation
   fs = 20000;
   ts = 1/fs;
   t1 = 0:ts:1;
   beep1 = sin(2*pi*330*t1);
   t2 = 0:ts:0.5;
   beep2 = sin(2*pi*330*t2);
   t3 = 0:ts:0.333;
   beep3 = sin(2*pi*330*t3);

   for i = 1:nc
       pause(5);
       sound(beep1,fs);
       pause(3);        
       sound(beep3,fs);
       dataf1(1+(500*secr)*(i-1):(500*secr)*i,:) = MatNICEEGRecord(secr,8,host); % The user is looking to the rectangle on the left and EEG signal is being recorded
       sound(beep2,fs);
       sound(beep2,fs);
       pause(7);
       sound(beep3,fs);
       baseline(1+(500*secb)*(i-1):(500*secb)*i,:) = MatNICEEGRecord(secb,8,host); % The user is looking to the dark screen and EEG signal is being recorded
       sound(beep3,fs);
       sound(beep3,fs);
       sound(beep3,fs);
       pause(7);
       sound(beep3,fs);
       dataf2(1+(500*secr)*(i-1):(500*secr)*i,:) = MatNICEEGRecord(secr,8,host); % The user is looking to the rectangle on the right and EEG signal is being recorded
       sound(beep3,2*fs);
   end;   

   %Training and testing the classifier

   [W1,W2,W3,acc1,acc2,acc3,acc11,acc22,acc33] = LDA4(dataf1,baseline,dataf2); % Calculate the LDA weights using the data acquired before

% Saving the results on handles

handles.W1 = W1;    
handles.W2 = W2;
handles.W3 = W3;
handles.acc = [acc1 acc11; acc2 acc22; acc3 acc33];
handles.dataf1 = dataf1;
handles.dataf2 = dataf2;
handles.baseline = baseline;
guidata(hObject, handles);  % Update handles structure

t = handles.uitable1;
accuracy = handles.acc;
set(t,'Data', accuracy);    % Write the results on the table


function edit6_Callback(hObject, eventdata, handles) % Get the number of seconds to control the robot
% hObject    handle to edit6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit6 as text
%        str2double(get(hObject,'String')) returns contents of edit6 as a double

period = str2double(get(hObject,'String')); % Get the number of seconds to control the robot and save on period
handles.period = period;                    % Save period on handles
guidata(hObject, handles);                  % Update handles structure


% --- Executes during object creation, after setting all properties.
function edit6_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
   set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton4.
function pushbutton4_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Setting the initial parameters
period = handles.period;    % recording period (seconds)
W1 = handles.W1;
W2 = handles.W2;
W3 = handles.W3;
host = handles.ipadress;
robot = handles.robot;
socket = handles.socket;
%motorPower = handles.motorPower;
n_channel = 8;  %number of channels

nsecond = 2; % period of data used to control the robot
maxi = 2*nsecond; %maximum value of cont, that determines the number of samples on the input matrix
   
ip   = host;       % default to home ip address
p    = 1;          % duration of the pause before acquisition (seconds, min 1s)
N    = n_channel;  % number of channels
bps  = 4;          % bytes per sample
Ns   = 250;        % number of samples to read in the buffer each time (between 20 and 250)
buff = bps*Ns*N;   % number of bytes to read in total for all channels
sps = 500;         % number of samples per second
collectN = ceil(period*sps/Ns);  %number of biffers to attain collection time


% Create TCP/IP connection
import java.net.Socket
import java.io.*
import java.nio.*
sock = Socket(ip,1234);
sock.setSoTimeout(10000); % Wait 10 seconds for a response
ch = channels.Channels.newChannel(sock.getInputStream);
bytes = ByteBuffer.allocate(buff);
bytes.order(ByteOrder.BIG_ENDIAN);

%disp('Please wait...')
pause(p);

eeg = zeros(N,collectN*Ns); % All the data recorded during the cotrolling task
   
signal = zeros(maxi*Ns,N);  % Data used for take one decision
n_bytes = 0;
counter = 0;
cont = 0;
handles.decision = {};  % Save the decisions 
% To skip the problem with the motors (see documentation)
robot.outputPower(0,Device.MotorB,10);
robot.outputPower(0,Device.MotorD,10);
robot.outputStart(0,Device.MotorB);
robot.outputStart(0,Device.MotorD);
pause(0.1)
robot.outputStopAll()
pause(0.1)
% Start to move the robot forward
robot.outputStart(0,Device.MotorB);
robot.outputStart(0,Device.MotorD);
while 1
   n_bytes = n_bytes + ch.read(bytes);
if n_bytes < buff % if the buffer is not completely full, keep streaming data
continue;
end
n_bytes = 0;    % if the buffer is completely full, record the data on the output matrix
bytes.rewind;
   for k = 1 : Ns
       for j = 1 : N
           eeg(j,k+counter*Ns) = bytes.getInt;    
       end;
   end;
   eeg1 = eeg'; 
   if cont ~= maxi  %input is not completely full
       for k = 1 : Ns
           for j = 1 : N
               signal(k+cont*Ns,j) = eeg1(k+counter*Ns,j);
           end;
       end;    
       cont = cont+1;
   else
       % Decision
       n = length(signal); 
       L = 2;
       nw = n/(500*L);
       Fs = 500;
       % f1
       a=round(6.4*500*L/Fs);  
       b=round(7*500*L/Fs);
       c=round(13.1*500*L/Fs);
       d=round(13.7*500*L/Fs);
       e=round(8.3*500*L/Fs);
       f=round(8.9*500*L/Fs);
       g=round(16.9*500*L/Fs);
       h=round(17.5*500*L/Fs);
       % Calculate the features for each channel
       % Channel 1   
       dec_1 = zeros(nw,4);
       for i = 1:nw
           x = signal((i-1)*500*L+1:i*(500*L),1);
           X = fft(x);
           absX = abs(X).^2;
           dec_1(i,1) = mean(absX(a:b));
           dec_1(i,2) = mean(absX(c:d));
           dec_1(i,3) = mean(absX(e:f));
           dec_1(i,4) = mean(absX(g:h));
       end;
       % Channel 3
       dec_3 = zeros(nw,4);
       for i = 1:nw
           x = signal((i-1)*500*L+1:i*(500*L),3);
           X = fft(x);
           absX = abs(X).^2;
           dec_3(i,1) = mean(absX(a:b));
           dec_3(i,2) = mean(absX(c:d));
           dec_3(i,3) = mean(absX(e:f));
           dec_3(i,4) = mean(absX(g:h));
       end;

       dec=[dec_1 dec_3];
       
       % f1 x baseline
       L111 = [zeros(1,1) dec] *W1';
       P111= exp(L111)./repmat(sum(exp(L111),2),[1 2]);
       [~,maxInd1] = max(P111,[],2);
       m1 = mode(maxInd1(1:length(maxInd1),:));
       
       % f2 x baseline
       L222 = [zeros(1,1) dec] *W2';
       P222 = exp(L222)./repmat(sum(exp(L222),2),[1 2]);
       [~,maxInd2] = max(P222,[],2);
       m2 = mode(maxInd2(1:length(maxInd2),:));
       
       % f2 x f1
       L333 = [zeros(1,1) dec] *W3';
       P333 = exp(L333)./repmat(sum(exp(L333),2),[1 2]);
       [~,maxInd3] = max(P333,[],2);
       m3 = mode(maxInd3(1:length(maxInd3),:));
       
       if m1 == 1 && m2 == 1
           decision = 0; % Baseline
       elseif (m1 == 1 && m2 == 2) ||  (m1 == 2 && m2 == 1) || (m1 == 2 && m2 == 2)
           if m3 == 1
               decision = 1; % f1
           else    
               decision = 2; % f2
           end;   
       end;
       guidata(hObject, handles);
       %decision = Decision4(signal,W1,W2,W3); %decision
       if decision == 2    % moving the robot to the right
           display('Turn to the right');    
           robot.outputPower(0,Device.MotorB,20);
           pause(2)
           robot.outputPower(0,Device.MotorB,10);
           handles.decision{end+1} = 'Turn to the right';
       elseif decision == 1 %moving the robot to the left
           display('Turn to the left');   
           robot.outputPower(0,Device.MotorD,20);
           pause(2)
           robot.outputPower(0,Device.MotorD,10);
           handles.decision{end+1} = 'Turn to the left';
       end;
       cont = 0;
       guidata(hObject, handles);
       decisionst = handles.decision;               
       set(handles.listbox2,'String', decisionst);  % Refresh the output listbox
        
  end     
  counter = counter+1;
  bytes.rewind;
  if counter == collectN
       break 
  end
end
%disp('Done.')
robot.outputStopAll();
dataf1 = handles.dataf1;
dataf2 = handles.dataf2;
baseline = handles.baseline;

save('data6.mat','W1','W2','W3','eeg1','dataf1','dataf2','baseline');   % Saving all the data acquired

%Clean up the connection
sock.close;
%disp('Connection closed.')
robot.outputStopAll();
%set(handles.pushbutton1,'Enable','on');
%ret = MatNICStopEEG(socket);


% --- Executes on selection change in listbox2.
function listbox2_Callback(hObject, eventdata, handles)
% hObject    handle to listbox2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns listbox2 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from listbox2
decisionst = handles.decision;
set(hObject,'String', decisionst);

% --- Executes during object creation, after setting all properties.
function listbox2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to listbox2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
   set(hObject,'BackgroundColor','white');
end

% --- Executes on button press in pushbutton8.
function pushbutton8_Callback(hObject, eventdata, handles)     % Instructions button
% hObject    handle to pushbutton8 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
message = sprintf(' 1-Set headband and electrodes on your head \n 2-Turn on Enobio \n 3-Plug robots USB cable \n 4-Turn on the robot \n 5-Run Stimulus.m on another computer \n 6-Open NIC \n 7-The IP adress is on the right \n bottom corner on NICs window');
instruc = msgbox(message,'Instructions','warn');

