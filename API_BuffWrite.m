close all; clear all; clc; 

% Call the API, get important sampling rates, set mode

syn = SynapseAPI(); 

SN = syn.getSamplingRates(); 
SN_RZ6 = SN.RZ6_1; %rz6 sampling rate for later caluclations 

Preview = 1; %set to 0 if you want to record

USER_GIZMO = 'BuffMstr1'; %name of your buffer gizmo


%% Make some signals for later

seconds = 1; %make a 1 second tone
fs = SN_RZ6; %your system's sampling rate for correct output

tone_length = seconds*fs; %length (in samples) of tone based on desired time
t = (0:1/fs:seconds); %time vector
f1 = 10; %sin freq
f2 = 100; 
f3 = 1000; 

sin1 = sin(2*pi*f1*t); %1 V p-p sin wave
sin2 = sin(2*pi*f2*t); %1 V p-p sin wave
sin3 = sin(2*pi*f3*t); %1 V p-p sin wave

%% Runtime

if syn.getMode() < 2
   if Preview
    syn.setModeStr('Preview');  
   else
    syn.setModeStr('Record'); 
   end
end

pause(5); %wait 

syn.setParameterValue(USER_GIZMO,'SignalLength',numel(sin1)); %set the schmitt2 length
syn.setParameterValues(USER_GIZMO,'signal',sin1); %write the array. Notice the 's' in Values
syn.setParameterValue('PulseGen1','Enable',1); %Trigger pulse gen

pause(20); %obviously this is a brutish way of doing all this

syn.setParameterValue(USER_GIZMO,'SignalLength',numel(sin2)); %set the schmitt2 length
syn.setParameterValues(USER_GIZMO,'signal',sin2); %write the array. Notice the 's' in Values
syn.setParameterValue('PulseGen1','Enable',1); %Trigger pulse gen

pause(20)

syn.setModeStr('Idle'); 


