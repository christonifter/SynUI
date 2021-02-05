syn = SynapseAPI(); 

SN = syn.getSamplingRates(); 
fs = SN.RZ6_1; 
dur = 60;
t = (1:dur*fs)./fs;
signal = sin(2*pi*2*t);
snipms = 30000;
npulses = dur*1E3/snipms;
snipsamps = round(30*fs);
snipsigs = reshape(signal, npulses, snipsamps);

syn.setMode(1);
% syn.setParameterValue('xDMR_audio', 'BufferLengthMS', snipms);
% syn.setParameterValue('xDMR_audio', 'TrainRepeats', npulses);
tic
syn.setParameterValues('xDMR_audio', 'buffer1', snipsigs(1,:));
syn.setParameterValues('xDMR_audio', 'buffer2', snipsigs(2,:));
toc
syn.setMode(2);
% r = rateControl(1E3/snipms);
% for i = 1:npulses
%     i
%     waitfor(r);
%     if mod(i, 2)
%         tic
%         toc
%     else
%         syn.setParameterValues('xDMR_audio', 'buffer2', snipsigs(i,:));
%     end    
% end

% syn.setMode(0);