initialize_experiment

%%

stimulators.left.port.BytesAvailable
[~] = fread(stimulators.left.port, stimulators.left.port.BytesAvailable); % Clear the readable bytes?

%%
stimulators.left.port.BytesAvailable
%%
startByte = 'FE';
commandLength = '02';
%commandBytes = ['01'; dec2hex(20)]; % Sets intensity to 20
commandBytes = ['01'; dec2hex(35)]; % Sets intensity to 20
endByte = 'FF';

checkSum = stimulators.left.calcCRC(commandBytes);
if length(checkSum)==1
    checkSum = strcat('0',checkSum);
end
sentControlCommand = hex2dec([startByte;commandLength;commandBytes;checkSum;endByte]); %Create the command line
fwrite(stimulators.left.port, sentControlCommand);
pause(0.5)
readData = fread(stimulators.left.port, stimulators.left.port.BytesAvailable)
%%
readData = fread(stimulators.left.port, 8);