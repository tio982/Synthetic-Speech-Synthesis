% Load MP3 file and check that it plays
[audioData, sampleRate] = audioread("had_f.wav");
% [audioData, sampleRate] = audioread("had_m.wav");
% sound(audioData, sampleRate);

% % Segment Isolation male
% segmentStart = 0.03; 
% segmentLength = 0.2; 

% Segment Isolation female
segmentStart = 0.05; 
segmentLength = 0.2; 

startTime = ceil(segmentStart * sampleRate); 
endTime = floor((segmentStart + segmentLength)* sampleRate);

endTime = min(endTime, length(audioData));
isolatedSegment = audioData(startTime:endTime);
isolatedSegment = isolatedSegment.*hamming(length(isolatedSegment));

%Pre-emphasis filter
preEmphasis = [1 -0.95];
% preEmphasis = [1 -0.97]; male
isolatedSegment = filter(preEmphasis, 1, isolatedSegment);

lpcOrder = 25;
[lpcCoefficients, ~] = lpc(isolatedSegment, lpcOrder);
processedSpeech = filter(1, lpcCoefficients, isolatedSegment);

%Formant estimates using roots of lpc 
rts = roots(lpcCoefficients);
rts = rts(imag(rts)>= 0);
rts = rts(abs(abs(rts)-1) < 0.06);
% rts = rts(abs(abs(rts)-1) < 0.03); male
angz = atan2(imag(rts), real(rts));
formantFrequencies = angz * (sampleRate/(2*pi));
sortedFormantFrequencies = sort(formantFrequencies);
% formantFrequencies = sortedFormantFrequencies(sortedFormantFrequencies > 350);
% male
formantFrequencies = sortedFormantFrequencies(sortedFormantFrequencies > 100); 

format short g;
disp(sortedFormantFrequencies);


[autoCorrelationValues, autoCorrelationLags] = xcorr(isolatedSegment, "coeff");

minLag = ceil(sampleRate/ 270);
% minLag = ceil(sampleRate/ 185); male


zeroLagIndex = find(autoCorrelationLags==0);

[peak,peakIndex] = max(autoCorrelationValues(zeroLagIndex + minLag:end));

peakLag = autoCorrelationLags(zeroLagIndex + minLag -1 + peakIndex);

if peakLag > 0
    meanFundamentalFrequency = sampleRate / peakLag;
else
    meanFundamentalFrequency = NaN;
end
disp('Mean Fundamental Frequency in Hz:');
disp(meanFundamentalFrequency);


f = meanFundamentalFrequency;
t = (0:length(isolatedSegment)-1) / sampleRate;
impulseSignal = zeros(size(t));
impulseInterval = round(sampleRate / f);
impulseSignal(1:impulseInterval:end) = 1;


windowLength = round(0.025 * sampleRate);
hammingWindow = hamming(windowLength);
numPulses = floor((length(impulseSignal) - windowLength) / impulseInterval) + 1;
windowedImpulse = zeros(size(impulseSignal));



for i = 1:numPulses
    startTime = (i-1) * impulseInterval + 1;
    endTime = startTime + windowLength - 1;
    if endTime > length(impulseSignal)
        endTime = length(impulseSignal);
    end
    currentWindowLength = endTime - startTime + 1;
    if currentWindowLength == windowLength
        currentWindow = hammingWindow;
    else
        currentWindow = hamming(currentWindowLength);
    end
    windowedImpulse(startTime:endTime) = windowedImpulse(startTime:endTime) + currentWindow(:).';
end


% LPC filter
allPoleFilter = filter(1, lpcCoefficients, windowedImpulse);
allPoleFilter = allPoleFilter / max(abs(allPoleFilter));
sound(allPoleFilter, sampleRate);

audiowrite('synthesised_had_f.wav', allPoleFilter, sampleRate);


%Processed speech waveforms and the original waveforms 
t = (0:length(isolatedSegment) -1)/sampleRate;
figure; 
subplot(2,1,1);
plot(t,isolatedSegment);
title("Initial Sample");
xlabel("Time(s)");
ylabel("Amplitude");
subplot(2,1,2);
plot(t, processedSpeech);
title("Processsed Speech");
xlabel("Time(s)");
ylabel("Amplitude");


% Frequency response of the LPC filter
[envResponse, envVector] = freqz(1,lpcCoefficients,1024,sampleRate);
freqDomain = fft(isolatedSegment, length(envResponse));
fftVector = (0:length(freqDomain)-1) * sampleRate / length(freqDomain);


% Spectogram envelope
figure;
plot(envVector, 20*log10(abs(envResponse)));
hold on;
plot(fftVector, 20*log10(abs(freqDomain)));
title("Spectrum Envelope of Isolated Segment 'had' ")
xlabel("Frequency(kHz)");
ylabel("Amplitude(dB)");














