%Load MP3 file and check that it plays
% [audioData, sampleRate] = audioread("heed_f.wav");
[audioData, sampleRate] = audioread("had_m.wav");
% sound(audioData, sampleRate);

segmentStart = 0.1; 
segmentLength = 0.2; 

startTime = ceil(segmentStart * sampleRate); 
endTime = floor((segmentStart + segmentLength)* sampleRate);

endTime = min(endTime, length(audioData));
isolatedSegment = audioData(startTime:endTime);


preEmphasis = [1 -0.85];
isolatedSegment = filter(preEmphasis, 1, isolatedSegment);


lpcOrder = 30;
[lpcCoefficients, ~] = lpc(isolatedSegment, lpcOrder);


%Formant estimates using roots of lpc 
rts = roots(lpcCoefficients);
rts = rts(imag(rts)>= 0);
rts = rts(abs(abs(rts)-1) < 0.013);
angz = atan2(imag(rts), real(rts));
formantFrequencies = angz * (sampleRate/(2*pi));
sortedFormantFrequencies = sort(formantFrequencies);
% formantFrequencies = sortedFormantFrequencies(sortedFormantFrequencies > 350 & formantFrequencies < 3000);
formantFrequencies = sortedFormantFrequencies(sortedFormantFrequencies > 100 & formantFrequencies < 900);

format short g;
disp(sortedFormantFrequencies);

%using freqz, calculate the frequency response of the LPC filter
% [envResponse, envVector] = freqz(1,lpcCoefficients,1024,sampleRate);

[autoCorrelationValues, autoCorrelationLags] = xcorr(isolatedSegment, "coeff");

minLag = ceil(sampleRate/ 265);

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


% %Plot the spectogram of the original speech
% figure;
% subplot(2,1,1);
% spectrogram(audioData, 256, 250, 256,sampleRate,"yaxis");
% title("Spectogram of Intial Sample")

% Fundamental frequencies using autocorrelation methods
% [autoCorrelationValues, autoCorrelationLags] = xcorr(isolatedSegment, "coeff");
% minLag = ceil(sampleRate/155); %set to filter out male noise
% minLag = ceil(sampleRate/85); %set to filter out female noise
% maxLag = floor(sampleRate/255);
% positiveLagsIndex = find(autoCorrelationLags(autoCorrelationLags >= minLag & maxLag >= autoCorrelationLags));
% positiveLags = autoCorrelationValues(positiveLagsIndex);
% positiveAutoCorrelationValues = autoCorrelationLags(positiveLagsIndex);
% 
% 
% [maxValue, maxIndex] = max(positiveAutoCorrelationValues);
% periodLag = positiveLags(maxIndex);
% 
% meanFundamentalFrequency = sampleRate/periodLag;
% disp("Mean Fundamental Frequency in Hz: ");
% disp(meanFundamentalFrequency);
% 
% 
% % Frequency Domain
% freqDomain = fft(isolatedSegment, length(envResponse));
% fftVector = (0:length(freqDomain)-1) * sampleRate / length(freqDomain);

% plot the spectogram envelope
% figure;
% plot(envVector, 20*log10(abs(envResponse)));
% hold on;
% plot(fftVector, 20*log10(abs(freqDomain)));
% title("Spectrum Envelope of Isolated Segment")
% xlabel("Frequency(Hz)");
% ylabel("Amplitude(dB)");
% 
% figure;
% plot(autoCorrelationLags, autoCorrelationValues);
% title('Autocorrelation of the Signal');
% xlabel('Lag');
% ylabel('Autocorrelation value');

% periodic impulse signal from fundamental frequency(use lpc coefficients
% on filter)
%get signal to sound like a vowel(aa)

f = meanFundamentalFrequency;
if isnan(f) || f <= 0
    disp("invalid");
    return;
end 

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



% Frequency response of LPC filter
[h, w] = freqz(1, lpcCoefficients, 1024, sampleRate);
plot(w, 20*log10(abs(h)));
title('LPC Spectrum');
xlabel('Frequency (Hz)');
ylabel('Magnitude (dB)');

% % Spectrum of the original signal
% [originalSpectrum, f] = periodogram(isolatedSegment, hamming(length(isolatedSegment)), 1024, sampleRate);
% hold on;
% plot(f, 20*log10(originalSpectrum), 'r');
% legend('LPC Spectrum', 'Original Spectrum');

disp(autoCorrelationValues);
plot(autoCorrelationLags, autoCorrelationValues);
title('Full Autocorrelation');
xlabel('Lags');
ylabel('Autocorrelation Value');


% % Calculate the autocorrelation
% [autoCorrelationValues, autoCorrelationLags] = xcorr(isolatedSegment, 'coeff');

% Plot the autocorrelation
% figure; % Creates a new figure window
% plot(autoCorrelationLags, autoCorrelationValues);
% title('Autocorrelation of the Signal');
% xlabel('Lag');
% ylabel('Autocorrelation value');
% grid on; % Adds a grid for easier reading

% Frequency Domain
freqDomain = fft(isolatedSegment, length(envResponse));
fftVector = (0:length(freqDomain)-1) * sampleRate / length(freqDomain);










