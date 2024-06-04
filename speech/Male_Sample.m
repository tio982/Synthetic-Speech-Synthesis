% Load MP3 file and check that it plays 
[audioData, sampleRate] = audioread("had_m.wav");
sound(audioData, sampleRate)