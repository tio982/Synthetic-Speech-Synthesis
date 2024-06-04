% % Load MP3 file and check that it plays 
% [audioData, sampleRate] = audioread("ee.wav");
% 
% disp(size(audioData))
% sound(audioData, sampleRate)
% 
% disp(sampleRate)

player = audioplayer(audioData, sampleRate);
play(player);
