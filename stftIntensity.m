function [I01,I02,I03,I12,I13,I23,iI01,iI02,iI03,iI12,iI13,iI23] = stftIntensity(mic0, mic1, mic2, mic3, fs, dr, freqMin, freqMax, atmDensity,nfft)
%% インテンシティ算出

%STFT
[s0,freqArray,t] = spectrogram(mic0,hann(nfft),nfft/2,nfft,fs);
[s1,f,t] = spectrogram(mic1,hann(nfft),nfft/2,nfft,fs);
[s2,f,t] = spectrogram(mic2,hann(nfft),nfft/2,nfft,fs);
[s3,f,t] = spectrogram(mic3,hann(nfft),nfft/2,nfft,fs);
% fftMic0 = fft(mic0,sampleLength) / sampleLength;
% fftMic1 = fft(mic1,sampleLength) / sampleLength;
% fftMic2 = fft(mic2,sampleLength) / sampleLength;
% fftMic3 = fft(mic3,sampleLength) / sampleLength;

%両側クロススペクトル
S01 = conj(s0) .* s1 / (nfft^2);
S02 = conj(s0) .* s2 / (nfft^2);
S03 = conj(s0) .* s3 / (nfft^2);
S12 = conj(s1) .* s2 / (nfft^2);
S13 = conj(s1) .* s3 / (nfft^2);
S23 = conj(s2) .* s3 / (nfft^2);

%% 初期設定　計算

%周波数帯域
df = fs / nfft;
fftIndexMin = ceil(freqMin / df) + 1;
fftIndexMax = ceil(freqMax / df);


%% 片側クロススペクトル
G01 = S01(1:nfft/2+1,:);
G02 = S02(1:nfft/2+1,:);
G03 = S03(1:nfft/2+1,:);
G12 = S12(1:nfft/2+1,:);
G13 = S13(1:nfft/2+1,:);
G23 = S23(1:nfft/2+1,:);
G01(2:end,:) = 4 * G01(2:end,:);
G02(2:end,:) = 4 * G02(2:end,:);
G03(2:end,:) = 4 * G03(2:end,:);
G12(2:end,:) = 4 * G12(2:end,:);
G13(2:end,:) = 4 * G13(2:end,:);
G23(2:end,:) = 4 * G23(2:end,:);

%fftした時の周波数値の行列

%2点間インテンシティ算出
iI01 = zeros(length(t),1);
iI02 = zeros(length(t),1);
iI03 = zeros(length(t),1);
iI12 = zeros(length(t),1);
iI13 = zeros(length(t),1);
iI23 = zeros(length(t),1);
for tt = 1:length(t)
iI01(tt) = sum(imag(G01(fftIndexMin:fftIndexMax,tt)) ./ freqArray(fftIndexMin:fftIndexMax)) / (2 * pi * atmDensity * dr);
iI02(tt) = sum(imag(G02(fftIndexMin:fftIndexMax,tt)) ./ freqArray(fftIndexMin:fftIndexMax)) / (2 * pi * atmDensity * dr);
iI03(tt) = sum(imag(G03(fftIndexMin:fftIndexMax,tt)) ./ freqArray(fftIndexMin:fftIndexMax)) / (2 * pi * atmDensity * dr);
iI12(tt) = sum(imag(G12(fftIndexMin:fftIndexMax,tt)) ./ freqArray(fftIndexMin:fftIndexMax)) / (2 * pi * atmDensity * dr);
iI13(tt) = sum(imag(G13(fftIndexMin:fftIndexMax,tt)) ./ freqArray(fftIndexMin:fftIndexMax)) / (2 * pi * atmDensity * dr);
iI23(tt) = sum(imag(G23(fftIndexMin:fftIndexMax,tt)) ./ freqArray(fftIndexMin:fftIndexMax)) / (2 * pi * atmDensity * dr);
end
I01 = sum(iI01)/length(mic1);
I02 = sum(iI02)/length(mic1);
I03 = sum(iI03)/length(mic1);
I12 = sum(iI12)/length(mic1);
I13 = sum(iI13)/length(mic1);
I23 = sum(iI23)/length(mic1);
end
