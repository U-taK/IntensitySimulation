function [localIntensity, intensityLevel] = CrossSpectrumMethod(mic0, mic1, mic2, mic3, fs, dr, freqMin, freqMax, atmDensity,arrayLeng)
%% �����ݒ�@�v�Z

%�T���v����
sampleLength = arrayLeng;

%���g���ш�
df = fs / sampleLength;
fftIndexMin = ceil(freqMin / df) + 1;
fftIndexMax = ceil(freqMax / df);


%% �C���e���V�e�B�Z�o

%FFT
fftMic0 = fft(mic0,sampleLength) / sampleLength;
fftMic1 = fft(mic1,sampleLength) / sampleLength;
fftMic2 = fft(mic2,sampleLength) / sampleLength;
fftMic3 = fft(mic3,sampleLength) / sampleLength;

%�����N���X�X�y�N�g��
S01 = conj(fftMic0) .* fftMic1;
S02 = conj(fftMic0) .* fftMic2;
S03 = conj(fftMic0) .* fftMic3;
S12 = conj(fftMic1) .* fftMic2;
S13 = conj(fftMic1) .* fftMic3;
S23 = conj(fftMic2) .* fftMic3;

%�Б��N���X�X�y�N�g��
G01 = S01(1:sampleLength/2+1);
G02 = S02(1:sampleLength/2+1);
G03 = S03(1:sampleLength/2+1);
G12 = S12(1:sampleLength/2+1);
G13 = S13(1:sampleLength/2+1);
G23 = S23(1:sampleLength/2+1);
G01(2:end-1) = 2 * G01(2:end-1);
G02(2:end-1) = 2 * G02(2:end-1);
G03(2:end-1) = 2 * G03(2:end-1);
G12(2:end-1) = 2 * G12(2:end-1);
G13(2:end-1) = 2 * G13(2:end-1);
G23(2:end-1) = 2 * G23(2:end-1);

%fft�������̎��g���l�̍s��
freqArray = (0:sampleLength/2)' * df;

%2�_�ԃC���e���V�e�B�Z�o
I01 = sum(imag(G01(fftIndexMin:fftIndexMax)) ./ freqArray(fftIndexMin:fftIndexMax),1) / (2 * pi * atmDensity * dr);
I02 = sum(imag(G02(fftIndexMin:fftIndexMax)) ./ freqArray(fftIndexMin:fftIndexMax),1) / (2 * pi * atmDensity * dr);
I03 = sum(imag(G03(fftIndexMin:fftIndexMax)) ./ freqArray(fftIndexMin:fftIndexMax),1) / (2 * pi * atmDensity * dr);
I12 = sum(imag(G12(fftIndexMin:fftIndexMax)) ./ freqArray(fftIndexMin:fftIndexMax),1) / (2 * pi * atmDensity * dr);
I13 = sum(imag(G13(fftIndexMin:fftIndexMax)) ./ freqArray(fftIndexMin:fftIndexMax),1) / (2 * pi * atmDensity * dr);
I23 = sum(imag(G23(fftIndexMin:fftIndexMax)) ./ freqArray(fftIndexMin:fftIndexMax),1) / (2 * pi * atmDensity * dr);

%�E��n���[���h���W�֕ϊ�(Unity���ƍ���n���W�ɂȂ邩��y,z���t)
Ix = -1/4 * (I01 - I02 + I23 - I13 - (2 *I12));
Iy = -1/(4 * sqrt(3)) * (-I01 - I02 + (2 * I03) + (3 * I13) + (3 * I23));
Iz = -1/sqrt(6) * (I01 + I02 + I03);
% Ix = I01 - I02 + I12 - I13 + I23;
% Iz = 1/(2*sqrt(3))*(-I01-I02+2*I03+I13+I23);
% Iy = 1/sqrt(6)*(I01 + I02 + I13);

localIntensity = [Ix, Iy, Iz];
intensityLevel = 10 * log10(norm(localIntensity) / 10^(-12));
end
