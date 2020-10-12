function [intensity,lv, IIntensity] = InstantIntensity(sig1, sig2)
fs = 44100;
% bw = 1/3;
%直接法の計算
%バンドパスフィルタ
% float fs … サンプリング周波数
% float freq … カットオフ周波数
% float bw   … 帯域幅
% omega = 2 * pi *  c_freq/fs;
% alpha = sin(omega) * sinh(log(2) / 2 * bw * omega / sin(omega));
 
% a0 =  1 + alpha;
% a1 = (-2 * cos(omega))/a0;
% a2 = (1 - alpha)/a0;
% b0 = alpha/a0;
% b1 =  0;
% b2 = -alpha/a0;
% a = [1,a1,a2]';
% b = [b0,b1,b2]';
% BPsig1 = filter(b,a,sig1);
% BPsig2 = filter(b,a,sig2);

% 大気密度とマイクロホン間隔
rho =  1.1923;
r = 0.05;
q = zeros(length(sig1),1);
% 粒子速度
u = zeros(length(sig1),1);

tau = 1/fs;
for t = 1:length(sig1)
u(t) = -1/(rho*r) * sum(sig1(1:t)-sig2(1:t))*tau;
q(t) = sum(sig1(1:t)-sig2(1:t));
end
% % 位相補正
% if ~isempty(Q)
% U = fft(u);
% U = U.*Q;
% % u = ifft(U,'symmetric');
% u = ifft(U);
% end
intensity = sum((sig1+sig2).*u(1:length(sig1))/2)/length(sig1);
% 瞬時音響インテンシティ
IIntensity = zeros(length(sig1),1);
for n = 1:length(sig1)
IIntensity(n) = (sig1(n)+sig2(n))*u(n)/2;
end
i0 = 10^-12;
lv = 10 * log10(abs(intensity)/i0);
end