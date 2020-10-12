function [localIntensity, intensityLevel, tDInstantIntensity, InstantIntensityLevel] ...
    = threeDIntensity(mic0,mic1,mic2,mic3,method)

switch method
    %単純な瞬時音響インテンシティの推定
    case 0
        [I01,~,iI01] = InstantIntensity(mic0,mic1);
        [I02,~,iI02] = InstantIntensity(mic0,mic2);
        [I03,~,iI03] = InstantIntensity(mic0,mic3);
        [I12,~,iI12] = InstantIntensity(mic1,mic2);
        [I13,~,iI13] = InstantIntensity(mic1,mic3);
        [I23,~,iI23] = InstantIntensity(mic2,mic3);
        %STFTを使った音響インテンシティの推定
    case 1
        [I01,I02,I03,I12,I13,I23,iI01,iI02,iI03,iI12,iI13,iI23] = stftIntensity(mic0, mic1, mic2, mic3, 44100, 0.05, 707, 1414, 1.1923,512);
end
%右手系ワールド座標へ変換(Unityだと左手系座標になるからy,zが逆)
Ix = -1/4 * (I01 - I02 + I23 - I13 - (2 *I12));
Iz = -1/sqrt(6) * (I01 + I02 + I03);
Iy = -1/(4 * sqrt(3)) * (-I01 - I02 + (2 * I03) + (3 * I13) + (3 * I23));
localIntensity = [Ix, Iy, Iz];
intensityLevel = 10 * log10(norm(localIntensity) / 10^(-12));

% 瞬時音響インテンシティを算出
iIx = -1/4 * (iI01 - iI02 + iI23 - iI13 - (2 *iI12));
iIz = -1/sqrt(6) * (iI01 + iI02 + iI03);
iIy = -1/(4 * sqrt(3)) * (-iI01 - iI02 + (2 * iI03) + (3 * iI13) + (3 * iI23));
% 周期分の足し算
% T = round(44100/c_freq);
% aiIx = zeros(round(length(iI01)/T),1);
% aiIy = zeros(round(length(iI01)/T),1);
% aiIz = zeros(round(length(iI01)/T),1);
% for j = 1:length(iI01)/T
%     aiIx(j) = mean(iIx((j-1)*T+1:j*T));
%     aiIy(j) = mean(iIy((j-1)*T+1:j*T));
%     aiIz(j) = mean(iIz((j-1)*T+1:j*T));
% end
tDInstantIntensity = [iIx,iIy,iIz];
% sampleLength = length(mic0);
InstantIntensityLevel = zeros(length(iI01),1);
for i = 1:length(iI01)
InstantIntensityLevel(i) = 10 * log10(norm(tDInstantIntensity(i,:)) / 10^(-12));
end
end