%マイク1がz軸手前にあり、マイク2,3,4がz軸奥にあるからz軸正方向に矢印が向くはず
mymap = jet(256);
liMax = 115;
liMin = 65;

%パラメータ設定
fs = 44100;
nfft = 4096;
f = 1000;%対象周波数
t = 0:1/fs:(nfft-1)/fs;
c = 340;
dr = 0.05;

%マイク1
m1 = dr*[0,0,-sqrt(6)/4]';
phi1 = dr*-sqrt(6)/4/c;
sig1 = sin(2*pi*f*(t-phi1))';

%マイク2,3,4
phi2 = dr*sqrt(6)/12/c;
m2 = dr*[-1/2, -sqrt(3)/6, sqrt(6)/12]';
m3 = dr*[1/2, -sqrt(3)/6, sqrt(6)/12]';
m4 = dr*[0, 1/sqrt(3), sqrt(6)/12]';
sig234 = sin(2*pi*f*(t-phi2))';

%音圧調整
p0 = 2*10^-5;
SPL = 20*log10(rms(sig1)/p0);
%誤差調整用理想パラメータ
rho =  1.1923;
eideal = 10*log10(400/(rho*c));


mic = [m1,m2,m3,m4];
x = mic(1,:);
y = mic(2,:);
z = mic(3,:);


c_freq = 1000;
%時間領域における単純な瞬時音響インテンシティの推定
[dI, dIlv, II, IIlv] = threeDIntensity(sig1,sig234,sig234,sig234,0);
[cI, cIlv] = CrossSpectrumMethod(sig1, sig234,sig234, sig234, fs, dr, 707, 1414, 1.1923,nfft);
[oI, oIlv] = OkawaCrossSpectrum(sig1, sig234,sig234, sig234, fs, dr, 707, 1414, 1.1923,nfft);

%STFTを用いた周波数領域における瞬時音響インテンシティの推定
%% 時間平均プロット
% コーン作成
U = dI(1);
V = dI(2);
W = dI(3);
UVW = [U;V;W] ./ sqrt(U.^2 + V.^2 + W.^2); % 正規化

figure
    
    conescale =  1 * (dIlv - liMin) / (liMax - liMin);
    colorindex = round(256 * (dIlv - liMin) / (liMax - liMin));
    
    if (colorindex > 0) && (colorindex < 257)
        hcone =  coneplot(0,0,0,UVW(1),UVW(2),UVW(3),conescale,'nointerp');
      hcone.EdgeColor = 'none';
      hcone.FaceColor = mymap(colorindex,:);
    else
        hcone =  coneplot(0,0,0,UVW(1),UVW(2),UVW(3),conescale,'nointerp');
      hcone.EdgeColor = 'none';
      hcone.FaceColor = [0 0 0];
    end


% 座標軸の設定
%axis(axislim2)
 axis equal;
xlim([-1 1])
ylim([-1 1])
zlim([-1 1])
xlabel({'x [m]'});
ylabel({'y [m]'});
zlabel({'z [m]'});

set(gca,'fontsize',16);
grid on

% カラーバー作成
colormap(mymap)
caxis([liMin liMax]);
colorbar

% 光源
lgt = camlight('headlight');
hcone.DiffuseStrength = 0.1;
 
% 視点移動
view(70,45);
camlight(lgt,'headlight')

%% 瞬時プロット
% コーン作成
U = II(:,1);
V = II(:,2);
W = II(:,3);
UVWi = [U,V,W] ./ sqrt(U.^2 + V.^2 + W.^2); % 正規化
% figure,
% hold on
% for n = 1:40
% % delete(hcone)
%     conescale =  0.1 * (IIlv(n) - liMin) / (liMax - liMin);
%     colorindex = round(256 * (IIlv(n) - liMin) / (liMax - liMin));
%     
%     if (colorindex > 0) && (colorindex < 257)
%         hcone =  coneplot(0,0,0,UVWi(n,1),UVWi(n,2),UVWi(n,3),conescale,'nointerp');
%       hcone.EdgeColor = 'none';
%       hcone.FaceColor = mymap(colorindex,:);
%     else
%         hcone =  coneplot(0,0,0,UVWi(n,1),UVWi(n,2),UVWi(n,3),conescale,'nointerp');
%       hcone.EdgeColor = 'none';
%       hcone.FaceColor = [0 0 0];
%     end
% 
% % grid on
% 
% hcone.DiffuseStrength = 0.1;
% 
% % 光源
% if n == 1
% lgt = camlight('headlight');
% 
%  
% % 視点移動
% view(70,45);
% camlight(lgt,'headlight')
% end
% end
% % scatter3(0,0,-1)
% scatter3(x,y,z)
% xlabel('x');
% ylabel('y');
% zlabel('z');
%  axis equal;
% xlim([-0.1 0.1])
% ylim([-0.1 0.1])
% zlim([-0.1 0.1])
% xlabel({'x [m]'});
% ylabel({'y [m]'});
% zlabel({'z [m]'});
% titletxt = [' InstantIntensity : ', num2str(n),':',num2str(IIlv(n))];
% 
% title(titletxt);
% set(gca,'fontsize',16);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure,
    % 座標軸の設定
%axis(axislim2)

xlabel({'x [m]'});
ylabel({'y [m]'});
zlabel({'z [m]'});
titletxt = ['DirectMethod'];

title(titletxt);
set(gca,'fontsize',16);
% カラーバー作成
colormap(mymap)
caxis([liMin liMax]);
colorbar
% 光源

lgt = camlight('headlight');

 
% 視点移動
view(70,45);
camlight(lgt,'headlight')

for n = 1:length(IIlv)
% delete(hcone)
    conescale =  0.1 * (IIlv(n) - liMin) / (liMax - liMin);
    colorindex = round(256 * (IIlv(n) - liMin) / (liMax - liMin));
    
    if (colorindex > 0) && (colorindex < 257)
        hcone =  coneplot(0,0,0,UVWi(n,1),UVWi(n,2),UVWi(n,3),conescale,'nointerp');
      hcone.EdgeColor = 'none';
      hcone.FaceColor = mymap(colorindex,:);
    else
        hcone =  coneplot(0,0,0,UVWi(n,1),UVWi(n,2),UVWi(n,3),conescale,'nointerp');
      hcone.EdgeColor = 'none';
      hcone.FaceColor = [0 0 0];
    end

% grid on

hcone.DiffuseStrength = 0.1;

% 光源
if n == 1
lgt = camlight('headlight');

 
% 視点移動
view(70,45);
camlight(lgt,'headlight')
end
pause(0.001)%0.01秒ごとに更新
   n 
% scatter3(0,0,-1)
scatter3(x,y,z)
xlabel('x');
ylabel('y');
zlabel('z');
 axis equal;
xlim([-0.1 0.1])
ylim([-0.1 0.1])
zlim([-0.1 0.1])
xlabel({'x [m]'});
ylabel({'y [m]'});
zlabel({'z [m]'});
titletxt = [' InstantIntensity : ', num2str(n),':',num2str(IIlv(n))];

title(titletxt);
set(gca,'fontsize',16);
end

