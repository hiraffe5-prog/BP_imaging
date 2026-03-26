function [ rour3dB, roua3dB, max_Pr, max_Pa, Ia , Ir ] = Reso_PSLR_ISLR_3dB(Ga_in,range_pixel,azimuth_pixel,upr,upa,flag)
%%%%% 求矩阵方位向和距离向的峰值旁瓣比与积分旁瓣比
%%%%% Ga_in:输入矩阵,设第一维为方位向，第二维为距离向
%%%%% range_pixel,azimuth_pixel：输入矩阵像素单元
%%%%% upr,upa：升采样倍数
%%%%% Pa,Ia:方位向PSLR,ISLR
%%%%% Pr,Ir:距离向PSLR,ISLR
%%%%% flag:是否绘图

Ga_in=abs(Ga_in);
[Na_in,Nr_in]=size(Ga_in);
if Na_in == 1 || Nr_in == 1
    error('请输入矩阵')
end
%%%%%%%%%图像升采样%%%%%%%%%%%%%%%
%距离升采样
% upr=4; %%距离向升采样系数
 Nr_up=Nr_in*upr;
 nz=Nr_up-Nr_in;
 ym_f=upr*fft(Ga_in,Nr_in,2);
 sig_upr=zeros(Na_in,Nr_up);
 sig_upr=[ym_f(:,1:Nr_in/2),zeros(Na_in,nz),ym_f(:,(Nr_in/2+1):Nr_in)];
 sig_upr=ifft(sig_upr,[],2);
 
 %方位升采样
% upa=4; %%距离向升采样系数
 Na_up=Na_in*upa;
 nz=Na_up-Na_in;
 ym_f=upa*fft(sig_upr,Na_in,1);
 sig_upa=zeros(Na_up,Nr_up);
 sig_upa=[ym_f(1:Na_in/2,:);zeros(nz,Nr_up);ym_f((Na_in/2+1):Na_in,:)];
 Ga=abs(ifft(sig_upa,[],1));

%升采样
% X=linspace(1,Nr_in,Nr_in);
% Y=linspace(1,Na_in,Na_in);
% 
% XI=linspace(1,Nr_in,upr*Nr_in);
% YI=linspace(1,Na_in,upa*Na_in);
% 
% Ga=interp2(X,Y,Ga_in,XI,YI,'cubic');

range_pixel=range_pixel/upr;
azimuth_pixel=azimuth_pixel/upa;
[Na,Nr]=size(Ga);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[max_val,RI]=max(max(Ga));%%RI是距离向最大值出现的位置
[~,AI]=max(Ga(:,RI));%%AI是方位向最大值出现的位置
max_val_3dB=10^(-3/20)*max_val;

%% 二维图
% x=((0:Nr-1)-Nr/2)*range_pixel;
% y=((0:Na-1)-Na/2)*azimuth_pixel;
% figure,imagesc(Ga)
% xlabel('Range(m)'),ylabel('Azimuth(m)')



Pa=20*log10(abs(Ga(:,RI))./abs(Ga(AI,RI)));
Pr=20*log10(abs(Ga(AI,:))./abs(Ga(AI,RI)));

%% 求解峰值旁瓣比
[r_pks, r_locs] = findpeaks(Pr);
[max_r_pks, max_r_locs] = max(r_pks);
max_Pr = max(r_pks(max_r_locs-1),r_pks(max_r_locs+1));

[a_pks, a_locs] = findpeaks(Pa);
[max_a_pks, max_a_locs] = max(a_pks);
max_Pa = max(a_pks(max_a_locs-1),a_pks(max_a_locs+1));

if(flag)
           
%% 方位向
Pa=20*log10(abs(Ga(:,RI))./abs(Ga(AI,RI)));
figure;
subplot(1,2,1);
plot(abs(Ga(:,RI)));
axis tight;
title('方位向切片'),grid on
subplot(1,2,2);
plot(Pa);
axis([1 Na -50 0]);
title('方位向峰值旁瓣比(db)'),grid on

%% 距离向
Pr=20*log10(abs(Ga(AI,:))./abs(Ga(AI,RI)));
figure;
subplot(1,2,1);
plot(abs(Ga(AI,:)));
axis tight;
title('距离向切片'),grid on
subplot(1,2,2)
plot(Pr);
axis([1 Nr -50 0]);
title('距离向峰值旁瓣比(db)'),grid on

end



%%%%%%%求方位向积分旁瓣比%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i=1:(Na-AI-3)
    %if Ga(AI+i,RI) > Ga(AI+i+1,RI) && Ga(AI+i+1,RI) < Ga(AI+i+2,RI)      %%%直到找到第一个旁瓣,注意防止索引溢出
    if  Ga(AI+i+1,RI) >= max_val_3dB &&  Ga(AI+i+2,RI) < max_val_3dB
        side_down=AI+i+1;
        break;
    end
end

for i=1:(AI-3)
    %if Ga(AI-i,RI)> Ga(AI-i-1,RI)&& Ga(AI-i-1,RI)<Ga(AI-i-2,RI)      %%%直到找到第一个旁瓣
    if Ga(AI-i-1,RI) >= max_val_3dB && Ga(AI-i-2,RI) < max_val_3dB
        side_up=AI-i-1;
        break;
    end
end



%% 主瓣零点对应4db宽度
azi_main_lobe=abs(side_up-side_down);


%%计算主瓣能量

mainEa=sum((Ga(AI-azi_main_lobe:AI+azi_main_lobe,RI)).^2);

%%计算总能量
%Ea=sum(Ga(1:Na,RI).^2);
Ea=sum(Ga(AI-10*azi_main_lobe:AI+10*azi_main_lobe,RI).^2);


%%计算旁瓣能量
sideEa=Ea-mainEa;
Ia=-10*log10(abs(mainEa/sideEa))

%%%%%%%求方位分辨率%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
roua3dB=abs(side_up-side_down)*azimuth_pixel;

%%%%%%%求距离向积分旁瓣比%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i=1:(Nr-RI-3)
    %if Ga(AI,RI+i) > Ga(AI,RI+i+1) & Ga(AI,RI+i+1) < Ga(AI,RI+i+2)      %%%直到找到第一个旁瓣
    if Ga(AI,RI+i+1) >= max_val_3dB &&  Ga(AI,RI+i+2) < max_val_3dB 
        side_right=RI+i+1;
        break;
    end
end

for i=1:(RI-3)
    %if Ga(AI,RI-i)> Ga(AI,RI-i-1)&& Ga(AI,RI-i-1)<Ga(AI,RI-i-2)      %%%直到找到第一个旁瓣
    if Ga(AI,RI-i-1) >= max_val_3dB && Ga(AI,RI-i-2) < max_val_3dB
        side_left=RI-i-1;
        break;
    end
end




%% 主瓣零点对应4db宽度
r_main_lobe=abs(side_left-side_right);%10倍

%%计算主瓣能量
mainEr=sum(Ga(AI,RI-r_main_lobe:RI+r_main_lobe).^2);


%%计算总能量
%Er=sum(Ga(AI,1:Nr).^2);
Er=sum(Ga(AI,RI-10*r_main_lobe:RI+10*r_main_lobe).^2);    %十个旁瓣

%%计算旁瓣能量
sideEr=Er-mainEr;
Ir=-10*log10(abs(mainEr/sideEr))

%%%%%%求距离分辨率%%%%%%%%%%%%%%%%%%%%
rour3dB=abs(side_right-side_left)*range_pixel;


%%%%%%%求二维积分旁瓣比%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%计算主瓣能量
mainE2D=sum(sum(Ga(side_up:side_down,side_left:side_right).^2));

%%计算总能量
E2D=sum(sum(Ga.^2));

%%计算旁瓣能量
sideE2D=E2D-mainE2D;
I2D=-10*log10(abs(mainE2D/sideE2D))

%save P15   Pr Pa rour4dB roua4dB  Ir Ia
