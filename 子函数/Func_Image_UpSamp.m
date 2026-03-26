function [ Ga_out ] = Func_Image_UpSamp(Ga_in,upr,upa)
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
 Ga_out = abs(ifft(sig_upa,[],1));
 
 
end