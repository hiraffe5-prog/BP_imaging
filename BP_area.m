% BP对月成像处理
clc,clear,close all
addpath 'C:\Users\eric\Desktop\月球面目标仿真代码\子函数'
%% 参数设置
c = 299792458;
R = 1737400;          % 月球参考半径(m)

Tcoh = 20;
PRF  = 20;
fc   = 3.25e9;
Tp   = 3e-4;          % 脉冲持续时间
B    = 2e6;
k    = B / Tp;
fs   = 4 * B;
Tr   = 1 / PRF;
Na   = round(PRF * Tcoh);
Na   = Na + mod(Na, 2);   % 保证Na为偶数
Tcoh = Na * Tr;
lambda = c / fc;

%% 坐标拟合
EARTH_emp_path = 'C:\Users\eric\Desktop\月球面目标仿真代码\输入文件\WGC_StateVector_20260327093305_radar_J2000.csv';

% 读取雷达在 J2000 月心惯性系（MCI）下的星历
[~, XYZ] = nasa_csv_read(EARTH_emp_path);

% 三维坐标（km -> m）
x_Radar_MCI = XYZ(:,1) * 1000;
y_Radar_MCI = XYZ(:,2) * 1000;
z_Radar_MCI = XYZ(:,3) * 1000;

% 星历采样间隔
dt_ephem = 1;              % s

% 原始星历点数
Ttotal = length(x_Radar_MCI);

% 原始星历时间轴（单位：s）
t_raw = 0:dt_ephem:(Ttotal-1)*dt_ephem;

% 对三个坐标分量分别进行三阶多项式拟合
x_par_radar = polyfit(t_raw, x_Radar_MCI, 3);   % [d_x, c_x, b_x, a_x]
y_par_radar = polyfit(t_raw, y_Radar_MCI, 3);   % [d_y, c_y, b_y, a_y]
z_par_radar = polyfit(t_raw, z_Radar_MCI, 3);   % [d_z, c_z, b_z, a_z]

% 雷达慢时间轴
ta = (0:Na-1)*Tr;

% 在慢时间轴上计算拟合后的雷达轨迹
x_radar_ta = polyval(x_par_radar, ta);
y_radar_ta = polyval(y_par_radar, ta);
z_radar_ta = polyval(z_par_radar, ta);

% 组合成雷达在慢时间轴上的三维位置
p_r_ta = [x_radar_ta(:), y_radar_ta(:), z_radar_ta(:)];

%% 欧拉角拟合
Euler_emp_path = 'C:\Users\eric\Desktop\月球面目标仿真代码\输入文件\WGC_FrameTransformation_20260327104919.csv';

% 读取欧拉角数据
[~, Euler_Angles] = nasa_csv_euler_read(Euler_emp_path);

% 欧拉角：度 -> 弧度
psi_raw   = Euler_Angles(:,1) * pi / 180;   % Angle 3
theta_raw = Euler_Angles(:,2) * pi / 180;   % Angle 2
phi_raw   = Euler_Angles(:,3) * pi / 180;   % Angle 1

% 三阶多项式拟合
psi_par   = polyfit(t_raw, psi_raw.', 3);
theta_par = polyfit(t_raw, theta_raw.', 3);
phi_par   = polyfit(t_raw, phi_raw.', 3);

% 提取系数
d_psi = psi_par(1);     c_psi = psi_par(2);     b_psi = psi_par(3);     a_psi = psi_par(4);
d_theta = theta_par(1); c_theta = theta_par(2); b_theta = theta_par(3); a_theta = theta_par(4);
d_phi = phi_par(1);     c_phi = phi_par(2);     b_phi = phi_par(3);     a_phi = phi_par(4);

%% 数据预处理
%[Tp,B,k,Na,Tr,Nr,fs,r,xm]=pre_process();
[fc,Tp,B,lambda,k,Na,Tr,PRF,Nr,fs,r,xm,Tcoh,Tar_xyz] = echo_generate_wyy(x_par_radar, y_par_radar, z_par_radar, ...
                                                        a_psi, b_psi, c_psi, d_psi, ...
                                                        a_theta, b_theta, c_theta, d_theta, ...
                                                        a_phi, b_phi, c_phi, d_phi);
%[fc,Tp,B,lambda,k,Na,Tr,PRF,Nr,fs,r,S,Tcoh,Tar_xyz]=echo_generate_test();
% [fc,Tp,B,lambda,k,Na,Tr,PRF,Nr,fs,r,xm,DEM,lon_lat_range]=echo_generate_widescene();

figure,imagesc(r,1:Na,abs(xm))
xlabel('距离/m'),ylabel('方位向(采样点)'),title('原始回波')

% %% 读取星历
% %STK
% trace_ref=[0,0]*pi/180; %星历参考点
% str='Target-Target1-To-Target-SJZ AER_2021-10-08T05-30-00Z_2021-10-08T06-30-00Z-1s.csv';
% %str='Target-Target1-To-Target-SJZ AER.csv';
% interal=1;              %星历更新间隔
% pos_xyz=sc2mc(trace_ref,str,interal,Tr,Na);

%NASA
% str='Target-Target1-To-Target-Qujing AER_0912_1s.csv';
% pos_xyz=csvread(str,1,1);

%% 脉压
ref_t=(0:1/fs:(Tp-1/fs)) - Tp/2;
ref=exp(1i*pi*k*ref_t.^2);
%ref=ref.*hamming(length(ref)).';%距离向加窗

ref=conj(fft(ref,Nr,2));
ym=bsxfun(@times,fft(xm,Nr,2),ref);
ym=ifft(ym,Nr,2);

figure;imagesc(r,1:Na,abs(ym));
%colormap(1-gray(64));
xlabel('距离/m'),ylabel('方位向(采样点)'),title('距离脉压')

%% BP处理
%% 升采样
 nup=4; %%距离向升采样系数
 Nr_up=Nr*nup;
 nz=Nr_up-Nr;
 dt=1/nup/fs;
 dr= dt*c/2;
 ym_f=nup*fft(ym,Nr,2);
 sig_up=zeros(Na,Nr_up);
 sig_up=[ym_f(:,1:Nr/2),zeros(Na,nz),ym_f(:,(Nr/2+1):Nr)];
 sig=ifft(sig_up,[],2);
 r1 = r(1) + (0:Nr_up-1)*dr;
%sig=ym;

figure;
imagesc(r1,1:Na,abs(sig));
%colormap(1-gray(64));
xlabel('距离/m'),ylabel('方位向(采样点)'),title('距离插值(4倍)')

%% 步骤二 对成像区域网格化 ，注意是否带DEM,DEM插值等
ParaCBP.Rmin = r(1);
ParaCBP.deltaR = dr;
ParaCBP.Na = Na;
ParaCBP.NrInterp = Nr_up;
ParaCBP.Lambda = lambda;

x_pos = pos_xyz(1,1:ParaCBP.Na);
y_pos = pos_xyz(2,1:ParaCBP.Na);
z_pos = pos_xyz(3,1:ParaCBP.Na);

% lon_range=[-13.5,-8.5];
% lat_range=[-40.5,-45.5];

lon_range=[-11.5,-10.5];
lat_range=[-42.5,-43.5];

lon_lat_range=[lon_range,lat_range];
del_lat = 0.005;
del_lon = 0.005;
lon = linspace(lon_lat_range(1), lon_lat_range(2), abs(lon_lat_range(2)-lon_lat_range(1))/del_lon);
lat = linspace(lon_lat_range(3), lon_lat_range(4), abs(lon_lat_range(4)-lon_lat_range(3))/del_lat);
[lon2D,lat2D] = meshgrid(lon,lat);

globalDEM_path = 'C:\Users\eric\Desktop\wyy_RD定位\读取DEM\DEM_global_moon.mat';
[DEM_1]=DEM_crop_fun(globalDEM_path,lon_lat_range); 
[DEM_interp]=DEM_resample(DEM_1,lon_lat_range,lon,lat);
%DEM_interp = 0;
XX = (R+DEM_interp).*cosd(lat2D).*cosd(lon2D);
YY = (R+DEM_interp).*cosd(lat2D).*sind(lon2D);
ZZ = (R+DEM_interp).*sind(lat2D);
Scene_X = XX;
Scene_Y = YY;
Scene_Z = ZZ;

[amsum,amsum_new,rp_middle] = Func_CBP_2D_GPU( sig(1:ParaCBP.Na,:),ParaCBP,Scene_X,Scene_Y,Scene_Z,x_pos,y_pos,z_pos);

% 显示
figure,imagesc(lon,lat,abs(amsum))
set(gca,'YDir','normal'),xlabel('经度/deg'),ylabel('纬度/deg'),title('BP成像')



