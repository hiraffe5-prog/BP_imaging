function [fc,Tp,B,lambda,k,Na,Tr,PRF,Nr,fs,r,xm,Tcoh,Tar_xyz] = echo_generate_wyy(x_par_radar, y_par_radar, z_par_radar, ...
                                                                a_psi, b_psi, c_psi, d_psi, ...
                                                                a_theta, b_theta, c_theta, d_theta, ...
                                                                a_phi, b_phi, c_phi, d_phi)
%% 单点目标回波生成（第谷坑中心）
% 目标中心经纬度：[-11°, -43°]
% 目标高程：-4470 m
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

%% 设置单点目标（第谷坑中心）
TarNum = 1;

% 目标经纬度（弧度）
% Tar = [lon, lat]
Tar = [-11, -43] * pi / 180;

% 目标高程（相对月球参考球）
h_tar = -4470;   % m

% 经纬高 -> MCMF（月固坐标系）
lon_tar = Tar(1);
lat_tar = Tar(2);
r_tar   = R + h_tar;

Tar_xyz = [ ...
    r_tar * cos(lat_tar) * cos(lon_tar), ...
    r_tar * cos(lat_tar) * sin(lon_tar), ...
    r_tar * sin(lat_tar) ...
    ];

%% 非停走误差,场景空变性小于波长，可直接在斜距上补偿
% 计算双程斜距
[R_two_way, Rm, t1_vec, t2_vec, t3_vec, R1_vec, R2_vec, R_stopgo] = ...
    calc_two_way_range_nonstop(Tar_xyz, Na, Tr, ...
                               x_par_radar, y_par_radar, z_par_radar, ...
                               a_psi, b_psi, c_psi, d_psi, ...
                               a_theta, b_theta, c_theta, d_theta, ...
                               a_phi, b_phi, c_phi, d_phi, c);

Rmax = max(Rm);
Rmin = min(Rm);

%% 产生回波
disp('生成回波信号')

win1 = 0.1*Tp;
win2 = Tp;

Nr = round((2 * Rmax / c - 2 * Rmin / c + win1 + win2) * fs);
% Nr = 2^nextpow2(Nr);
Nr = Nr + mod(Nr,2);

tn = linspace(-win1 + 2 * Rmin / c, win2 + 2 * Rmax / c, Nr);
r  = tn * c / 2;
fs=Nr/(2*Rmax/c-2*Rmin/c+win1+win2);

xm = zeros(Na, Nr);

% 单点目标双程时延
td = 2 * Rm / c;    % Na x 1

% 构造“快时间 - 时延”矩阵
td_tn = ones(Na, 1) * tn - td * ones(1, Nr);

% 基带LFM回波包络
x = rectpuls(td_tn - Tp / 2, Tp) .* exp(1i * pi * k * (td_tn - Tp / 2).^2);

% 加入载频传播相位
xm = x .* (exp(-2j * pi * fc * td) * ones(1, Nr));
disp('已生成回波信号')
%% 投影（如后续需要可打开）
% R_ref = Rm(Na/2);
% RD_proj(Tar_xyz, pos_xyz, R_ref);

end
