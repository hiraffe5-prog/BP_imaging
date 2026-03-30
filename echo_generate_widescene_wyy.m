function [fc,Tp,B,lambda,k,Na,Tr,PRF,Nr,fs,r,xm,DEM,lon_lat_range,Tar_xyz,sigma,Rm_all] = ...
    echo_generate_widescene_wyy(pos_xyz, ...
                            x_par_radar, y_par_radar, z_par_radar, ...
                            a_psi, b_psi, c_psi, d_psi, ...
                            a_theta, b_theta, c_theta, d_theta, ...
                            a_phi, b_phi, c_phi, d_phi)
%% 宽场景非停走回波生成
% 说明：
% 1) 不再进行站心坐标系 -> 月固坐标系转换，pos_xyz 直接作为已知输入
% 2) 雷达/目标传播几何在 MCI(J2000) 中处理
% 3) 目标点在 MCMF 中固定，先生成场景中所有散射点的 MCMF 坐标
% 4) 对所有散射点计算非停走等效单程斜距 Rm=(R1+R2)/2
% 5) 再按原 echo_generate_widescene 的方式生成宽场景回波

%% 参数设置（保持原 echo_generate_widescene 参数）
c = 299792458;
R = 1737400;
addpath 'C:\Users\eric\Desktop\wyy_RD定位\读取DEM'
Tcoh = 120;
PRF  = 20;
fc   = 3.1e9;
Tp   = 1e-4;      % 脉冲持续时间
B    = 0.5e6;
k    = B / Tp;
fs   = 2 * B;
Tr   = 1 / PRF;
Na   = round(PRF * Tcoh);
Na   = Na + mod(Na,2);   % 保证 Na 为偶数
Tcoh = Na * Tr;
lambda = c / fc;

%% 设置测绘区域（保持原 echo_generate_widescene 风格）
globalDEM_path = 'DEM_global_moon.mat';
lon_lat_range = [-13,-9,-41,-45];
% lon_lat_range = [-12,-10,-42,-44];
% lon_lat_range = [-15,-7,-39,-47];
del_lat = 0.02;
del_lon = 0.02;

lon = linspace(lon_lat_range(1), lon_lat_range(2), ...
    abs(lon_lat_range(2)-lon_lat_range(1))/del_lon);
lat = linspace(lon_lat_range(3), lon_lat_range(4), ...
    abs(lon_lat_range(4)-lon_lat_range(3))/del_lat);

[DEM_1] = DEM_crop_fun(globalDEM_path, lon_lat_range);
[DEM]   = DEM_resample(DEM_1, lon_lat_range, lon, lat);

figure, imagesc(lon, lat, DEM)
set(gca,'YDir','normal')
xlabel('经度/deg'), ylabel('纬度/deg'), title('高程图')

%% 场景网格点（MCMF）
[row, col] = size(DEM);

del_lat = abs(lon_lat_range(3)-lon_lat_range(4)) / row;
del_lon = abs(lon_lat_range(1)-lon_lat_range(2)) / col;

B_vec = (lon_lat_range(3):-del_lat:lon_lat_range(4)+del_lat) * pi/180;   % 纬度(rad)
L_vec = (lon_lat_range(1): del_lon: lon_lat_range(2)-del_lon) * pi/180;  % 经度(rad)

TarNum = length(B_vec) * length(L_vec);
Tar_xyz = zeros(TarNum, 3);
% 
h1 = waitbar(0, '散射点位置计算中...');
for n = 1:length(B_vec)
    waitbar(n/length(B_vec), h1);
    for m = 1:length(L_vec)
        Tar_xyz((n-1)*length(L_vec)+m, :) = (R + DEM(n,m)) * ...
            [cos(B_vec(n))*cos(L_vec(m)), ...
             cos(B_vec(n))*sin(L_vec(m)), ...
             sin(B_vec(n))];
    end
end
close(h1);

%% 散射系数
% 这里仍沿用原函数中的 scatterpara
[incidence, sigma] = scatterpara(R, pos_xyz, DEM, L_vec, B_vec);
sigma = sigma.';
sigma = sigma(:);

figure, imagesc(L_vec*180/pi, B_vec*180/pi, abs(cos(incidence)))
xlabel('经度/deg'), ylabel('纬度/deg')
set(gca,'YDir','normal')

%% 组织非停走参数
ParaNSG.Na = Na;
ParaNSG.Tr = Tr;
ParaNSG.c0 = c;

ParaNSG.x_par_radar = x_par_radar;
ParaNSG.y_par_radar = y_par_radar;
ParaNSG.z_par_radar = z_par_radar;

ParaNSG.a_psi   = a_psi;
ParaNSG.b_psi   = b_psi;
ParaNSG.c_psi   = c_psi;
ParaNSG.d_psi   = d_psi;

ParaNSG.a_theta = a_theta;
ParaNSG.b_theta = b_theta;
ParaNSG.c_theta = c_theta;
ParaNSG.d_theta = d_theta;

ParaNSG.a_phi   = a_phi;
ParaNSG.b_phi   = b_phi;
ParaNSG.c_phi   = c_phi;
ParaNSG.d_phi   = d_phi;

%% 计算所有散射点的非停走等效单程斜距
disp('非停走等效单程斜距计算中...')

xp2 = reshape(Tar_xyz(:,1), length(B_vec), length(L_vec));
yp2 = reshape(Tar_xyz(:,2), length(B_vec), length(L_vec));
zp2 = reshape(Tar_xyz(:,3), length(B_vec), length(L_vec));

Rm_all = calc_Rm_all_nonstopgo(xp2, yp2, zp2, ParaNSG);   % [row, col, Na]

% 拉平成 [Na, TarNum]，便于后续逐点生成回波
Rm = zeros(Na, TarNum);
for kslow = 1:Na
    tmp = Rm_all(:,:,kslow);
    Rm(kslow,:) = tmp(:).';
end

Rmax = max(Rm(:));
Rmin = min(Rm(:));

%% 产生回波（保留原函数思路，只是把 Rm 换成非停走等效单程斜距）
disp('生成回波信号...')

win1 = 0.1 * Tp;
win2 = Tp;
Nr = round((2*Rmax/c - 2*Rmin/c + win1 + win2) * fs);
Nr = Nr + mod(Nr,2);

tn = linspace(-win1 + 2*Rmin/c, win2 + 2*Rmax/c, Nr);
fs = Nr / (2*Rmax/c - 2*Rmin/c + win1 + win2);
r  = tn * c / 2;

xm = zeros(Na, Nr);

tng = gpuArray(ones(Na,1) * tn);
point_pha = zeros(1, TarNum);

for m = 1:TarNum
    if mod(m, max(1, floor(TarNum/200))) == 0 || m == TarNum
        fprintf('回波生成进度：%d / %d (%.2f%%)\n', m, TarNum, 100*m/TarNum);
    end

    td = gpuArray(2 * Rm(:,m) / c);
    td_tn = tng - td * ones(1, Nr);
    x = rectpuls(td_tn - Tp/2, Tp) .* exp(1i*pi*k*(td_tn - Tp/2).^2);

    xm = xm + sigma(m) * x .* ...
        (exp(-2j*pi*fc*td) .* exp(2j*pi*point_pha(m)) * ones(1, Nr));
end

xm = gather(xm);

disp('已生成回波信号')

end


%% ============================================================
% 计算所有网格点、所有PRT的非停走等效单程斜距
% 输出：
%   Rm_all(i,j,k) 表示第k个PRT时，第(i,j)个网格点的等效单程斜距
%% ============================================================
function Rm_all = calc_Rm_all_nonstopgo(xp2, yp2, zp2, ParaCBP)

[Ny, Nx] = size(xp2);
Na = ParaCBP.Na;
Np = Ny * Nx;

xp_vec = xp2(:);
yp_vec = yp2(:);
zp_vec = zp2(:);

Rm_mat = zeros(Np, Na);

% 进度条
h = waitbar(0, '开始计算 Rm_all ...', 'Name', '非停走斜距计算');
dq = parallel.pool.DataQueue;
numDone = 0;
updateStep = max(1, floor(Np / 200));

afterEach(dq, @updateWaitbar);

    function updateWaitbar(n)
        numDone = numDone + n;
        if isvalid(h)
            waitbar(min(numDone/Np,1), h, ...
                sprintf('Rm计算中: %d / %d (%.2f%%)', ...
                numDone, Np, min(numDone/Np,1)*100));
        end
    end

parfor p = 1:Np
    Tar_xyz = [xp_vec(p); yp_vec(p); zp_vec(p)];
    Rm_mat(p,:) = calc_one_pixel_Rm_nonstopgo(Tar_xyz, ParaCBP).';
    if mod(p, updateStep) == 0
        send(dq, updateStep);
    end
end

remain = mod(Np, updateStep);
if remain ~= 0
    send(dq, remain);
end

if isvalid(h)
    close(h);
end

Rm_all = zeros(Ny, Nx, Na);
for k = 1:Na
    Rm_all(:,:,k) = reshape(Rm_mat(:,k), Ny, Nx);
end

end


%% ============================================================
% 单个散射点的非停走等效单程斜距
%% ============================================================
function Rm_vec = calc_one_pixel_Rm_nonstopgo(Tar_xyz, ParaCBP)

Na = ParaCBP.Na;
Tr = ParaCBP.Tr;
c0 = ParaCBP.c0;

x_par_radar = ParaCBP.x_par_radar;
y_par_radar = ParaCBP.y_par_radar;
z_par_radar = ParaCBP.z_par_radar;

a_psi   = ParaCBP.a_psi;
b_psi   = ParaCBP.b_psi;
c_psi   = ParaCBP.c_psi;
d_psi   = ParaCBP.d_psi;

a_theta = ParaCBP.a_theta;
b_theta = ParaCBP.b_theta;
c_theta = ParaCBP.c_theta;
d_theta = ParaCBP.d_theta;

a_phi   = ParaCBP.a_phi;
b_phi   = ParaCBP.b_phi;
c_phi   = ParaCBP.c_phi;
d_phi   = ParaCBP.d_phi;

t1_vec = (0:Na-1).' * Tr;
Rm_vec = zeros(Na,1);

psi_fun   = @(t) d_psi   * t.^3 + c_psi   * t.^2 + b_psi   * t + a_psi;
theta_fun = @(t) d_theta * t.^3 + c_theta * t.^2 + b_theta * t + a_theta;
phi_fun   = @(t) d_phi   * t.^3 + c_phi   * t.^2 + b_phi   * t + a_phi;

Rz = @(ang) [ cos(ang),  sin(ang), 0;
             -sin(ang),  cos(ang), 0;
                   0   ,      0   , 1 ];

Rx = @(ang) [1,     0     ,      0    ;
             0, cos(ang),  sin(ang);
             0, -sin(ang), cos(ang)];

tran = @(t) Rz(psi_fun(t)) * Rx(theta_fun(t)) * Rz(phi_fun(t));

xr_fun = @(t) x_par_radar(1)*t.^3 + x_par_radar(2)*t.^2 + x_par_radar(3)*t + x_par_radar(4);
yr_fun = @(t) y_par_radar(1)*t.^3 + y_par_radar(2)*t.^2 + y_par_radar(3)*t + y_par_radar(4);
zr_fun = @(t) z_par_radar(1)*t.^3 + z_par_radar(2)*t.^2 + z_par_radar(3)*t + z_par_radar(4);
pr_fun = @(t) [xr_fun(t); yr_fun(t); zr_fun(t)];

pt_fun = @(t) tran(t) * Tar_xyz;

tau12_prev = 1.3;
tau23_prev = 1.3;

for k = 1:Na
    t1 = t1_vec(k);
    pr_t1 = pr_fun(t1);

    % ---------- 求 t2 ----------
    fun_t2 = @(t) norm(pt_fun(t) - pr_t1) - c0 * (t - t1);

    if k == 1
        R1_init = norm(pt_fun(t1) - pr_t1);
        t2_init = t1 + R1_init / c0;
    else
        t2_init = t1 + tau12_prev;
    end

    t2 = local_fzero_with_expand(fun_t2, t2_init, t1, t1 + 5.0);

    pt_t2 = pt_fun(t2);
    R1 = norm(pt_t2 - pr_t1);

    % ---------- 求 t3 ----------
    fun_t3 = @(t) norm(pr_fun(t) - pt_t2) - c0 * (t - t2);

    if k == 1
        t3_init = t2 + R1 / c0;
    else
        t3_init = t2 + tau23_prev;
    end

    t3 = local_fzero_with_expand(fun_t3, t3_init, t2, t2 + 5.0);

    pr_t3 = pr_fun(t3);
    R2 = norm(pr_t3 - pt_t2);

    Rm_vec(k) = (R1 + R2) / 2;

    tau12_prev = t2 - t1;
    tau23_prev = t3 - t2;
end

end


%% ============================================================
% 带区间扩展的 fzero 求根
%% ============================================================
function root = local_fzero_with_expand(fun, x0, left_min, right_max)

try
    root = fzero(fun, x0);
    if isreal(root) && root >= left_min && root <= right_max
        return;
    end
catch
end

step_list = [1e-6, 1e-5, 1e-4, 1e-3, 1e-2, 5e-2, 1e-1, 2e-1, 5e-1, 1.0];

for s = step_list
    a = max(left_min, x0 - s);
    b = min(right_max, x0 + s);

    fa = fun(a);
    fb = fun(b);

    if ~isnan(fa) && ~isnan(fb) && isfinite(fa) && isfinite(fb)
        if fa == 0
            root = a;
            return;
        elseif fb == 0
            root = b;
            return;
        elseif fa * fb < 0
            root = fzero(fun, [a, b]);
            return;
        end
    end
end

Ngrid = 2000;
grid = linspace(left_min, right_max, Ngrid);
vals = nan(size(grid));

for i = 1:Ngrid
    vals(i) = fun(grid(i));
end

idx = find(isfinite(vals(1:end-1)) & isfinite(vals(2:end)) & ...
           vals(1:end-1).*vals(2:end) <= 0, 1);

if ~isempty(idx)
    root = fzero(fun, [grid(idx), grid(idx+1)]);
    return;
end

error('fzero 求根失败：在给定区间内没有找到可靠根。');
end
