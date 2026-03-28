function [amsum,amsum_new,rp_middle] = Func_CBP_2D_GPU_nsg(s1,ParaCBP,xp2,yp2,zp2,x_pos,y_pos,z_pos)
% 非停走BP
% 说明：
% 1) 这里对每个网格点、每个PRT，先求非停走条件下的等效单程斜距 Rm = (R1+R2)/2
% 2) 再按照原BP代码的思路，在第i个方位向上取 rp = Rm(:,:,i) 进行距离索引与相参积累
% 3) 相位补偿时，取中间方位向时刻对应的 rp_middle = Rm(:,:,round(Na/2)+1)
%
% 需要在 ParaCBP 中额外提供以下字段：
%   ParaCBP.c0
%   ParaCBP.Tr
%   ParaCBP.x_par_radar
%   ParaCBP.y_par_radar
%   ParaCBP.z_par_radar
%   ParaCBP.a_psi, ParaCBP.b_psi, ParaCBP.c_psi, ParaCBP.d_psi
%   ParaCBP.a_theta, ParaCBP.b_theta, ParaCBP.c_theta, ParaCBP.d_theta
%   ParaCBP.a_phi, ParaCBP.b_phi, ParaCBP.c_phi, ParaCBP.d_phi

r_start = ParaCBP.Rmin;
delta_r = ParaCBP.deltaR;
Na      = ParaCBP.Na;
Nr_up   = ParaCBP.NrInterp;
lambda  = ParaCBP.Lambda;

%% ============================================================
% 先在CPU上计算所有网格点、所有PRT对应的非停走等效单程斜距
% Rm_all 大小与 xp2 一致，并在第3维上对应慢时间
%% ============================================================
disp('开始计算所有网格点的非停走等效单程斜距 Rm ...');

Rm_all = calc_Rm_all_nonstopgo(xp2, yp2, zp2, ParaCBP);   % [Ny, Nx, Na]


%% 转到GPU做BP积累
s1     = gpuArray(s1);
Rm_all = gpuArray(Rm_all);

h = waitbar(0,'1','Name','bp处理');

amsum = gpuArray.zeros(size(xp2));   % 存储幅度

for i = 1:Na
    s_temp = s1(i,:);

    % 当前方位向（当前PRT）下，所有网格点的非停走等效单程斜距
    rp = Rm_all(:,:,i);

    n_round = round((rp - r_start) / delta_r + 1);
    n_round(n_round < 1)    = 1;
    n_round(n_round > Nr_up)= Nr_up;

    amsum = amsum + s_temp(n_round) .* exp(1j * 4 * pi * rp / lambda);

    waitbar(i/Na, h, sprintf('%6.2f %%', i/Na*100));
end
delete(h)

%% 相位补偿（保相）
midIdx = round(Na/2) + 1;
rp_middle = Rm_all(:,:,midIdx);
amsum_new = amsum .* exp(-1i * 4 * pi * rp_middle / lambda);

amsum     = gather(amsum);
amsum_new = gather(amsum_new);
rp_middle = gather(rp_middle);

end

function Rm_all = calc_Rm_all_nonstopgo(xp2, yp2, zp2, ParaCBP)

[Ny, Nx] = size(xp2);
Na = ParaCBP.Na;
Np = Ny * Nx;

% 拉平成列向量，便于 parfor
xp_vec = xp2(:);
yp_vec = yp2(:);
zp_vec = zp2(:);

% 每个像素对应一条 Na×1 的 Rm 曲线
% 为了适配 parfor，这里先存成 Np×Na
Rm_mat = zeros(Np, Na);

% ========= 新增：parfor进度条 =========
h = waitbar(0, '开始计算 Rm_all ...', 'Name', '非停走斜距计算');
dq = parallel.pool.DataQueue;
numDone = 0;
updateStep = max(1, floor(Np / 200));   % 大约更新200次

afterEach(dq, @updateWaitbar);

    function updateWaitbar(n)
        numDone = numDone + n;
        waitbar(min(numDone / Np, 1), h, ...
            sprintf('Rm计算中: %d / %d (%.2f%%)', ...
            numDone, Np, min(numDone / Np, 1) * 100));
    end
% ====================================

parfor p = 1:Np
    Tar_xyz = [xp_vec(p); yp_vec(p); zp_vec(p)];
    Rm_mat(p, :) = calc_one_pixel_Rm_nonstopgo(Tar_xyz, ParaCBP).';
    if mod(p, updateStep) == 0
        send(dq, updateStep);
    end
end
% 补上最后不足 updateStep 的部分
remain = mod(Np, updateStep);
if remain ~= 0
    send(dq, remain);
end
% 关闭进度条
if isvalid(h)
    close(h);
end
% 再恢复成 Ny×Nx×Na
Rm_all = zeros(Ny, Nx, Na);
for k = 1:Na
    Rm_all(:,:,k) = reshape(Rm_mat(:,k), Ny, Nx);
end

end

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

% ---------- 欧拉角函数 ----------
psi_fun   = @(t) d_psi   * t.^3 + c_psi   * t.^2 + b_psi   * t + a_psi;
theta_fun = @(t) d_theta * t.^3 + c_theta * t.^2 + b_theta * t + a_theta;
phi_fun   = @(t) d_phi   * t.^3 + c_phi   * t.^2 + b_phi   * t + a_phi;

% ---------- 基本旋转矩阵 ----------
Rz = @(ang) [ cos(ang),  sin(ang), 0;
             -sin(ang),  cos(ang), 0;
                   0   ,      0   , 1 ];

Rx = @(ang) [1,     0     ,      0    ;
             0, cos(ang),  sin(ang);
             0, -sin(ang), cos(ang)];

% ---------- MCMF -> MCI ----------
tran = @(t) Rz(psi_fun(t)) * Rx(theta_fun(t)) * Rz(phi_fun(t));

% ---------- 雷达 MCI 坐标函数 ----------
xr_fun = @(t) x_par_radar(1)*t.^3 + x_par_radar(2)*t.^2 + x_par_radar(3)*t + x_par_radar(4);
yr_fun = @(t) y_par_radar(1)*t.^3 + y_par_radar(2)*t.^2 + y_par_radar(3)*t + y_par_radar(4);
zr_fun = @(t) z_par_radar(1)*t.^3 + z_par_radar(2)*t.^2 + z_par_radar(3)*t + z_par_radar(4);
pr_fun = @(t) [xr_fun(t); yr_fun(t); zr_fun(t)];

% ---------- 目标在 MCI 下的位置函数 ----------
pt_fun = @(t) tran(t) * Tar_xyz;

tau12_prev = 1.3;
tau23_prev = 1.3;

for k = 1:Na
    t1 = t1_vec(k);
    pr_t1 = pr_fun(t1);

    %% ---------- 求 t2 ----------
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

    %% ---------- 求 t3 ----------
    fun_t3 = @(t) norm(pr_fun(t) - pt_t2) - c0 * (t - t2);

    if k == 1
        t3_init = t2 + R1 / c0;
    else
        t3_init = t2 + tau23_prev;
    end

    t3 = local_fzero_with_expand(fun_t3, t3_init, t2, t2 + 5.0);

    pr_t3 = pr_fun(t3);
    R2 = norm(pr_t3 - pt_t2);

    %% ---------- 等效单程斜距 ----------
    Rm_vec(k) = (R1 + R2) / 2;

    tau12_prev = t2 - t1;
    tau23_prev = t3 - t2;
end

end

%% ============================================================
% 子函数：带区间扩展的 fzero 求根
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

idx = find(isfinite(vals(1:end-1)) & isfinite(vals(2:end)) & vals(1:end-1).*vals(2:end) <= 0, 1);

if ~isempty(idx)
    root = fzero(fun, [grid(idx), grid(idx+1)]);
    return;
end

error('fzero 求根失败：在给定区间内没有找到可靠根。');
end
