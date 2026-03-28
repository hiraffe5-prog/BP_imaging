function [R_two_way, Rm_equiv, t1_vec, t2_vec, t3_vec, R1_vec, R2_vec, R_stopgo] = ...
    calc_two_way_range_nonstop(Tar_xyz, Na, Tr, ...
                               x_par_radar, y_par_radar, z_par_radar, ...
                               a_psi, b_psi, c_psi, d_psi, ...
                               a_theta, b_theta, c_theta, d_theta, ...
                               a_phi, b_phi, c_phi, d_phi, c0)
% ============================================================
% 功能：
%   计算单个目标点 Tar_xyz 在整个观测历程中的非停走双程斜距
%
% 输入：
%   Tar_xyz      : 目标点在 MCMF(月固系)下的坐标 [1x3] 或 [3x1]
%   Na           : 脉冲数
%   Tr           : PRT
%   x_par_radar  : 雷达 MCI 坐标 x(t) 的三次拟合系数 [1x4]
%   y_par_radar  : 雷达 MCI 坐标 y(t) 的三次拟合系数 [1x4]
%   z_par_radar  : 雷达 MCI 坐标 z(t) 的三次拟合系数 [1x4]
%   a_psi...d_psi, a_theta...d_theta, a_phi...d_phi
%                : 欧拉角三次拟合系数
%   c0           : 光速，默认 299792458
%
% 输出：
%   R_two_way    : 每个发射时刻对应的双程斜距 [Na x 1]
%   Rm_equiv     : 若仍沿用停走模型成像，可取 R_two_way/2 作为等效单程斜距 [Na x 1]
%   t1_vec       : 发射时刻 [Na x 1]
%   t2_vec       : 到达目标时刻 [Na x 1]
%   t3_vec       : 接收回波时刻 [Na x 1]
%   R1_vec       : 去程斜距 [Na x 1]
%   R2_vec       : 回程斜距 [Na x 1]
%
% 说明：
%   1) Tar_xyz 是 MCMF 固定坐标；
%   2) 目标在 MCI 中的位置由 pt_fun(t)=tran(t)*Tar_xyz 给出；
%   3) 雷达在 MCI 中的位置由三次多项式拟合系数给出。
% ============================================================

if nargin < 19 || isempty(c0)
    c0 = 299792458;
end

Tar_xyz = Tar_xyz(:);   % 强制转成 3x1 列向量

%% ---------------- 慢时间 ----------------
t1_vec = (0:Na-1).' * Tr;

%% ---------------- 欧拉角函数 ----------------
psi_fun   = @(t) d_psi   * t.^3 + c_psi   * t.^2 + b_psi   * t + a_psi;
theta_fun = @(t) d_theta * t.^3 + c_theta * t.^2 + b_theta * t + a_theta;
phi_fun   = @(t) d_phi   * t.^3 + c_phi   * t.^2 + b_phi   * t + a_phi;

%% ---------------- 基本旋转矩阵 ----------------
Rz = @(ang) [ cos(ang),  sin(ang), 0;
             -sin(ang),  cos(ang), 0;
                   0   ,      0   , 1 ];

Rx = @(ang) [1,     0     ,      0    ;
             0, cos(ang),  sin(ang);
             0, -sin(ang), cos(ang)];

%% ---------------- MCMF -> MCI 旋转矩阵 ----------------
tran = @(t) Rz(psi_fun(t)) * Rx(theta_fun(t)) * Rz(phi_fun(t));

%% ---------------- 目标在 t 时刻的 MCI 坐标函数 ----------------
pt_fun = @(t) tran(t) * Tar_xyz;

%% ---------------- 雷达在 t 时刻的 MCI 坐标函数 ----------------
xr_fun = @(t) x_par_radar(1)*t.^3 + x_par_radar(2)*t.^2 + x_par_radar(3)*t + x_par_radar(4);
yr_fun = @(t) y_par_radar(1)*t.^3 + y_par_radar(2)*t.^2 + y_par_radar(3)*t + y_par_radar(4);
zr_fun = @(t) z_par_radar(1)*t.^3 + z_par_radar(2)*t.^2 + z_par_radar(3)*t + z_par_radar(4);

pr_fun = @(t) [xr_fun(t); yr_fun(t); zr_fun(t)];

%% ---------------- 结果预分配 ----------------
t2_vec = zeros(Na,1);
t3_vec = zeros(Na,1);
R1_vec = zeros(Na,1);
R2_vec = zeros(Na,1);
R_two_way = zeros(Na,1);
pr_t1_all = zeros(Na,3);   % 每个t1时刻雷达MCI坐标
pt_t1_all = zeros(Na,3);   % 每个t1时刻目标MCI坐标
R_stopgo  = zeros(Na,1);   % 停走假设单程斜距

%% ---------------- 初始传播时延猜测 ----------------
% 地月单程一般 ~1.2~1.4 s，这里给一个保守的初值与搜索范围
tau12_prev = 1.3;   % 去程初值（秒）
tau23_prev = 1.3;   % 回程初值（秒）

%% ================= 主循环：遍历每个发射时刻 =================
for k = 1:Na
    t1 = t1_vec(k);
    pr_t1 = pr_fun(t1);     % 雷达在t1时刻的MCI坐标
    pt_t1 = pt_fun(t1);     % 目标在t1时刻的MCI坐标

    % 保存
    pr_t1_all(k,:) = pr_t1.';
    pt_t1_all(k,:) = pt_t1.';
    
    % 停走假设单程斜距
    R_stopgo(k) = norm(pr_t1 - pt_t1);

    %% ---------- 第一步：求 t2 ----------
    % 方程：|| pt(t2) - pr(t1) || = c0 * (t2 - t1)
    fun_t2 = @(t) norm(pt_fun(t) - pr_t1) - c0 * (t - t1);

    % 初值：优先用上一脉冲传播时延；第一脉冲则用几何近似
    if k == 1
        R1_init = norm(pt_fun(t1) - pr_t1);
        t2_init = t1 + R1_init / c0;
    else
        t2_init = t1 + tau12_prev;
    end
    % R1_init = norm(pt_fun(t1) - pr_t1);
    % t2_init = t1 + R1_init / c0;

    % 用带搜索区间扩展的 fzero，避免单点初值不稳
    t2 = local_fzero_with_expand(fun_t2, t2_init, t1, t1 + 5.0);

    pt_t2 = pt_fun(t2);
    R1 = norm(pt_t2 - pr_t1);

    %% ---------- 第二步：求 t3 ----------
    % 方程：|| pr(t3) - pt(t2) || = c0 * (t3 - t2)
    fun_t3 = @(t) norm(pr_fun(t) - pt_t2) - c0 * (t - t2);

    if k == 1
        t3_init = t2 + R1 / c0;
    else
        t3_init = t2 + tau23_prev;
    end

    % % ---------- t3 初值 ----------
    % t3_init = t2 + R1 / c0;
    t3 = local_fzero_with_expand(fun_t3, t3_init, t2, t2 + 5.0);

    pr_t3 = pr_fun(t3);
    R2 = norm(pr_t3 - pt_t2);

    %% ---------- 保存 ----------
    t2_vec(k) = t2;
    t3_vec(k) = t3;
    R1_vec(k) = R1;
    R2_vec(k) = R2;
    R_two_way(k) = R1 + R2;

    %% ---------- 更新下一脉冲初值 ----------
    tau12_prev = t2 - t1;
    tau23_prev = t3 - t2;
end

%% 若仍想沿用停走模型的成像框架，可取等效单程斜距
Rm_equiv = R_two_way / 2;

end


%% ============================================================
% 子函数：带区间扩展的 fzero 求根
% 目的：比直接 fzero(fun, x0) 更稳一点
%% ============================================================
function root = local_fzero_with_expand(fun, x0, left_min, right_max)

% 先尝试单点初值
try
    root = fzero(fun, x0);
    if isreal(root) && root >= left_min && root <= right_max
        return;
    end
catch
end

% 再尝试在 x0 附近找变号区间
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

% 如果还不行，就在整个允许区间内粗搜索变号
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
