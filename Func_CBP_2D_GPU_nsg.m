function [amsum, amsum_new, rp_middle] = Func_CBP_2D_GPU_nsg_2(sig, ParaCBP, Scene_X, Scene_Y, Scene_Z)
% 非“停-走”二维CBP（GPU积累版）
% 改写要点：
% 1) 不再一次性计算并存储 Rm_all(:,:,1:Na)
% 2) 改为每个慢时间 i 单独计算当前 rp(:,:,i)，然后立即参与BP积累
% 3) 单像素函数只求第 i 个慢时间对应的等效单程斜距
%
% 输入：
%   sig      : 距离压缩后的回波数据，大小 [Na, Nr_up]
%   ParaCBP  : 参数结构体
%   Scene_X  : 成像网格X坐标 [Ny, Nx]
%   Scene_Y  : 成像网格Y坐标 [Ny, Nx]
%   Scene_Z  : 成像网格Z坐标 [Ny, Nx]
%
% 输出：
%   amsum    : 复图像（聚焦结果）
%   amsum_new: 相位补偿后的图像
%   rp_middle: 中间慢时间对应的等效单程斜距 [Ny, Nx]
%
% 注意：
% 1) 这里默认 ParaCBP 中已有下列字段：
%    Na, Tr, c0, lambda, r_start, delta_r
%    x_par_radar, y_par_radar, z_par_radar
%    a_psi,b_psi,c_psi,d_psi
%    a_theta,b_theta,c_theta,d_theta
%    a_phi,b_phi,c_phi,d_phi
%
% 2) 如果你的原始 ParaCBP 字段名与这里不一致，需要把下面对应位置改成你的字段名。
%
% 3) 如果你的原始 amsum_new 不是 abs(amsum)，只需要改最后几行输出部分即可。

    %% ---------------- 基本参数 ----------------
    r_start = ParaCBP.Rmin;
    delta_r = ParaCBP.deltaR;
    Na      = ParaCBP.Na;
    Nr_up   = ParaCBP.NrInterp;
    lambda  = ParaCBP.Lambda;
    
    xp2 = Scene_X;
    yp2 = Scene_Y;
    zp2 = Scene_Z;

    %% ---------------- 转GPU ----------------
    s1 = gpuArray(sig);

    % 复数累加器
    amsum_gpu = gpuArray.zeros(size(xp2));   % 存储幅度

    %% ---------------- 进度条 ----------------
    h = waitbar(0, '初始化中...', 'Name', '非停走CBP处理中');
    t_start = tic;

    cleanupObj = onCleanup(@() safeCloseWaitbar(h));

    %% ---------------- 主循环：逐慢时间计算并积累 ----------------
    for i = 1:Na
        % 当前慢时间对应的一条距离向数据
        s_temp = s1(i, :);   % [1, Nr_up]

        % 仅计算当前慢时间 i 对应的所有像素等效单程斜距 rp(:,:,i)
        rp = calc_Rm_one_slowtime_nonstopgo(xp2, yp2, zp2, ParaCBP, i);   % [Ny, Nx]

        % 送GPU
        rp_gpu = gpuArray(rp);

        % 斜距 -> 距离向采样点索引
        n_round = round((rp_gpu - r_start) ./ delta_r + 1);

        % 越界裁剪
        n_round(n_round < 1)     = 1;
        n_round(n_round > Nr_up) = Nr_up;

        % BP积累
        amsum_gpu = amsum_gpu + s_temp(n_round) .* exp(1j * 4 * pi * rp_gpu / lambda);

        % 更新进度条
        elapsed_sec = toc(t_start);
        ratio = i / Na;
        if ratio > 0
            remain_sec = elapsed_sec * (1 - ratio) / ratio;
        else
            remain_sec = inf;
        end

        msg = sprintf(['慢时间进度: %d / %d (%.2f%%)\n' ...
                       '已耗时: %s\n' ...
                       '预计剩余: %s'], ...
                       i, Na, ratio*100, ...
                       sec2str(elapsed_sec), sec2str(remain_sec));

        if isgraphics(h)
            waitbar(ratio, h, msg);
        end
        drawnow limitrate;
    end
    %% 相位补偿（保相）

    i_mid = round(Na/2);
    rp_middle = calc_Rm_one_slowtime_nonstopgo(xp2, yp2, zp2, ParaCBP, i_mid);
    amsum_new = amsum_gpu .* exp(-1i * 4 * pi * rp_middle / lambda);

    amsum = gather(amsum_gpu);   % 复图像
    amsum_new = gather(amsum_new);
    rp_middle = gather(rp_middle);

end


%% ========================================================================
function rp = calc_Rm_one_slowtime_nonstopgo(xp2, yp2, zp2, ParaCBP, i_slow)
% 计算“当前第 i_slow 个慢时间”对应的所有网格点等效单程斜距
% 输出：
%   rp [Ny, Nx]

    [Ny, Nx] = size(xp2);
    Np = Ny * Nx;

    xp_vec = xp2(:);
    yp_vec = yp2(:);
    zp_vec = zp2(:);

    rp_vec = zeros(Np, 1, 'like', xp_vec);

    % parfor 按像素并行
    parfor p = 1:Np
        Tar_xyz = [xp_vec(p); yp_vec(p); zp_vec(p)];
        rp_vec(p) = calc_one_pixel_Rm_one_slowtime(Tar_xyz, ParaCBP, i_slow);
    end

    rp = reshape(rp_vec, Ny, Nx);
end


%% ========================================================================
function Rm_val = calc_one_pixel_Rm_one_slowtime(Tar_xyz, ParaCBP, i_slow)
% 单像素、单慢时间：只求当前 i_slow 对应的等效单程斜距
%
% 三个关键时刻：
%   t1 : 发射时刻（当前慢时间）
%   t2 : 到达目标时刻
%   t3 : 返回雷达时刻
%
% 定义：
%   R1 = |pt(t2) - pr(t1)|
%   R2 = |pr(t3) - pt(t2)|
%   Rm = (R1 + R2) / 2

    %% ----- 读参数 -----
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

    %% ----- 当前慢时间发射时刻 -----
    t1 = (i_slow - 1) * Tr;

    %% ----- 欧拉角函数 -----
    psi_fun   = @(t) d_psi   * t.^3 + c_psi   * t.^2 + b_psi   * t + a_psi;
    theta_fun = @(t) d_theta * t.^3 + c_theta * t.^2 + b_theta * t + a_theta;
    phi_fun   = @(t) d_phi   * t.^3 + c_phi   * t.^2 + b_phi   * t + a_phi;

    %% ----- 坐标旋转 -----
    Rz = @(ang) [ cos(ang),  sin(ang), 0;
                 -sin(ang),  cos(ang), 0;
                       0   ,      0   , 1 ];

    Rx = @(ang) [1,     0     ,      0    ;
                 0, cos(ang),  sin(ang);
                 0, -sin(ang), cos(ang)];

    tran = @(t) Rz(psi_fun(t)) * Rx(theta_fun(t)) * Rz(phi_fun(t));

    %% ----- 雷达轨迹函数 -----
    xr_fun = @(t) x_par_radar(1)*t.^3 + x_par_radar(2)*t.^2 + x_par_radar(3)*t + x_par_radar(4);
    yr_fun = @(t) y_par_radar(1)*t.^3 + y_par_radar(2)*t.^2 + y_par_radar(3)*t + y_par_radar(4);
    zr_fun = @(t) z_par_radar(1)*t.^3 + z_par_radar(2)*t.^2 + z_par_radar(3)*t + z_par_radar(4);

    pr_fun = @(t) [xr_fun(t); yr_fun(t); zr_fun(t)];
    pt_fun = @(t) tran(t) * Tar_xyz;

    %% ----- 已知发射时刻雷达位置 -----
    pr_t1 = pr_fun(t1);

    %% ===================== 求 t2 =====================
    % 去程传播约束：
    % |pt(t2) - pr(t1)| = c0 * (t2 - t1)
    fun_t2 = @(t) norm(pt_fun(t) - pr_t1) - c0 * (t - t1);

    % 初值：用 t1 时刻雷达到目标的近似直线距离
    R1_init = norm(pt_fun(t1) - pr_t1);
    t2_init = t1 + R1_init / c0;

    t2 = local_fzero_with_expand(fun_t2, t2_init, t1, t1 + 5.0);

    pt_t2 = pt_fun(t2);
    R1 = norm(pt_t2 - pr_t1);

    %% ===================== 求 t3 =====================
    % 回程传播约束：
    % |pr(t3) - pt(t2)| = c0 * (t3 - t2)
    fun_t3 = @(t) norm(pr_fun(t) - pt_t2) - c0 * (t - t2);

    t3_init = t2 + R1 / c0;

    t3 = local_fzero_with_expand(fun_t3, t3_init, t2, t2 + 5.0);

    pr_t3 = pr_fun(t3);
    R2 = norm(pr_t3 - pt_t2);

    % ----- 等效单程斜距 -----
    Rm_val = (R1 + R2) / 2;

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


%% ========================================================================
function str = sec2str(sec_in)
% 秒数转字符串
    if ~isfinite(sec_in)
        str = '未知';
        return;
    end

    sec_in = max(sec_in, 0);

    hh = floor(sec_in / 3600);
    mm = floor((sec_in - hh*3600) / 60);
    ss = sec_in - hh*3600 - mm*60;

    if hh > 0
        str = sprintf('%02d:%02d:%05.2f', hh, mm, ss);
    else
        str = sprintf('%02d:%05.2f', mm, ss);
    end
end


%% ========================================================================
function safeCloseWaitbar(h)
% 安全关闭进度条
    if ~isempty(h) && isgraphics(h)
        close(h);
    end
end
