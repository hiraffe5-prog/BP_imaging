%% 单点目标二维基带回波仿真 —— 完整代码
clear; clc; close all;

%% ———— 参数定义 —————————————————————————————
c     = 3e8;                % 光速 [m/s]
lambda= 0.05;            % 波长 [m]
fc    = c/lambda;           % 载波频率 [Hz]
H     = 3000;               % 飞行高度 [m]
% Adeg  = 1.5;                % 波束全宽 [deg]
% A     = Adeg*pi/180;        % 波束全宽 [rad]
A     = 0.025;        % 波束全宽 [rad]
Rg0   = 15200;              % 场景中心地面距离 [m]
Xo    = 100;                % 距离向半宽 [m]
Yo    = 50;                 % 方位向半宽 [m]
Tp    = 20e-6;              % 脉冲宽度 [s]
B     = 70e6;               % 发射带宽 [Hz]
fs    = 500e6;              % 距离向采样率 [Hz]
PRF   = 800;                % 脉冲重复频率 [Hz]
Fa    = PRF;                % 方位向采样率
v0    = 100;                % 平台速度 [m/s]
Kr    = B / Tp;             % 距离调频率 [Hz/s]
A0    = 1;                  % 目标反射幅度

%% ———— 几何量计算 —————————————————————————————
% 中心斜距
R0   = sqrt(Rg0^2 + H^2);
Ka  = 2*v0^2/(lambda*R0);   % 方位向调频率

% 合成孔径长度 & 合成孔径慢时间
Ls   = A * R0;              
Ta   = Ls / v0;             
% 方位向带宽 B_a
Ba    = Ka * Ta;            % =2*v0*A/lambda
pa=v0/Ba;                   % 方位向分辨率
% 最近/最远斜距
Rmin = sqrt((Rg0 - Xo)^2 + H^2);
Rmax = sqrt((Rg0 + Xo)^2 + H^2 + (Ls/2)^2);

% 平台沿方位向需行走的总路程
rm   = Ls + 2*Yo;

%% ———— 时间轴离散 —————————————————————————————
% 慢时间 tm：覆盖合成孔径和场景
tm = 0 : 1/PRF : rm/v0 - 1/PRF;     
y  = -rm/2 + v0 * tm;               % 对应方位坐标 [m]

% 快时间 tr：覆盖距离门
tr = (2*Rmin/c - Tp/2) : 1/fs : (2*Rmax/c + Tp/2 - 1/fs);

Na = numel(tm);
Nr = numel(tr);

%% ———— 目标位置定义 ————————————————————————————
% 单点目标就在场景中心地面投影
targetX = Rg0;    % 地面距离向坐标 [m]
targetY = 0;      % 方位向坐标 [m]

%% ———— 回波生成 —————————————————————————————
s_echo = zeros(Na, Nr);

for ia = 1:Na
    % 瞬时斜距
    R_ta = sqrt( targetX^2 + (targetY - y(ia))^2 + H^2 );
    % 时延
    tau  = tr - 2*R_ta/c;
    % 距离门限
    wr   = rectpuls(tau, Tp);
    % 方位门限
    waz  = double( abs(targetY - y(ia)) < Ls/2 );
    % 回波：载波相位 + LFM 相位
    s_echo(ia,:) = A0 ...
                   * waz ...
                   .* wr ...
                   .* exp(-1j*4*pi*fc*R_ta/c) ...
                   .* exp(1j*pi*Kr*tau.^2);
end

%% ———— 坐标转换 & 绘图 ————————————————————————
% 地面距离向坐标
r = sqrt( (tr*c/2).^2 - H^2 );
[Rg,Ya] = meshgrid(r, y);

% figure;
% pcolor(r, y, real(s_echo));
% shading flat;
% colormap(gray);
% xlabel('距离向 (m)');
% ylabel('方位向 (m)');
% title('单点目标基带回波实部灰度图');

%% 距离压缩

% 距离频率轴（和 MATLAB 的FFT顺序一致）
f_tau = (-Nr/2:Nr/2-1) * (fs/Nr);         % Hz

% 距离向匹配滤波器 Hr(f_tau) = rect(f_tau/B) * exp(j*pi*f_tau^2/Kr)
Hr = (abs(f_tau) <= B/2) .* exp(1j*pi*(f_tau.^2)/Kr);
Hr = fftshift(Hr);                        % 变到与fft一致的顺序

% 距离向FFT -> 乘Hr -> IFFT
S_fd   = fft(s_echo, [], 2);              % 沿距离向做FFT
S_mfd  = ifft( S_fd .* repmat(Hr, size(s_echo,1), 1), [], 2 );
signal_matched = S_mfd;                   % 距离压缩结果（线性卷积等效）

% figure;
% pcolor(r, y, abs(signal_matched));
% shading flat;     % 去掉网格线
% colormap(gray);
% xlim([Rg0-Xo, Rg0+Xo]);
% xlabel('距离向 (m)');
% ylabel('方位向 (m)');
% title('距离压缩后回波数据');

%% —— 距离维插值（频域补零，L=10）——仅到插值与绘图 —— 
L = 4;                        % 升采样倍数
Nr_i = L * Nr;                 % 插值后的距离向点数

% 1) 频域补零（与参考代码一致的拆分点）
S_fd = fft(signal_matched, [], 2);              % Na×Nr
mid  = floor((Nr + 1)/2);                       % 拆分点：1…mid | mid+1…end
S_fd_pad = [ S_fd(:,1:mid), ...
             zeros(Na, (L-1)*Nr), ...
             S_fd(:,mid+1:end) ];               % Na×(L*Nr)

% 2) IFFT 得到距离插值后的回波
signal_interp = ifft(S_fd_pad, Nr_i, 2);        % Na×(L*Nr)

% 3) 构造插值后的快时间/地距坐标并绘图
tr_i = linspace(min(tr), max(tr), Nr_i);        % 与参考一致：同时间跨度、更细采样
r_i  = sqrt( (tr_i * c/2).^2 - H^2 );           % 地面距离坐标（与原流程一致）

% figure;
% pcolor(r_i, y, abs(signal_interp));             % 或 imagesc(r_i, y, 255-abs(signal_interp))
% shading flat; colormap(gray);
% xlim([Rg0 - Xo, Rg0 + Xo]);
% xlabel('距离向 (m)');
% ylabel('方位向 (m)');
% title(sprintf('距离插值 (L = %d) 后的距离压缩回波', L));

% % 1) 选一个 PRT （慢时间）索引
% ia0 = round(size(signal_interp,1)/2);  % 用合成孔径中心这一行
% 
% %% 2) 提取幅度并归一化转 dB
% vec     = abs( signal_interp(ia0,:) );    % 幅度
% vec     = vec / max(vec);                   % 归一化到 0…1
% vec_db  = 20*log10( vec );                  % 转为 dB
% 
% % 对应的斜距坐标（直接用 slant range）
% range_slant = tr_i * c/2;   
% %% 3) 找 –3 dB 主瓣宽度（采样点数和时间）
% % 找到峰值索引（应当是 0 dB 对应）
% [~, idx0]     = max(vec_db);
% R_peak = sqrt((tr_i(idx0)*c/2)^2 - H^2);
% fprintf('地距峰值 = %.3f m, 与 Rg0 差 %.3f m\n', R_peak, R_peak - Rg0);
% % 4) 在左右两侧各找第一个穿过 -3 dB 的点
% % 左侧
% j = find(vec_db(1:idx0) < -3, 1, 'last');
% x1 = interp1(vec_db(j:j+1), range_slant(j:j+1), -3);  
% % 右侧
% k = idx0 - 1 + find(vec_db(idx0:end) < -3, 1, 'first');
% x2 = interp1(vec_db(k-1:k), range_slant(k-1:k), -3);
% 
% % 5) 斜距分辨力
% range_res = x2 - x1;          
% 
% % 6) 绘图
% figure; hold on;
% plot(range_slant, vec_db, 'b-', 'LineWidth',1.5);
% % –3 dB 水平线
% yline(-3, 'r--', '-3 dB');
% % 两个标记点
% plot([x1 x2], [-3 -3], 'ro', 'MarkerFaceColor','r');
% % 垂直线标示
% xline(x1, 'r:');
% xline(x2, 'r:');
% xlim([15450, 15550]);
% ylim([-40 0]);
% xlabel('斜距 R (m)');
% ylabel('幅度 (dB)');
% 
% title(sprintf('压缩脉冲主瓣 (\\DeltaR = %.2f m)', range_res));
% grid on;
% legend('脉压曲线','–3 dB','交叉点','Location','Best');
% 
% fprintf('–3 dB 交叉点在 R = %.3f m 和 R = %.3f m\n', x1, x2);
% fprintf('对应斜距分辨力 ΔR = %.3f m\n', range_res);
% 
% %% ===== 1) 找主峰位置 =====
% [~, idx0] = max(vec);
% 
% %% ===== 2) 找主峰左右最近的“极小值”作为主瓣边界 =====
% % 为防止噪声导致假极小，可对“寻找极小值用的曲线”做轻微平滑，
% % 但能量一律在原曲线上算（不要用平滑后的）。
% az_for_min = movmean(vec, 3);           % 轻微平滑，只用于找极小值
% 
% % 用 islocalmin 找极小值（若版本不支持，可改用 findpeaks(-az_for_min)）
% isMin  = islocalmin(vec,'FlatSelection','center');
% minIdx = find(isMin);
% 
% leftMin  = minIdx(minIdx < idx0);
% rightMin = minIdx(minIdx > idx0);
% 
% if isempty(leftMin) || isempty(rightMin)
%     error('未在主峰两侧检测到极小值，请扩大方位范围或加大平滑窗口。');
% end
% 
% iL = leftMin(end);      % 主峰左侧最近的极小值
% iR = rightMin(1);       % 主峰右侧最近的极小值
% 
% % 为避免将极小值点本身纳入主瓣（通常幅度接近零），可适度收缩 1 个样点
% iL_use = iL; 
% iR_use = iR;
% 
% %% ===== 3) 计算 ISLR（积分旁瓣比） =====
% P = vec.^2;                     % 功率
% E_total = sum(P);                    % 总能量（等间隔采样时比例系数可省）
% idx_main = iL_use:iR_use;            % 主瓣样点索引
% E_main  = sum(P(idx_main));          % 主瓣能量
% E_side  = E_total - E_main;          % 旁瓣总能量
% ISLR_dB = 10*log10(E_side / E_main);
% 
% %% ===== 4) 打印与（可选）标注信息 =====
% fprintf('主峰索引 = %d，对应方位 y = %.4f m\n', idx0, range_slant(idx0));
% fprintf('主瓣边界（极小值法）：左 = %.4f m，右 = %.4f m\n', range_slant(iL), range_slant(iR));
% fprintf('主瓣样点区间用于积分：[%d, %d]\n', iL_use, iR_use);
% fprintf('ISLR = %.3f dB\n', ISLR_dB);
% 

%% —— 网格剖分（±50 区间）——
az_center = 0;                                         % 场景中心方位（与你的 y=0 对齐）
az_grid = (az_center - 50) : 0.05 : (az_center + 50);  % 方位网格
rg_grid = (Rg0 - 50) : 0.1  : (Rg0 + 50);              % 地距网格
Na_pix  = numel(az_grid);
Nr_pix  = numel(rg_grid);

%% —— 预计算 —— 
Na        = numel(y);                % 慢时间脉冲数（行数）
dt_i      = tr_i(2) - tr_i(1);       % 插值后时间分辨率
img_bp    = zeros(Na_pix, Nr_pix);   % BP复值图（方位×地距）

%% —— BP：对每个像素，沿距离徙动曲线取样+相位补偿+相参累加 —— 
for iaz = 1:Na_pix
    az_pix = az_grid(iaz);

    % 天线方位门限（矩形窗；也可换成加权，如cos或Hamming）
    in_beam = abs(az_pix - y) < (Ls/2);              % 逻辑向量，长度 Na_pulses
    pulse_idx = find(in_beam);                       % 参与累加的脉冲下标
    if isempty(pulse_idx)
        continue;                                    % 此像素不被照明，跳过
    end

    % 为了少重复索引，加个本地变量
    y_used = y(pulse_idx);

    for jrg = 1:Nr_pix
        rg_pix = rg_grid(jrg);

        % 该像素对所有参与脉冲的瞬时斜距 R(ia) 与回波到达时延 tau(ia)
        R   = sqrt( rg_pix.^2 + (az_pix - y_used).^2 + H^2 );  % 1×Npulses
        tau = 2 * R / c;

        % 把 tau 映射到 signal_interp 的列索引（最近邻）
        idx = round( (tau - tr_i(1)) / dt_i ) + 1;             % 1×Npulses
        valid = (idx >= 1) & (idx <= Nr_i);                    % 边界保护
        if ~any(valid), continue; end

        % 从 signal_interp(row=脉冲索引, col=对应距离索引) 取样
        lin_id = sub2ind([Na, Nr_i], pulse_idx(valid), idx(valid));
        rd     = signal_interp(lin_id);

        % 载频相位补偿（把 e^{-j 4π fc R/c} 拉回同相位以便相参叠加）
        rd = rd .* exp(1j * 4*pi*fc * R(valid) / c);

        % （可选）距离权重/余弦加权等：
        % w = ones(size(rd)); % 或者 cos 窗、Hamming 窗等
        % rd = rd .* w;

        % 相参累加
        img_bp(iaz, jrg) = sum(rd);

        % （可选）做个能量归一化，避免像素照明脉冲数不等导致亮度不均：
        % img_bp(iaz, jrg) = img_bp(iaz, jrg) / nnz(valid);
    end
end

% 结果矩阵：img_bp (方位×地距)，与 (az_grid, rg_grid) 对应
% 
% figure;
% mesh(rg_grid,az_grid,abs(img_bp));                       
% xlabel('地距/m');
% ylabel('方位/m');
% title('SAR图像(插值后）');
% colormap(gray)

% figure;
% mesh(rg_grid,az_grid,abs(img_bp)); 
% view(0,90);
% % xlim([Rg0-Xo, Rg0+Xo]);
% ylim([-10 10]); 
% xlabel('距离向 (m)');
% ylabel('方位向 (m)');
% title('SAR成像结果');

% —— 这里假定你已有变量：img_bp (Na_pix×Nr_pix), az_grid, rg_grid, Rg0, lambda, R0, Ls ——

% 1) 用网格的地距轴来找最靠近 Rg0 的列
[~, ir0] = min(abs(rg_grid - Rg0));      % ir0 是 img_bp 的列号（<= length(rg_grid)）

% 2) 提取这一列的幅度，并归一化 + 转为 dB
az_axis       = az_grid(:);              % 方位像素坐标 (Na_pix×1)
az_profile    = abs(img_bp(:, ir0));     % Na_pix×1
az_profile    = az_profile / max(az_profile + eps);
az_profile_db = 20*log10(az_profile + eps);

% 3) 找 –3 dB 主瓣宽度（加入边界保护）
[~, idx_peak] = max(az_profile_db);      % 峰值索引（在方位像素轴上）

% 左侧第一个穿过 -3 dB 的点
j = find(az_profile_db(1:idx_peak) <= -3, 1, 'last');
if ~isempty(j) && j < idx_peak
    y_left = interp1(az_profile_db(j:j+1), az_axis(j:j+1), -3);
else
    y_left = NaN;
end

% 右侧第一个穿过 -3 dB 的点
k_rel = find(az_profile_db(idx_peak:end) <= -3, 1, 'first');
if ~isempty(k_rel) && (idx_peak + k_rel - 1) > idx_peak
    k = idx_peak + k_rel - 1;
    y_right = interp1(az_profile_db([k-1 k]), az_axis([k-1 k]), -3);
else
    y_right = NaN;
end

% 方位向分辨率（可能出现 NaN，视主瓣是否完整穿越 -3 dB）
az_resolution = y_right - y_left;

% 理论方位分辨率（与你前面推导等价：0.886*lambda/(2*A)；用 Ls=A*R0 写成下面形式）
az_res_theory = 0.886 * lambda * R0 / (2 * Ls);
% 4) 绘图（横轴用像素方位轴）
figure; hold on;
plot(az_axis, az_profile_db, 'b-', 'LineWidth', 1.5);
yline(-3, 'r--', '-3 dB');
if ~isnan(y_left) && ~isnan(y_right)
    plot([y_left y_right], [-3 -3], 'ro', 'MarkerFaceColor', 'r');
end
grid on;
xlabel('方位向 (m)');
ylabel('幅度 (dB)');
title(sprintf('方位向压缩主瓣 (r \\approx %.1f m)', rg_grid(ir0)));
xlim([-10 10]);           % 需要可自行调整
ylim([-60 0]);
legend('脉压曲线','–3 dB','交叉点','Location','Best');
% 5) 打印结果
if ~isnan(az_resolution)
    fprintf('–3 dB 交叉点在 y = %.3f m 和 y = %.3f m\n', y_left, y_right);
    fprintf('测量方位分辨力 Δy = %.3f m\n', az_resolution);
    fprintf('理论方位分辨力 Δy_theory = %.3f m\n', az_res_theory);
    fprintf('误差 = %.3f m (%.2f%%)\n', az_resolution - az_res_theory, ...
        (az_resolution - az_res_theory)/az_res_theory*100);
else
    fprintf('未能找到完整的 -3 dB 交叉点，请检查主瓣是否在显示窗口内或提升动态范围。\n');
end

% 主瓣位置（在像素方位轴上）
[~, ia_pk] = max(az_profile);
fprintf('主瓣位于 y = %.3f m\n', az_axis(ia_pk));

%% ===== 1) 找主峰位置 =====
[~, idx0] = max(az_profile);

%% ===== 2) 找主峰左右最近的“极小值”作为主瓣边界 =====
% 为防止噪声导致假极小，可对“寻找极小值用的曲线”做轻微平滑，
% 但能量一律在原曲线上算（不要用平滑后的）。
az_for_min = movmean(az_profile, 3);           % 轻微平滑，只用于找极小值

% 用 islocalmin 找极小值（若版本不支持，可改用 findpeaks(-az_for_min)）
isMin  = islocalmin(az_profile,'FlatSelection','center');
minIdx = find(isMin);

leftMin  = minIdx(minIdx < idx0);
rightMin = minIdx(minIdx > idx0);

if isempty(leftMin) || isempty(rightMin)
    error('未在主峰两侧检测到极小值，请扩大方位范围或加大平滑窗口。');
end

iL = leftMin(end);      % 主峰左侧最近的极小值
iR = rightMin(1);       % 主峰右侧最近的极小值

% 为避免将极小值点本身纳入主瓣（通常幅度接近零），可适度收缩 1 个样点
iL_use = min(iL+1, idx0); 
iR_use = max(iR-1, idx0);

%% ===== 3) 计算 ISLR（积分旁瓣比） =====
P = az_profile.^2;                   % 功率
E_total = sum(P);                    % 总能量（等间隔采样时比例系数可省）
idx_main = iL_use:iR_use;            % 主瓣样点索引
E_main  = sum(P(idx_main));          % 主瓣能量
E_side  = E_total - E_main;          % 旁瓣总能量
ISLR_dB = 10*log10(E_side / E_main);

%% ===== 4) 打印与（可选）标注信息 =====
fprintf('主峰索引 = %d，对应方位 y = %.4f m\n', idx0, az_axis(idx0));
fprintf('主瓣边界（极小值法）：左 = %.4f m，右 = %.4f m\n', az_axis(iL), az_axis(iR));
fprintf('主瓣样点区间用于积分：[%d, %d]\n', iL_use, iR_use);
fprintf('ISLR = %.3f dB\n', ISLR_dB);
% 
