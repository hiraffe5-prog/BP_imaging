%########################   DEM到RD图      #############################%
%########################         Arthur: zgw        #############################%
%########################          Date：2021.6.9         #############################%
%########################     DEM截取    #############################%
%——————————————————————更新状态——————————————————————%
% 版 本 号：1.0
% 近期修改：无
% 执 行 人：张光伟
% 审 核 人：无
%————————————————————————————————————————————————%
%——————————————————————函数介绍——————————————————————%
% 函数功能： 基于RD图的反向地理解码
%
%————————————————————————————————————————————————%
%                    参数名             参数定义              数据尺寸          数据类型
% 输入参数：          pos_fit           补偿时刻，定为Na/2           
%                    lon_lat_range      经纬度范围,deg
%                    DEM                高程
%                    N_DEM_interp       插值系数
%                    Tar_ref            场景参考点，弧度，与成像参考一致，注意是否带高度，一般不带
%                    master_r           图像距离轴
%                    master_fd          图像多普勒轴
%                    y_rd               RD图
%                    pos_xyz            平台月固坐标
%                    Tr,Na,fc           脉冲重复周期、脉冲数、载频
%————————————————————————————————————————————————%
%                    参数名             参数定义              数据尺寸          数据类型
% 输出参数：          y_geo_samp        地理解码图像   
%                    DEM_interp         DEM插值结果
%                     lon               经度坐标轴
%                     lat               纬度坐标轴
%————————————————————————————————————————————————%


function [y_rd_ref] = RD_geocoding3(y_bp,pos_fit,lon_lat_range,del_lon,del_lat,DEM,pos_xyz,Tar_ref_xyz,master_r,master_fd,Tr,Na,fc)
%pos_fit(补偿时刻，定为Na/2),lon_lat_range(经纬度范围),DEM,N_DEM_interp（插值倍数）,Tar_ref（参考点）,master_r,master_fd,y_rd（图像域）,...
%pos_xyz（月固航迹）,Tr,Na,fc
   %% 基本参数
    c = 299792458;
    lambda = c/fc;
    R=1737400;

    %% 航迹速度
    %pos_xyz=sc2mc(Track_refer,str,interval,Tr,Na);
    
    %__计算绝对速度矢量__%
    st=(0:Na-1)*Tr;
    px = polyfit(st,pos_xyz(1,1:Na),3);  %三维坐标对应的多项式系数
    py = polyfit(st,pos_xyz(2,1:Na),3);
    pz = polyfit(st,pos_xyz(3,1:Na),3);
    
    Vx = 3*px(1)*st(1,pos_fit)^2 + 2*px(2)*st(1,pos_fit) + px(3); %三维速度
    Vy = 3*py(1)*st(1,pos_fit)^2 + 2*py(2)*st(1,pos_fit) + py(3);
    Vz = 3*pz(1)*st(1,pos_fit)^2 + 2*pz(2)*st(1,pos_fit) + pz(3);
    
    velocity_tar = [Vx,Vy,Vz];
   

   
   %% DEM升采样
    % N_interp为插值精度
   
    lon = linspace(lon_lat_range(1), lon_lat_range(2), abs(lon_lat_range(2)-lon_lat_range(1))/del_lon);
    lat = linspace(lon_lat_range(3), lon_lat_range(4), abs(lon_lat_range(4)-lon_lat_range(3))/del_lat);
    Nrow = length(lat);
    Ncol = length(lon);
    
    Vm = zeros(Nrow,Ncol);          %声明速度矩阵（存）
    Rm = zeros(Nrow,Ncol);          %声明主图像斜距（存）
    %% 计算R，V,fd   
    for ii = 1:length(lat)
        ii/length(lat)
        for jj = 1:length(lon)
%             waitbar(((ii-1)*length(lon)+jj)/969/973/(N_interp^2) , wh1, ['斜距已计算 ', num2str(((ii-1)*length(lon)+jj)/969/973/(N_interp^2)*100, '%.2f'), '%' ]);
            Tar_x = (R + DEM(ii,jj))*cosd(lat(ii))*cosd(lon(jj));   %目标点对应的三维坐标
            Tar_y = (R + DEM(ii,jj))*cosd(lat(ii))*sind(lon(jj));
            Tar_z = (R + DEM(ii,jj))*sind(lat(ii)); 
            R_vector = [(pos_xyz(1,pos_fit) - Tar_x).',...
                (pos_xyz(2,pos_fit) - Tar_y).',(pos_xyz(3,pos_fit) - Tar_z).'];    %斜距    
            
            Vm(ii, jj) = dot(R_vector, velocity_tar) / norm(R_vector);     %雷达视线方向的速度投影
            Rm(ii, jj) = norm(R_vector);    

        end
    end
   
    %参考点信息,和成像保持一致，注意是否带高度，一般不带
    %Tar_ref_xyz=R*[cos(Tar_ref(2))*cos(Tar_ref(1)),cos(Tar_ref(2))*sin(Tar_ref(1)),sin(Tar_ref(2))]; %参考点三维坐标
    R_vector_ref = [(pos_xyz(1,pos_fit) - Tar_ref_xyz(1)).',...
                (pos_xyz(2,pos_fit) - Tar_ref_xyz(2)).',(pos_xyz(3,pos_fit) - Tar_ref_xyz(3)).'];
    Vm_ref= dot(R_vector_ref, velocity_tar) / norm(R_vector_ref);                                   %参考点速度
    fd_ref = - 2*Vm_ref/lambda;
    Rm_ref = norm(R_vector_ref);    
    
    %其他点信息
    fd = - 2*Vm/lambda;  % 绝对多普勒，适用于地理编码多普勒方程   
    fd = fd - fd_ref; % 相对多普勒，适用于重采样至距离多普勒平面
    
%     % 采样
%     disp('正在进行图像重采样')
%     % 重整投影范围坐标轴
%     % 计算投影后距离和多普勒最大值与最小值以及坐标轴间隔
%     Rm_max = max(max(Rm)); Rm_min = min(min(Rm));
%     fd_max = max(max(fd)); fd_min = min(min(fd));
%     delta_Rm = master_r(2)-master_r(1); delta_fd = master_fd(2)-master_fd(1);
%     % 构建新坐标轴
%     master_r = Rm_min:delta_Rm:Rm_max; master_fd = fd_min:delta_fd:fd_max;
%     %y_geo_samp = interp2(master_r, master_fd, y_rd, Rm, fd, 'spline');   
%     [ master_R, master_Fd] = meshgrid( master_r, master_fd);
%     % y_rd_ref = interp2(Rm(1,:), fd(:,1), y_bp, master_R, master_Fd, 'cubic'); 
%     % tmp = abs(y_bp); tmp = tmp-min(min(tmp)); tmp = tmp/max(max(tmp));
%     % y_rd_ref = griddata(Rm(:), fd(:),abs(y_bp(:)), master_R, master_Fd,'cubic');  
%     F = scatteredInterpolant(Rm(:),fd(:),abs(y_bp(:)),'linear');
%     y_rd_ref = F(master_R,master_Fd);
%     %'linear'、'nearest'、'cubic'、'makima' 或 'spline'
%     figure; imagesc(y_rd_ref)
%     disp('已完成图像重采样')
    % 绘制反投到RD域的bp成像结果
    figure; scatter(Rm(:),fd(:),[],abs(y_bp(:)),'.');

end
