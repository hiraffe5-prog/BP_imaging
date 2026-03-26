function [fc,Tp,B,lambda,k,Na,Tr,PRF,Nr,fs,r,xm,DEM,lon_lat_range]=echo_generate_widescene()
 %% 参数设置
 c=299792458;
 R=1737400;
 Tcoh= 120;
 PRF= 20;
 fc=3.1e9;
 Tp=1e-4;    %脉冲持续时间
 B= 0.5e6;
 k= B/Tp;
 fs= 2*B;
 c= 299792458;
 Tr= 1/PRF;
 Na=round(PRF*Tcoh);
 Na=Na+mod(Na,2);
 Tcoh=Na*Tr;
 lambda=c/fc;

 
%% 站心坐标系到月固坐标系转换
trace_ref=[0,0]*pi/180; %星历参考点
str='Target-Target1-To-Target-SJZ AER_2021-10-08T05-30-00Z_2021-10-08T06-30-00Z-1s.csv';
%str='Target-MOON-To-Target-SJZ AER_2022-05-14T11-00-00Z_2022-05-14T12-00-00Z-1s.csv';
%str='Target-Target1-To-Target-SJZ AER.csv';
interal=1;              %星历更新间隔
pos_xyz=sc2mc(trace_ref,str,interal,Tr,Na);

%% 航迹误差
% st1=(0:Na-1)*Tr;
% p_x=polyfit(st1,pos_xyz(1,:),3);
% p_y=polyfit(st1,pos_xyz(2,:),3);
% p_z=polyfit(st1,pos_xyz(3,:),3);

% p_x(end-2)=p_x(end-2)+0.2;
% p_y(end-2)=p_y(end-2)+0.1;
% p_z(end-2)=p_z(end-2)-0.2;

% pos_xyz1(1,:) = polyval(p_x,st1);
% pos_xyz1(2,:) = polyval(p_y,st1);
% pos_xyz1(3,:) = polyval(p_z,st1);


%% 设置测绘区域
%load('DEM_Tycho_500.mat');%包含DEM,lon_lat_range
%load('DEM_23_0_250m.mat');%包含DEM,lon_lat_range
%DEM=DEM_crop;
globalDEM_path = 'DEM_global_moon.mat';
lon_lat_range = [-13,-9,-41,-45];
del_lat = 0.02;
del_lon = 0.02;
lon = linspace(lon_lat_range(1), lon_lat_range(2), abs(lon_lat_range(2)-lon_lat_range(1))/del_lon);
lat = linspace(lon_lat_range(3), lon_lat_range(4), abs(lon_lat_range(4)-lon_lat_range(3))/del_lat);

[DEM_1]=DEM_crop_fun(globalDEM_path,lon_lat_range); 
[DEM]=DEM_resample(DEM_1,lon_lat_range,lon,lat);

figure,imagesc(lon,lat,DEM)
set(gca,'YDir','normal')
xlabel('经度/deg'),ylabel('纬度/deg'),title('高程图')

%% 点目标设置
 [row,col]=size(DEM);

 del_lat=abs(lon_lat_range(3)-lon_lat_range(4))/row;
 del_lon=abs(lon_lat_range(1)-lon_lat_range(2))/col;
 B_vec=(lon_lat_range(3):-del_lat:lon_lat_range(4)+del_lat)*pi/180;
 L_vec=(lon_lat_range(1):del_lon:lon_lat_range(2)-del_lon)*pi/180;
%B_vec=Tar_cen(2)+(0)*pi/180;
%L_vec=Tar_cen(1)+(5)*pi/180;
TarNum=length(B_vec)*length(L_vec);
Tar_xyz=zeros(TarNum,3);

h1=waitbar(0,"散射点位置计算中...");
for n=1:length(B_vec)
     waitbar(n/row);
    for m=1:length(L_vec)
        %Tar((n-1)*length(L_vec)+m,:)=[L_vec(m),B_vec(n)];
        Tar_xyz((n-1)*length(L_vec)+m,:)=(R+DEM(n,m))*...
            [cos(B_vec(n)).*cos(L_vec(m)),cos(B_vec(n)).*sin(L_vec(m)),sin(B_vec(n))];
    end
end
close(h1);
[incidence,sigma]=scatterpara(R,pos_xyz,DEM,L_vec,B_vec);
%sigma(round(row/2),round(col/2))=100*sigma(round(row/2),round(col/2));
sigma=sigma.';
sigma=sigma(:);
figure,imagesc(L_vec*180/pi,B_vec*180/pi,abs(cos(incidence)))
xlabel('经度/deg'),ylabel('纬度/deg')
set(gca,'YDir','normal')

%% 斜距计算
disp(['斜距计算中...'])

%加误差
% st=(0:Na-1)*Tr-Na*Tr/2;
% r_error = polyval([0.0000,0.01,0e3],st);
% r_error_p = polyfit(st,r_error,2);

Rm=zeros(Na,TarNum);
%Rm(:,1)=trace(:,3);
for n=1:TarNum
    n/TarNum
   % Rm(:,n)=sqrt((pos_xyz(1,:)-Tar_xyz(n,1)).^2+(pos_xyz(2,:)-Tar_xyz(n,2)).^2+(pos_xyz(3,:)-Tar_xyz(n,3)).^2)+r_error;
    Rm(:,n)=sqrt((pos_xyz(1,:)-Tar_xyz(n,1)).^2+(pos_xyz(2,:)-Tar_xyz(n,2)).^2+(pos_xyz(3,:)-Tar_xyz(n,3)).^2);

end

Rmax=max(max(Rm));
Rmin=min(min(Rm));
 
 %% 产生回波
disp(['生成回波信号'])

win1=0.1*Tp;
win2=Tp;
Nr=round((2*Rmax/c-2*Rmin/c+win1+win2)*fs);
%Nr=2^nextpow2(Nr);
Nr = Nr + mod(Nr,2);
tn=linspace(-win1+2*Rmin/c,win2+2*Rmax/c,Nr);
fs=Nr/(2*Rmax/c-2*Rmin/c+win1+win2);
r=tn*c/2;


xm=zeros(Na,Nr);
%Rmg=gpuArray(Rm);
tng=gpuArray(ones(Na,1)*tn);

%point_pha=randn(1,TarNum);
point_pha=zeros(1,TarNum);

 for m=1:TarNum
     m/TarNum
     td=gpuArray(2*Rm(:,m)/c);
     td_tn=tng-td*ones(1,Nr);
     x=rectpuls(td_tn-Tp/2,Tp).*exp(1i*pi*k*(td_tn-Tp/2).^2);
     xm=xm+sigma(m)*x.*(exp(-2j*pi*fc*td).*exp(2j*pi*point_pha(m))*ones(1,Nr));
 end
 raw1=gather(xm);
 xm=raw1;
 Para.xm=xm;
 Para.fs=fs;
 Para.Nr=Nr;
 Para.Na=Na;
 Para.tn=tn;
 Para.r=r;
 Para.pos_xyz=pos_xyz;
 Para.DEM=DEM;
 Para.lon_lat_range=lon_lat_range;
disp(['已生成回波信号'])

end
