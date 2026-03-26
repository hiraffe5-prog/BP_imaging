% BP对月成像处理
clc,clear,close all
%% 常数
R=1737400;
c=299792458;

%% 数据预处理
%[Tp,B,k,Na,Tr,Nr,fs,r,xm]=pre_process();
%[fc,Tp,B,lambda,k,Na,Tr,PRF,Nr,fs,r,xm]=echo_generate();
%[fc,Tp,B,lambda,k,Na,Tr,PRF,Nr,fs,r,S,Tcoh,Tar_xyz]=echo_generate_test();
[fc,Tp,B,lambda,k,Na,Tr,PRF,Nr,fs,r,xm,DEM,lon_lat_range]=echo_generate_widescene();

figure,imagesc(r,1:Na,abs(xm))
xlabel('距离/m'),ylabel('方位向(采样点)'),title('原始回波')

%% 读取星历
%STK
trace_ref=[0,0]*pi/180; %星历参考点
str='Target-Target1-To-Target-SJZ AER_2021-10-08T05-30-00Z_2021-10-08T06-30-00Z-1s.csv';
%str='Target-Target1-To-Target-SJZ AER.csv';
interal=1;              %星历更新间隔
pos_xyz=sc2mc(trace_ref,str,interal,Tr,Na);

%NASA
% str='Target-Target1-To-Target-Qujing AER_0912_1s.csv';
% pos_xyz=csvread(str,1,1);

%% 非停走误差,场景空变性小于波长，可直接在斜距上补偿



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


lon_range=[-13.5,-8.5];
lat_range=[-40.5,-45.5];
lon_lat_range=[lon_range,lat_range];
del_lat = 0.005;
del_lon = 0.005;
lon = linspace(lon_lat_range(1), lon_lat_range(2), abs(lon_lat_range(2)-lon_lat_range(1))/del_lon);
lat = linspace(lon_lat_range(3), lon_lat_range(4), abs(lon_lat_range(4)-lon_lat_range(3))/del_lat);
[lon2D,lat2D] = meshgrid(lon,lat);

globalDEM_path = 'DEM_global_moon.mat';
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



