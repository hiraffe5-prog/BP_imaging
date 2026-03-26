%修改了斜距计算方式，防止点数过多内存溢出
function [fc,Tp,B,lambda,k,Na,Tr,PRF,Nr,fs,r,xm,DEM,lon_lat_range]=echo_generate_widescene_biger()
tic 
%% 参数设置
 c=299792458;
 R=1737400;
 Tcoh=300;
 PRF=20;
 fc=3.25e9;
 Tp=3e-4;    %脉冲持续时间
 B=2e6;
 k=B/Tp;
 fs=4e6;
 c=299792458;
 Tr=1/PRF;
 Na=round(PRF*Tcoh);
 Na=Na+mod(Na,2);
 lambda=c/fc;
 v=lambda*PRF/2;%最大多普勒速度，盲速

 
%% 站心坐标系到月固坐标系转换
trace_ref=[0,0]*pi/180; %星历参考点
str='Target-Target1-To-Target-SJZ AER_2021-10-08T05-30-00Z_2021-10-08T06-30-00Z-1s.csv';
%str='Target-Target1-To-Target-SJZ AER.csv';
interal=1;              %星历更新间隔
pos_xyz=sc2mc(trace_ref,str,interal,Tr,Na);

%% 航迹误差
st1=(0:Na-1)*Tr;
p_x=polyfit(st1,pos_xyz(1,:),3);
p_y=polyfit(st1,pos_xyz(2,:),3);
p_z=polyfit(st1,pos_xyz(3,:),3);

% p_x(end-2)=p_x(end-2)+0.2;
% p_y(end-2)=p_y(end-2)+0.1;
% p_z(end-2)=p_z(end-2)-0.2;

pos_xyz1(1,:) = polyval(p_x,st1);
pos_xyz1(2,:) = polyval(p_y,st1);
pos_xyz1(3,:) = polyval(p_z,st1);


%% 设置测绘区域
%load('DEM_Tycho_500.mat');%包含DEM,lon_lat_range
load('DEM_23_0_250m.mat');%包含DEM,lon_lat_range
%DEM=DEM_crop;
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
sigma(round(row/2),round(col/2))=200*sigma(round(row/2),round(col/2));
sigma=sigma.';
sigma=sigma(:);
figure,imagesc(L_vec*180/pi,B_vec*180/pi,abs(cos(incidence)))

%% 斜距计算
disp(['斜距计算中...'])

%计算最大最小斜距
n=1;
Rm=sqrt((pos_xyz(1,:)-Tar_xyz(n,1)).^2+(pos_xyz(2,:)-Tar_xyz(n,2)).^2+(pos_xyz(3,:)-Tar_xyz(n,3)).^2);
Rmax=max(Rm);
Rmin=min(Rm);

for n=1:TarNum
    n/TarNum
    %Rm(:,n)=sqrt((pos_xyz(1,:)-Tar_xyz(n,1)).^2+(pos_xyz(2,:)-Tar_xyz(n,2)).^2+(pos_xyz(3,:)-Tar_xyz(n,3)).^2);
    Rm=sqrt((pos_xyz(1,:)-Tar_xyz(n,1)).^2+(pos_xyz(2,:)-Tar_xyz(n,2)).^2+(pos_xyz(3,:)-Tar_xyz(n,3)).^2);
    Rmax=max(Rmax,max(Rm));
    Rmin=min(Rmin,min(Rm));
end

% Rmax=max(max(Rm));
% Rmin=min(min(Rm));
 
 %% 产生回波
disp(['生成回波信号'])

win1=0.1*Tp;
win2=Tp;
Nr=round((2*Rmax/c-2*Rmin/c+win1+win2)*fs);
%Nr=2^nextpow2(Nr);
tn=linspace(-win1+2*Rmin/c,win2+2*Rmax/c,Nr);
%fs=Nr/(2*Rmax/c-2*Rmin/c+win1+win2);
r=tn*c/2;


xm=zeros(Na,Nr);
%Rmg=gpuArray(Rm);
tng=gpuArray(ones(Na,1)*tn);

%point_pha=randn(1,TarNum);
point_pha=zeros(1,TarNum);

 for m=1:TarNum
     m/TarNum
     Rm=sqrt((pos_xyz(1,:)-Tar_xyz(m,1)).^2+(pos_xyz(2,:)-Tar_xyz(m,2)).^2+(pos_xyz(3,:)-Tar_xyz(m,3)).^2);
     td=gpuArray(2*Rm.'/c);
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
figure,imagesc(abs(xm))
toc
save testpara_1022  fc c B Tp Tr  PRF fs k lambda r xm
end
