function [fc,Tp,B,lambda,k,Na,Tr,PRF,Nr,fs,r,xm,Tcoh,Tar_xyz]=echo_generate()
 %% 参数设置
 c=299792458;
 R=1737400;
 Tcoh=20;
 PRF=20;
 fc=3.25e9;
 Tp=3e-4;    %脉冲持续时间
 B=2e6;
 k=B/Tp;
 fs=4*B;
 Tr=1/PRF;
 Na=round(PRF*Tcoh);
 Na=Na+mod(Na,2);
 Tcoh=Na*Tr;
 lambda=c/fc;


 
%% 站心坐标系到月固坐标系转换 
trace_ref=[0,0]*pi/180; %星历参考点
str='Target-Target1-To-Target-SJZ AER_2021-10-08T05-30-00Z_2021-10-08T06-30-00Z-1s.csv';
interal=1;              %星历更新间隔
pos_xyz=sc2mc(trace_ref,str,interal,Tr,Na);

%% 航迹误差
% st1=(0:Na-1)*Tr;
% p_x=polyfit(st1,pos_xyz(1,:),3);
% p_y=polyfit(st1,pos_xyz(2,:),3);
% p_z=polyfit(st1,pos_xyz(3,:),3);
% 
% % p_x(end-2)=p_x(end-2)+0.2;
% % p_y(end-2)=p_y(end-2)+0.1;
% % p_z(end-2)=p_z(end-2)-0.2;
% 
% pos_xyz1(1,:) = polyval(p_x,st1);
% pos_xyz1(2,:) = polyval(p_y,st1);
% pos_xyz1(3,:) = polyval(p_z,st1);

%% 设置测绘区域
TarNum=1;
%小场景
% Tar=[-11,-43;
%      -16,-43;
%      -6,-43;
%      -11,-39;
%      -11,-47;
%             ]*pi/180;


Tar=[23,0;
     19,0;
     27,0;
     23,-4;
     23,4;
            ]*pi/180;



Tar_xyz=(R)*[cos(Tar(:,2)).*cos(Tar(:,1)),cos(Tar(:,2)).*sin(Tar(:,1)),sin(Tar(:,2))];


Rm=zeros(Na,TarNum);
for n=1:TarNum
    Rm(:,n)=sqrt((pos_xyz(1,:)-Tar_xyz(n,1)).^2+(pos_xyz(2,:)-Tar_xyz(n,2)).^2+(pos_xyz(3,:)-Tar_xyz(n,3)).^2);
end
Rmax=max(max(Rm));
Rmin=min(min(Rm));


 %% 产生回波
disp(['生成回波信号'])
win1=Tp;
win2=Tp;
Nr=round((2*Rmax/c-2*Rmin/c+win1+win2)*fs);
%Nr=2^nextpow2(Nr);
tn=linspace(-win1+2*Rmin/c,win2+2*Rmax/c,Nr);
r=tn*c/2;
%fs=Nr/(2*Rmax/c-2*Rmin/c+win1+win2);


%scapha=normrnd(0,1,1,TarNum);

xm=zeros(Na,Nr);
 for m=1:TarNum
     td=2*Rm(:,m)/c;
     td_tn=ones(Na,1)*tn-td*ones(1,Nr);
     x=rectpuls(td_tn-Tp/2,Tp).*exp(1i*pi*k*(td_tn-Tp/2).^2);
     xm=xm+x.*(exp(-2j*pi*fc*td)*ones(1,Nr));
 end

%% 投影
% R_ref=Rm(Na/2,1);
% RD_proj(Tar_xyz,pos_xyz,R_ref);
