function [Tranx,Trany,Tranz]=mc2dd(pos_xyz,nad)
%% 月固到延时多普勒 求沿x y z 旋转角度
nad1=[nad(1),nad(2),0];                           %x'在x'' y''平面（月固面）的投影
lambdasrp=acos(abs(nad1*[1 0 0].')/norm(nad1));   %投影与x''夹角
thetasrp=acos(abs(nad1*nad)/norm(nad1)/norm(nad));%投影与x'的夹角
if(nad(2)<0)                                     %若星下点在y''正侧，则由x''-x'需要绕z''逆时针转(角度为正值)
    lambdasrp=-lambdasrp;
end

if(nad(3)>0)                                     %若星下点在z''负侧，则由x''-x'需要绕y''逆时针转(角度为正值)      
    thetasrp=-thetasrp;
end
    
Trany=[cos(thetasrp) 0 -sin(thetasrp);0 1 0;sin(thetasrp) 0 cos(thetasrp)];
Tranz=[cos(lambdasrp) sin(lambdasrp) 0;-sin(lambdasrp) cos(lambdasrp) 0;0 0 1];

r1=Trany*Tranz*pos_xyz(:,1);%x'坐标系运动起点，这里算运行矢量，4分钟
r2=Trany*Tranz*pos_xyz(:,end);%x'坐标系运动终点
Run=r2-r1;%雷达运行矢量，便于求视旋转轴
z=cross(Run,[1,0,0]);%交叉积，垂直于两向量的向量,视轴方向

if(z(3)<=0)
    z=-z;
end    
gama=acos(abs(z*[0,0,1].')/norm(z));%，在x'y'z'坐标系下，z'方向[0,0,1],需要求z'到视轴z的旋转角度
if(z(2)>=0)
    gama=-gama;
end

Tranx=[1,0,0;0,cos(gama),sin(gama);0,-sin(gama),cos(gama)];
end
