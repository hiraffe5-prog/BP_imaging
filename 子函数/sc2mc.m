function pos_xyz=sc2mc(Track_refer,str,interal,Tr,Na) 
R=1737400;
Track_refer_xyz=R*[cos(Track_refer(2))*cos(Track_refer(1)),cos(Track_refer(2))*sin(Track_refer(1)),sin(Track_refer(2))];

pitch=-(pi/2-Track_refer(2));%俯仰向沿y轴顺时针转
azimuth=-Track_refer(1);

% 转换矩阵，Tx沿x轴转换
Ty=[cos(pitch) 0 -sin(pitch);0 1 0;sin(pitch) 0 cos(pitch)];
Tz=[cos(azimuth) sin(azimuth) 0;-sin(azimuth) cos(azimuth) 0;0 0 1];

% 月固坐标系轨道数据
%str='Target-Target1-To-Target-Qujing AER_0912_1s.csv';
%interal=1;  %更新间隔
trace=csvread(str,1,1);
len=size(trace,1);

trace=trace(:,end-2:end);
trace(:,3)=1000*trace(:,3);

%拟合

st0=(0:len-1)*interal;
st1=(0:Na-1)*Tr;
p=polyfit(st0,trace(1:len,2).',3);
trace1(2,:) = polyval(p,st1);

p=polyfit(st0,trace(1:len,1).',3);
trace1(1,:) = polyval(p,st1);

p=polyfit(st0,trace(1:len,3).',3);
trace1(3,:) = polyval(p,st1);

trace1(4,:)=trace1(3,:).*cos(trace1(2,:).*pi/180).*cos(pi-trace1(1,:).*pi/180);
trace1(5,:)=trace1(3,:).*cos(trace1(2,:).*pi/180).*sin(pi-trace1(1,:).*pi/180);
trace1(6,:)=trace1(3,:).*sin(trace1(2,:).*pi/180);
pos_xyz=bsxfun(@plus,(Tz*(Ty*trace1(4:6,:))),Track_refer_xyz.');

end
