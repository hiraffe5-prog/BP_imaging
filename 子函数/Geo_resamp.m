function [jingdu_mesh,weidu_mesh,sca_mesh]=Geo_resamp(lon_lat_range,del_lon,del_lat,L_mcmf,B_mcmf,result_glo)

down_sample_jingdu_line=L_mcmf(:);
down_sample_weidu_line=B_mcmf(:);
down_sample_height_line=abs(result_glo(:));

%通过插值方法获得均匀平面网格上的数据
F = TriScatteredInterp(down_sample_jingdu_line,down_sample_weidu_line,down_sample_height_line);
%F = scatteredInterpolant(down_sample_jingdu_line.',down_sample_weidu_line.',down_sample_height_line);
% 注意 : 纬度从上到下减小
%ylin = linspace(max(down_sample_weidu_line),min(down_sample_weidu_line),(max(down_sample_weidu_line)-min(down_sample_weidu_line))/del_lat);
ylin = linspace(lon_lat_range(3),lon_lat_range(4),abs(lon_lat_range(3)-lon_lat_range(4))/del_lat);

% 注意 : 经度从左到右增大 
%xlin = linspace(min(down_sample_jingdu_line),max(down_sample_jingdu_line),(max(down_sample_jingdu_line)-min(down_sample_jingdu_line))/del_lon);
xlin = linspace(lon_lat_range(1),lon_lat_range(2),abs(lon_lat_range(1)-lon_lat_range(2))/del_lon);

X_num = length(xlin);
Y_num = length(ylin);
mesh_size = [Y_num,X_num];
% 重采样网格
[jingdu_mesh,weidu_mesh] = meshgrid(xlin,ylin);
sca_mesh = zeros(Y_num,X_num);

% 逐行插值
% 注意Xlin向量方向
for ii = 1:Y_num
    ii
    YY_line = ones(1,X_num)*ylin(ii);
    sca_mesh(ii,:) = F(xlin,YY_line);
end


% figure,imagesc(xlin,ylin,sca_mesh);
% set(gca,'YDir','normal')
% %figure,contour(xlin,ylin,abs(sca_mesh),20);
% ylabel('纬度/度'),xlabel('经度/度'),title('成像结果')
% col=[0.3 0.3 0.3]+gray(128)*0.7;
% col(1,:)=[0,0,0];
% colormap(col)

end