%########################     DEM截取    #############################%
%——————————————————————更新状态——————————————————————%
% 版 本 号：1.0
% 近期修改：无
% 执 行 人：张光伟
% 审 核 人：无
%————————————————————————————————————————————————%
%——————————————————————函数介绍——————————————————————%
% 函数功能： DEM截取，这里未进行插值，500米间隔
%
%————————————————————————————————————————————————%
%                    参数名             参数定义              数据尺寸          数据类型
% 输入参数：          DEM_file_path      文件路径               
%                    lon_lat_range       经纬度范围
%————————————————————————————————————————————————%
%                    参数名             参数定义              数据尺寸          数据类型
% 输出参数：          DEM_crop         截取的DEM    
%————————————————————————————————————————————————%




function [DEM_crop]=DEM_crop_fun(DEM_file_path,lon_lat_range)
R=1737400;
%% 读取下载的DEM数据
%[DEM_moon,info] = geotiffread('globalDEM.tif');
load(DEM_file_path);
[row,col] = size(DEM_moon);  

L=linspace(-180,180,col);  %经纬度坐标
B=linspace(90,-90,row);
% figure,imagesc(L,B,DEM_moon)
% set(gca,'YDir','normal');

latcell=180/row;
loncell=360/col;
%% 确定区域边界经纬度范围  经度范围 从左至右；维度范围，从上至下
%lon_lat_range=[-16,-6,-39,-47];
%lon_lat_range=[19,27,4,-4];
%lon_lat_range=[16,22,-49,-53];


lat_start = (90-lon_lat_range(3))/180;     %经纬度起始，归一化
lat_end =(90-lon_lat_range(4))/180;
lon_start = abs(-180-lon_lat_range(1))/360;
lon_end =  abs(-180-lon_lat_range(2))/360;


%计算区域范围
lat_range=abs(lon_lat_range(4)-lon_lat_range(3))*pi/180*R;
RR=R*cos((lon_lat_range(3)+lon_lat_range(4))/2*pi/180);
lon_range=abs(lon_lat_range(1)-lon_lat_range(2))*pi/180*RR;

row_start = max(ceil(lat_start*row),1);
row_end = min(ceil(lat_end*row),row);
col_start = max(ceil(lon_start*col),1);
col_end = min(ceil(lon_end*col),col);
DEM_crop = DEM_moon(row_start:row_end, col_start:col_end);
DEM_crop=double(DEM_crop);
end