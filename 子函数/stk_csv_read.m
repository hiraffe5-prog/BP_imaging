function [Shike_Line_Str,AER]=stk_csv_read(FIle_Name_Full1)
%实现STK原始数据的读取

fid1 = fopen(FIle_Name_Full1,'r');
data_hangji = textscan(fid1,'%s','delimiter',',');
fclose(fid1);

%%
% 坐标及区间转换
%
% 坐标时刻
Shike_Line = char(data_hangji{1,1}(5:4:end));
%
Shike_Line_Str = datestr(Shike_Line,'yyyy-mm-ddTHH-MM-SS-FFFZ');
% 转为数字
% Shike_Line_Num = datenum(Shike_Line_Str);
% intervel=(Shike_Line_Num(3)-Shike_Line_Num(2))*3600*24;   % 时间单位：秒
% num_num=ceil(Tcoh/intervel)+1;
%
% 三维坐标转为数字
% 方位角
Fangwei_Line = str2num(char(data_hangji{1,1}(6:4:end)));
% 俯仰角
Fuyang_Line = str2num(char(data_hangji{1,1}(7:4:end)));
% 斜距
R_Line = str2num(char(data_hangji{1,1}(8:4:end)));
clear data_hangji
%
%%
% 方位角
Az = Fangwei_Line;
% 俯仰角不变
El = Fuyang_Line;
% 依据斜距计算坐标，单位为km
Range = R_Line;

AER=[Az,El,Range];
end