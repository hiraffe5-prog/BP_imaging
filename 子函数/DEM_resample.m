function [DEM_interp]=DEM_resample(DEM,lon_lat_range,Laxis,Baxis)

% N_interp为插值精度
[Nrow, Ncol] = size(DEM);
lon = linspace(lon_lat_range(1), lon_lat_range(2), Ncol);
lat = linspace(lon_lat_range(3), lon_lat_range(4), Nrow);
[lon1, lat1] = meshgrid(lon, lat);
[Laxis1, Baxis1] = meshgrid(Laxis,Baxis);
DEM_interp = interp2(lon1, lat1, DEM, Laxis1, Baxis1, 'linear');

% lin_X = linspace(1, Ncol, Ncol);
% lin_Y = linspace(1, Nrow, Nrow);
% [X, Y] = meshgrid(lin_X, lin_Y);
% lin_XI = linspace(1, Ncol, (Ncol-1)*N_DEM_interp+1);
% lin_YI = linspace(1, Nrow, (Nrow-1)*N_DEM_interp+1);
% [XI, YI] = meshgrid(lin_XI, lin_YI);
% DEM_interp = interp2(X, Y, DEM, XI, YI, 'linear');

end