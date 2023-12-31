 function [lat, lon] = Ras2Point(i, j, img_X0, img_Y0, cell_width, cell_height)
%POINT2RAS 此处显示有关此函数的摘要
%   此处显示详细说明
lat=img_X0+(j-1)*cell_width+0.5*cell_width;
lon=img_Y0-(i-1)*cell_height-0.5*cell_height;
end

