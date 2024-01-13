 function [i, j, altitude, net] = Point2Ras(X, Y, img_X0, img_Y0, cell_width, cell_height,dsm_array)

i = floor((img_Y0-Y)/cell_height)+1;
j = floor((X-img_X0)/cell_width)+1;

altitude=dsm_array(i,j);
net=0;
end

