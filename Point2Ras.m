 function [i, j, altitude, net] = Point2Ras(X, Y, img_X0, img_Y0, cell_width, cell_height,dsm_array)
%POINT2RAS 此处显示有关此函数的摘要
%   此处显示详细说明
i = floor((img_Y0-Y)/cell_height)+1;
j = floor((X-img_X0)/cell_width)+1;
if i<=0
    fprintf('i小于等于0，i为%d',i)
    Y
end
if j<=0
    
    fprintf('j小于等于0，j为%d',j)
    fprintf('X为%5.5f',X)
end
[h,w]=size(dsm_array);
if i>h
    fprintf('i大于h，i为%d',i)
    Y
end
if j>w
    
    fprintf('j大于w，j为%d',j)
    fprintf('X为%5.5f',X)
end
altitude=dsm_array(i,j);
net=0;
end

