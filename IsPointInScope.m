function [isIn] = IsPointInScope(point_x,point_y,user_x,user_y,scope)

dis=sqrt((point_x-user_x)*(point_x-user_x)+(point_y-user_y)*(point_y-user_y));

if dis<scope
    isIn=true;
else
    isIn=false;
end
end

