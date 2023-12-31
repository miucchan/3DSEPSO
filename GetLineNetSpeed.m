function [building_Penetrate_Distance,tree_Penetrate_Distance] = GetLineNetSpeed(X1, Y1, Z1, X2, Y2, Z2, img_X0, img_Y0, cell_width, cell_height, dsm_array, building)
%GETLINENETSPEED 此处显示有关此函数的摘要
%   此处显示详细说明
% 获得通信线穿过地物的属性，该函数在Line2Ras函数的基础上，
% 判断每根通信线中的栅格是否通过了不同地物，并返回每根通信线通过不同地物的穿透距离
% 其中features为代表地物的DSM
building_Penetrate_Distance = 0;  % 穿过建筑物的距离
tree_Penetrate_Distance = 0; % 穿过树林的距离
not_Building = 0;
%起始点默认为点1，最好将服务电台点设置为点1,不考虑坐标重合情况
if X1 == X2 && Y1 == Y2
    disp('return')
    return
end
[startRow,startCol] = Point2Ras(X1, Y1, img_X0, img_Y0, cell_width, cell_height,dsm_array);
[endRow,endCol] = Point2Ras(X2, Y2, img_X0, img_Y0, cell_width, cell_height,dsm_array);
Current_Z=Z1;

l = sqrt((X1-X2)*(X1-X2)+(Y1-Y2)*(Y1-Y2)+(Z1-Z2)*(Z1-Z2));  % 线段总长度

if abs((X1-X2)/cell_width) > abs((Y1-Y2)/cell_height)
    k = (Y2-Y1)/(X2-X1);  % 斜率

    d = abs(endCol-startCol);  % 相差的列数

    % 每个栅格平摊的高程和线段长度
    
    divide_Z = (Z1-Z2)/d;
    divide_S = l/d;
    

    if endCol > startCol
        for i =1 : d
            tmpCol = startCol+i;
            tmpRow = floor(startRow-i*k);
            
            Current_Z = Current_Z - divide_Z;
            % 如果DSM的栅格的高度大于或等于当前栅格的高度，可以认为线段通过了DSM该栅格的部分
            if dsm_array(tmpRow,tmpCol) >= Current_Z
                if building(tmpRow,tmpCol)~=not_Building
                    building_Penetrate_Distance = building_Penetrate_Distance+ divide_S;
                else
                    tree_Penetrate_Distance = tree_Penetrate_Distance+divide_S;
                end
            end
        end
    else
        for i =1 : d
            tmpCol = startCol-i;
            tmpRow = floor(startRow+i*k);

            Current_Z = Current_Z - divide_Z;

            if dsm_array(tmpRow,tmpCol) >= Current_Z
                if building(tmpRow,tmpCol)~=not_Building
                    building_Penetrate_Distance = building_Penetrate_Distance+ divide_S;
                else
                    tree_Penetrate_Distance = tree_Penetrate_Distance+divide_S;
                end
            end
        end
    end
else
    k = (X2-X1)/(Y2-Y1);  % 斜率

    d = abs(endRow-startRow);  % 相差的行数

    % 每个栅格平摊的高程和线段长度
    divide_Z = (Z1-Z2)/d;
    divide_S = l/d;
    
    if endRow > startRow
        for i = 1 : d
            tmpRow = startRow+i;
            tmpCol = floor(startCol-i*k);

            Current_Z = Current_Z - divide_Z;

            if dsm_array(tmpRow,tmpCol) >= Current_Z
                if building(tmpRow,tmpCol)~=not_Building
                    building_Penetrate_Distance = building_Penetrate_Distance+divide_S;
                else
                    tree_Penetrate_Distance = tree_Penetrate_Distance+divide_S;
                end
            end
        end
    else
        for i = 1 : d
            tmpRow = startRow-i;
            tmpCol = floor(startCol+i*k);

            Current_Z = Current_Z - divide_Z;

            if dsm_array(tmpRow,tmpCol) >= Current_Z
                if building(tmpRow,tmpCol)~=not_Building
                    building_Penetrate_Distance =building_Penetrate_Distance+ divide_S;
                else
                    tree_Penetrate_Distance =tree_Penetrate_Distance+ divide_S;

                end
            end
        end
    end
end
