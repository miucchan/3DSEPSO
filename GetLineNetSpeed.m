function [building_Penetrate_Distance,tree_Penetrate_Distance] = GetLineNetSpeed(X1, Y1, Z1, X2, Y2, Z2, img_X0, img_Y0, cell_width, cell_height, dsm_array, building)
%GETLINENETSPEED 

% Obtain the attribute of the communication line passing through the terrain, which is based on the Line2Ras function，
% Determine whether the grids in each communication line have passed through different terrain, and return the penetration distance of each communication line through different terrain
% Among them, features represent the DSM of the features
building_Penetrate_Distance = 0;  % Distance passing through buildings
tree_Penetrate_Distance = 0; % The distance through the forest
not_Building = 0;
%The default starting point is point 1, and it is best to set the service radio station point to point 1, without considering coordinate overlap
if X1 == X2 && Y1 == Y2
    disp('return')
    return
end
[startRow,startCol] = Point2Ras(X1, Y1, img_X0, img_Y0, cell_width, cell_height,dsm_array);
[endRow,endCol] = Point2Ras(X2, Y2, img_X0, img_Y0, cell_width, cell_height,dsm_array);
Current_Z=Z1;

l = sqrt((X1-X2)*(X1-X2)+(Y1-Y2)*(Y1-Y2)+(Z1-Z2)*(Z1-Z2));  % Total length of line segments

if abs((X1-X2)/cell_width) > abs((Y1-Y2)/cell_height)
    k = (Y2-Y1)/(X2-X1);  % slope

    d = abs(endCol-startCol);  % Number of columns with differences

    % The elevation and segment length of each grid tile
    
    divide_Z = (Z1-Z2)/d;
    divide_S = l/d;
    

    if endCol > startCol
        for i =1 : d
            tmpCol = startCol+i;
            tmpRow = floor(startRow-i*k);
            
            Current_Z = Current_Z - divide_Z;
            % If the height of the grid in DSM is greater than or equal to the height of the current grid, it can be considered that the line segment has passed through the part of the grid in DSM
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
    k = (X2-X1)/(Y2-Y1);  

    d = abs(endRow-startRow);  


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
