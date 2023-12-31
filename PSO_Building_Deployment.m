clc;
clear;

%% Read in image

[dsm_array, dsm_refmat] = readgeoraster("DSM.tif"); % Read data

[building_array, building_refmat] = readgeoraster("building.tif");
building_array = building_array(:, :, 1);
building_array = bwareaopen(building_array, 20); % Remove small buildings with an area less than 20 pixels
not_Building = 0; % Pixel value representing non-building
radio_Toground_Dis = 1; % Distance from the radio to the ground

cell_width = dsm_refmat.CellExtentInWorldX;
cell_height = dsm_refmat.CellExtentInWorldY;
img_X = dsm_refmat.XWorldLimits;
img_Y = dsm_refmat.YWorldLimits;

% Read coordinates for top, bottom, left, and right
img_X0 = img_X(1);
img_X1 = img_X(2);
img_Y0 = img_Y(2);
img_Y1 = img_Y(1);

[img_height, img_width] = size(building_array);

%% Binary Processing to Extract Buildings
[building_pos_array, building_Num] = bwlabel(building_array); % Label connected components in a 2D binary image
conncomp = bwconncomp(building_array); % Find connected components in the binary image and count them
building_info = struct();

point_Interval = 20; % Interval between generated points
% Generate a series of test points for calculating network speed within the image bounds
test_Points = [];
test_Points_Z = 15; % Elevation of test points, relative to the model elevation

for i = floor(img_Y0 - point_Interval): -point_Interval: floor(img_Y1 + point_Interval)
    for j = floor(img_X0 + point_Interval): point_Interval: floor(img_X1 - point_Interval)
        % The fourth column represents the maximum network speed the point can achieve,
        % and the fifth column indicates whether the point is inside a building (1 for inside, 0 for outside).
        [tmp_row, tmp_col] = Point2Ras(j, i, img_X0, img_Y0, cell_width, cell_height, dsm_array);
        if building_array(tmp_row, tmp_col) == not_Building
            test_Points = [test_Points; [j, i, test_Points_Z, 0, 0]];
        else
            test_Points = [test_Points; [j, i, test_Points_Z, 0, 1]];
        end
    end
end
[test_Points_Num, ~] = size(test_Points); % Number of test points
best_Points = []; % Record data of test points at the global best moment

% Network speed attenuation equation coefficients: Y = 41.2771 - 0.0149X_1 - 0.8211X_2 - 0.0003795X_1^2
b0 = 41.2771;
b1 = -0.0149;
b2 = -0.8211;
b3 = -0.0003795;
min_net_speed = 0; % Minimum network speed rate

for i = 1:building_Num
    building_pos = conncomp.PixelIdxList{i};
    [h, ~] = size(building_pos);
    building_info(i).id = i;
    building_info(i).pos = building_pos;
    tmp_max_altitude = 0;
    tmp_max_altitude_pos = 0;
    for j = 1:h
        if dsm_array(building_pos(j)) > tmp_max_altitude
            tmp_max_altitude = dsm_array(building_pos(j));
            tmp_max_altitude_pos = building_pos(j);
        end
    end
    building_info(i).max_altitude = tmp_max_altitude;

    max_altitude_pos_X = mod(tmp_max_altitude_pos, img_height);
    if max_altitude_pos_X == 0
        max_altitude_pos_X = img_height;
    end
    max_altitude_pos_Y = ceil(tmp_max_altitude_pos / img_height);

    tmp_max_altitude_pos = [max_altitude_pos_X, max_altitude_pos_Y];
    building_info(i).max_altitude_pos = tmp_max_altitude_pos;
    % Get the coordinates to place the radio at the highest position on the building
    [lat, lon] = Ras2Point(tmp_max_altitude_pos(1), tmp_max_altitude_pos(2), ...
        img_X0, img_Y0, cell_width, cell_height);

    for k = 1:test_Points_Num
        tespoint = test_Points(k, :);
        % Calculate penetration distance

        [building_dis, tree_dis] = GetLineNetSpeed(lat, lon, tmp_max_altitude + radio_Toground_Dis, ...
            tespoint(1), tespoint(2), tespoint(3), img_X0, img_Y0, cell_width, cell_height, dsm_array, building_array);
        net_speed = floor(b0 + b1 * tree_dis + b2 * building_dis + b3 * tree_dis * tree_dis); % Calculate network speed
        isLastPointPredicted = 0;
        last_building_dis = building_dis;
        if net_speed > min_net_speed
            test_Points(k, 4) = net_speed;
        end
    end

    building_info(i).NET = sum(test_Points(:, 4)) / test_Points_Num; % Calculate the network speed obtained when placing the radio at the highest position on each building
    building_info(i).Fitness = (sum(test_Points(:, 4)) / test_Points_Num) / (1 + std(test_Points(:, 4))); % Calculate the network speed obtained when placing the radio at the highest position on each building
    test_Points(:, 4) = 0;

end

% Arrange candidate buildings based on fitness obtained by placing the radio at the highest altitude
[b, index] = sort([building_info.Fitness], 'descend'); % Sort network speed values in descending order
building_info = building_info(index);

selectRatio = 0.3; % Select a proportion of buildings
selected_Building_Num = ceil(building_Num * selectRatio); % Calculate the number of selected buildings as radio placement points

selected_Building_List = building_info(1:selected_Building_Num); % List of selected buildings
selected_Building_Weight = zeros(1, selected_Building_Num);
for i = 1:selected_Building_Num
    pos_List = selected_Building_List(i).pos;
    tmp_list = [mod(pos_List, img_height), ceil(pos_List ./ img_height)];
    tmp_list(tmp_list == 0) = img_height;
    selected_Building_List(i).pos = tmp_list;
    selected_Building_List(i).orderNum = i;
    selected_Building_Weight(i) = selected_Building_List(i).Fitness;
end

disp('Building processing completed')


%% Initialize parameters before the start of iterations
serverNum = 3; % Number of radio stations
pNum = 150; % Number of particles
userNum = 20;
radioUserRatio = serverNum / userNum;

% Initialize particle positions
pX = zeros(pNum, serverNum);
pY = zeros(pNum, serverNum);
tmp_P_Pos_List = zeros(pNum, 2);

for i = 1:pNum
    randBuildingList = randperm(selected_Building_Num, serverNum);
    for j = 1:serverNum
        building_Index = randBuildingList(j);
        tmp_Building_Pos = selected_Building_List(building_Index).pos;
        [building_Area, ~] = size(tmp_Building_Pos);
        rand_Pos = randperm(building_Area, 1);
        tmp_P_Pos = tmp_Building_Pos(rand_Pos, :);
        tmp_P_Pos_List(i, 1) = tmp_P_Pos(1);
        tmp_P_Pos_List(i, 2) = tmp_P_Pos(2);
        [lat, lon] = Ras2Point(tmp_P_Pos(1), tmp_P_Pos(2), img_X0, img_Y0, cell_width, cell_height);
        pX(i, j) = lat;
        pY(i, j) = lon;
    end
end

pVx = zeros(pNum, serverNum); % Particle X velocities, initialized to 0
pVy = zeros(pNum, serverNum); % Particle Y velocities, initialized to 0
pBests = struct('pX', zeros(1, serverNum), 'pY', zeros(1, serverNum), 'NET', 0); % Historical best position for each particle
gBest = struct('pX', zeros(1, serverNum), 'pY', zeros(1, serverNum), 'NET', 0); % Global historical best position

Building_Distance = zeros(serverNum, test_Points_Num);
Tree_Distance = zeros(serverNum, test_Points_Num);


%% Calculate initial pbest and gbest

pp = 1;
num = 0;
max_fitness = 0;

for i = 1:pNum

    for j = 1:serverNum

        [tmp_row, tmp_col, pZ] = Point2Ras(pX(i, j), pY(i, j), img_X0, img_Y0, cell_width, cell_height, dsm_array);

        if building_array(tmp_row, tmp_col) ~= not_Building
            pZ = pZ + radio_Toground_Dis;
        else
            pZ = test_Points_Z + radio_Toground_Dis;
        end

        serverpoint = [pX(i, j), pY(i, j), pZ];

        for k = 1:test_Points_Num
            tespoint = test_Points(k, :);

            % Calculate penetration distance
            [building_dis, tree_dis] = GetLineNetSpeed(serverpoint(1), serverpoint(2), serverpoint(3), ...
                tespoint(1), tespoint(2), tespoint(3), img_X0, img_Y0, cell_width, cell_height, dsm_array, building_array);

            net_speed = floor(b0 + b1 * tree_dis + b2 * building_dis + b3 * tree_dis * tree_dis); % Calculate network speed
            isLastPointPredicted = 0;
            last_building_dis = building_dis;

            if net_speed > min_net_speed && radioUserRatio * net_speed > test_Points(k, 4)
                test_Points(k, 4) = radioUserRatio * net_speed;
            end

        end
    end

    % Calculate the weighted average of the maximum and average network speed for each point
    sum_speed = sum(test_Points(:, 4)) / test_Points_Num;
    fitness = sum_speed / (1 + std(test_Points(:, 4)));

    % Update pBests
    pBests(i).pX = pX(i, :);
    pBests(i).pY = pY(i, :);
    pBests(i).NET = sum_speed;
    pBests(i).Fitness = fitness;

    if fitness > max_fitness
        max_fitness = fitness;
        pp = i;
        best_Points = test_Points;
    end

    % Update gBest
    gBest.pX = pX(pp, :);
    gBest.pY = pY(pp, :);
    gBest.NET = pBests(pp).NET;
    gBest.Fitness = pBests(pp).Fitness;
    test_Points(:, 4) = 0;

end

%% Set parameters for velocity update formula

particleUpdateNum = 50;
x = 1:(10-1)/particleUpdateNum:10;
C1 = (-x.^5/200000 + 1); % C1 learning rate
C2 = (-exp(x.^-0.5)/1.5 + 2); % C2 learning rate
crossP = 0.1; % Crossover probability

% mutationP = 0.05:(0.005-0.05)/particleUpdateNum:0.005; % Mutation probability

iter = 0:2.3/particleUpdateNum:2.3;

%% Start the iteration

for n = 1:particleUpdateNum
    n
    weight = 0.4 + 0.5 * exp(-iter(n) * iter(n));
    for i = 1:pNum

        r1 = rand(1);
        crossRand = rand(1);
        r2 = rand(1);

        if crossRand > crossP
            pVx(i, :) = weight * pVx(i, :) + C1(n) * r1 * (pBests(i).pX - pX(i, :)) + C2(n) * r2 * (gBest.pX - pX(i, :)); % Calculate x-axis velocity
            pVy(i, :) = weight * pVy(i, :) + C1(n) * r1 * (pBests(i).pY - pY(i, :)) + C2(n) * r2 * (gBest.pY - pY(i, :)); % Calculate y-axis velocity

        else
            pVx(i, :) = r1 * (2 * gBest.pX - pBests(i).pX - pX(i, :)); % Crossover occurred, calculate x-axis velocity
            pVy(i, :) = r1 * (2 * gBest.pY - pBests(i).pY - pY(i, :)); % Calculate y-axis velocity
        end

        %         mutationRand=rand(1);
        %         if mutationRand<mutationP(n)
        %             tmp=ceil(rand(1)*serverNum);
        %             pX(i,tmp)=pX(i,tmp)+rand(1)*2-1;
        %             pY(i,tmp)=pY(i,tmp)+rand(1)*2-1;
        %
        %         end
        % Update particle position
        tmp_px = pX(i, :) + pVx(i, :);
        tmp_py = pY(i, :) + pVy(i, :);
        for pxnum = 1:serverNum
            if tmp_px(pxnum) < img_X0 + cell_width * 2
                tmp_px(pxnum) = img_X0 + cell_width * 2;
            end
            if tmp_px(pxnum) > img_X1 - cell_width * 2
                tmp_px(pxnum) = img_X1 - cell_width * 2;
            end
            if tmp_py(pxnum) < img_Y1 + cell_height * 2
                tmp_py(pxnum) = img_Y1 + cell_height * 2;
            end
            if tmp_py(pxnum) > img_Y0 - cell_height * 2
                tmp_py(pxnum) = img_Y0 - cell_height * 2;
            end
        end

        pX(i, :) = tmp_px;
        pY(i, :) = tmp_py;

        for j = 1:serverNum
            [tmp_row, tmp_col, pZ] = Point2Ras(pX(i, j), pY(i, j), img_X0, img_Y0, cell_width, cell_height, dsm_array);
            if building_array(tmp_row, tmp_col) ~= not_Building
                pZ = pZ + radio_Toground_Dis;
            else
                pZ = test_Points_Z + radio_Toground_Dis;
            end
            serverpoint = [pX(i, j), pY(i, j), pZ];

            for k = 1:test_Points_Num
                tespoint = test_Points(k, :);
                % Calculate penetration distance

                [building_dis, tree_dis] = GetLineNetSpeed(serverpoint(1), serverpoint(2), serverpoint(3), ...
                    tespoint(1), tespoint(2), tespoint(3), img_X0, img_Y0, cell_width, cell_height, dsm_array, building_array);
                net_speed = floor(b0 + b1 * tree_dis + b2 * building_dis + b3 * tree_dis * tree_dis); % Calculate net speed
                isLastPointPredicted = 0;
                last_building_dis = building_dis;
                if net_speed > min_net_speed && radioUserRatio * net_speed > test_Points(k, 4)
                    test_Points(k, 4) = radioUserRatio * net_speed;
                end

            end
        end
        sum_speed = sum(test_Points(:, 4)) / test_Points_Num;
        fitness = sum_speed / (1 + std(test_Points(:, 4)));

        % Update pbest
        if fitness > pBests(i).Fitness
            pBests(i).NET = sum_speed;
            pBests(i).Fitness = fitness;
            pBests(i).pX = pX(i, :);
            pBests(i).pY = pY(i, :);
        end
        % Update gbest
        if fitness > gBest.Fitness
            gBest.NET = sum_speed;
            gBest.Fitness = fitness;
            gBest.pX = pX(i, :);
            gBest.pY = pY(i, :);
        end
        test_Points(:, 4) = 0;
    end



end
disp(['fitness：' num2str(gBest.Fitness)]);
%% Output deployment location
col = {'X' 'Y' 'Z'};
pZ = zeros(serverNum, 1);
for i = 1:serverNum
    [tmp_row, tmp_col, h] = Point2Ras(gBest.pX(i), gBest.pY(i), img_X0, img_Y0, cell_width, cell_height, dsm_array);
    if building_array(tmp_row, tmp_col) ~= not_Building
        h = h + radio_Toground_Dis;
    else
        h = test_Points_Z + radio_Toground_Dis;
    end
    pZ(i) = h;
end

ANET_Pos_table = table((gBest.pX)', (gBest.pY)', pZ, 'VariableNames', col);
writetable(ANET_Pos_table, 'PSO_ANET_Pos.csv');
col = {'X' 'Y' 'Z' 'NET'};
result_table = table(best_Points(:, 1), best_Points(:, 2), best_Points(:, 3), best_Points(:, 4), 'VariableNames', col);

writetable(result_table, 'PSO_testPointsNet.csv');
