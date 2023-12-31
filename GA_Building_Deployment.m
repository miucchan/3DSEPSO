clc;
clear;


%% Read Images
tic;
[dsm_array, dsm_refmat] = readgeoraster("DSM.tif"); % Read data

[building_array, building_refmat] = readgeoraster("building.tif");
building_array = building_array(:,:,1);
building_array = bwareaopen(building_array, 20); % Remove small buildings, set to 20 pixels
not_Building = 0; % Pixel value representing non-building areas
radio_Toground_Dis = 1; % Distance from the radio device to the ground

cell_width = dsm_refmat.CellExtentInWorldX;
cell_height = dsm_refmat.CellExtentInWorldY;
img_X = dsm_refmat.XWorldLimits;
img_Y = dsm_refmat.YWorldLimits;

% Read coordinates (top, bottom, left, right)
img_X0 = img_X(1);
img_X1 = img_X(2);
img_Y0 = img_Y(2);
img_Y1 = img_Y(1);

[img_height, img_width] = size(building_array); % Get image dimensions


%% Binary Processing for Building Extraction
disp('Processing buildings...')
[building_pos_array, building_Num] = bwlabel(building_array); % Label connected components in a binary image
conncomp = bwconncomp(building_array); % Find connected components in a binary image and count them
building_info = struct(); % Structure to store information about building grid cells

point_Interval = 20; % Interval for generating test points
% Generate a series of test points within the image for calculating network speed
test_Points = [];
test_Points_Z = 15; % Elevation of test points, relative to the model's elevation

% Record user's planar coordinates to determine the weight of network speed for test points
for i = floor(img_Y0 - point_Interval):-point_Interval:floor(img_Y1 + point_Interval)
    for j = floor(img_X0 + point_Interval):point_Interval:floor(img_X1 - point_Interval)
        % The fourth column represents the maximum network speed the point can achieve,
        % and the fifth column represents whether the point is inside a building (1 for inside, 0 for outside).
        [tmp_row, tmp_col] = Point2Ras(j, i, img_X0, img_Y0, cell_width, cell_height, dsm_array);
        if building_array(tmp_row, tmp_col) == not_Building
            test_Points = [test_Points; [j, i, test_Points_Z, 0, 0]];
        else
            test_Points = [test_Points; [j, i, test_Points_Z, 0, 1]];
        end
    end
end

[test_Points_Num, ~] = size(test_Points); % Number of test points

% Network speed decay coefficients: Y = 41.2771 - 0.0149X_1 - 0.8211X_2 - 0.0003795X_1^2
b0 = 41.2771;
b1 = -0.0149;
b2 = -0.8211;
b3 = -0.0003795;
min_net_speed = 0; % Minimum network speed

for i = 1:building_Num
    building_pos = conncomp.PixelIdxList{i}; % Get all pixel positions representing the building
    [h, ~] = size(building_pos); % Get the total number of pixels representing the building
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
    building_info(i).max_altitude = tmp_max_altitude; % Get the maximum height of the building

    max_altitude_pos_X = mod(tmp_max_altitude_pos, img_height);
    if max_altitude_pos_X == 0
        max_altitude_pos_X = img_height;
    end
    max_altitude_pos_Y = ceil(tmp_max_altitude_pos / img_height);

    tmp_max_altitude_pos = [max_altitude_pos_X, max_altitude_pos_Y]; % Get the pixel position of the building's maximum height
    building_info(i).max_altitude_pos = tmp_max_altitude_pos;
    
    % Get the coordinates to place the radio device at the highest position on the building
    [lat, lon] = Ras2Point(tmp_max_altitude_pos(1), tmp_max_altitude_pos(2), ...
        img_X0, img_Y0, cell_width, cell_height);

    % Calculate the communication network speed between each radio device and test points
    for k = 1:test_Points_Num
        tespoint = test_Points(k,:);
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

    building_info(i).NET = sum(test_Points(:, 4)) / test_Points_Num; % Calculate the network speed obtained when placing the radio device at the highest position of each building
    building_info(i).Fitness = (sum(test_Points(:, 4)) / test_Points_Num) / (1 + std(test_Points(:, 4))); % Calculate the network speed obtained when placing the radio device at the highest position of each building
    test_Points(:, 4) = 0;
end

% Arrange candidate buildings based on the fitness obtained by placing the radio device at the highest height
[b, index] = sort([building_info.Fitness], 'descend'); % Sort network speed values in descending order
building_info = building_info(index);

selectRatio = 0.3; % Selection ratio of buildings
selected_Building_Num = ceil(building_Num * selectRatio); % Calculate the number of buildings selected as radio device placement points

selected_Building_List = building_info(1:selected_Building_Num); % List of selected buildings

% Get the list of selected buildings
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


%% Initialize Parameters Before Iteration
groupNum = 150; % Number of chromosomes
serverNum = 3; % Number of radio devices
gBest = struct('pX', zeros(1, serverNum), 'pY', zeros(1, serverNum), 'NET', 0); % Global best position in history
userNum = 20;
radioUserRatio = serverNum / userNum;

% Initialize chromosome positions
pX = zeros(groupNum, serverNum);
pY = zeros(groupNum, serverNum);
tmp_P_Pos_List = zeros(groupNum, 2);

for i = 1:groupNum
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


%% Calculate Initial pbest and gbest

disp('Initialization begins...')

% Set parameters for velocity update formula
groupEvolutionNum = 50; % Number of evolution iterations

costArray = zeros(1, groupEvolutionNum); % Store historical gbest values

groupMutationP = 0.1; % Chromosome mutation probability

lowFitnessGroup = [];
groupFitness = zeros(1, groupNum);
singleBestFitness=0;

%% Start Iteration
for groupEvolutionTimes = 1:groupEvolutionNum
    groupEvolutionTimes
    for groupIndex = 1:groupNum

        for j = 1:serverNum
            [tmp_row, tmp_col, pZ] = Point2Ras(pX(groupIndex, j), pY(groupIndex, j), ...
                img_X0, img_Y0, cell_width, cell_height, dsm_array);
            if building_array(tmp_row, tmp_col) ~= not_Building
                pZ = pZ + radio_Toground_Dis;
            else
                pZ = test_Points_Z + radio_Toground_Dis;
            end

            serverpoint = [pX(groupIndex, j), pY(groupIndex, j), pZ];

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
        if fitness > singleBestFitness
            singleBestFitness = fitness;
            gBest.pX = pX(groupIndex, :);
            gBest.pY = pY(groupIndex, :);
            best_Points = test_Points;
        end
        groupFitness(groupIndex) = fitness;
        test_Points(:, 4) = 0;
    end

    [b, index] = sort(groupFitness, 'descend');

    costArray(groupEvolutionTimes) = groupFitness(index(1));

    reverseIndex = flip(index); % Reverse ranking

    evolutionIndex = 0;
    evolutionTimes = 0;

    % Population evolution
    % Ratio of selected populations
    selectedRaio = 0.2;
    selectedNum = selectedRaio * groupNum;
    selectedGroupList = zeros(1, selectedNum);

    selectedIndex = 1;
    while selectedIndex <= selectedNum
        randBuildingList = randperm(selected_Building_Num, serverNum);
        selectedGroup = randsample(groupNum, 1, true, groupFitness);
        if ismember(selectedGroup, selectedGroupList)
            continue;
        end
        selectedGroupList(selectedIndex) = selectedGroup;
        selectedIndex = selectedIndex + 1;
    end
    fatherList = randperm(selectedNum);
    crossIndex = 1;
    while crossIndex < selectedNum
        father1 = selectedGroupList(fatherList(crossIndex));
        father2 = selectedGroupList(fatherList(crossIndex + 1));
        crossVec = randi([0, 1], 1, serverNum);
        crossPos = find(crossVec == 1);
        tmp_pX = pX(father1, :);
        pX(father1, crossPos) = pX(father2, crossPos);
        pX(father2, crossPos) = tmp_pX(crossPos);
        crossIndex = crossIndex + 2;
    end
    for mutateIndex = 1:groupNum
        for jj = 1:serverNum
            if rand(1) < groupMutationP
                xShifting = rand(1) * 100 - 50;
                yShifting = rand(1) * 100 - 50;
                pX(mutateIndex, jj) = pX(mutateIndex, jj) + xShifting;
                pY(mutateIndex, jj) = pY(mutateIndex, jj) + yShifting;

                if pX(mutateIndex, jj) <= img_X0
                    pX(mutateIndex, jj) = img_X0 + 0.5;
                end
                if pX(mutateIndex, jj) >= img_X1 - cell_width
                    pX(mutateIndex, jj) = img_X1 - cell_width - 0.5;
                end
                if pY(mutateIndex, jj) <= img_Y1
                    pY(mutateIndex, jj) = img_Y1 + 0.5;
                end
                if pY(mutateIndex, jj) >= img_Y0
                    pY(mutateIndex, jj) = img_Y0 - 0.5;
                end
            end
        end
    end

end
disp(['fitness：' num2str(singleBestFitness)]);

%% Output deployment location
col={'X' 'Y' 'Z'};
pZ=zeros(serverNum,1);
for i= 1:serverNum
    [tmp_row,tmp_col,h]=Point2Ras(gBest.pX(i), gBest.pY(i), img_X0, img_Y0, cell_width, cell_height,dsm_array);
    if building_array(tmp_row,tmp_col) ~= not_Building
        h=h+radio_Toground_Dis;
    else
        h=test_Points_Z+radio_Toground_Dis;
    end
    pZ(i)=h;
end

ANET_Pos_table=table((gBest.pX)',(gBest.pY)',pZ,'VariableNames',col);
writetable(ANET_Pos_table, 'GA_ANET_Pos.csv');
col={'X' 'Y' 'Z' 'NET'};

result_table=table(best_Points(:,1),best_Points(:,2),best_Points(:,3),best_Points(:,4),'VariableNames',col);

writetable(result_table, 'GA_testPointsNet.csv');

