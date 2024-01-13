clc;
clear;

%% Read in image
tic;
[dsm_array, dsm_refmat] = readgeoraster("DSM_mock.tif"); % Read data
[building_array, building_refmat] = readgeoraster("building_mock.tif");
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

%% Initialize parameters before the start of iterations
serverNum = 5; % Number of radio stations
pNum = 150; % Number of particles
userNum = 20;
radioUserRatio = serverNum / userNum;

% Initialize particle positions
pX = zeros(pNum, serverNum);
pY = zeros(pNum, serverNum);
tmp_P_Pos_List = zeros(pNum, 2);

for i = 1:pNum
    for j = 1:serverNum
        pX(i,j)=img_X0 + (img_X1-img_X0) * rand();
        pY(i,j)=img_Y1 + (img_Y0-img_Y1) * rand();
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

        pVx(i, :) = weight * pVx(i, :) + C1(n) * r1 * (pBests(i).pX - pX(i, :)) + C2(n) * r2 * (gBest.pX - pX(i, :)); % Calculate x-axis velocity
        pVy(i, :) = weight * pVy(i, :) + C1(n) * r1 * (pBests(i).pY - pY(i, :)) + C2(n) * r2 * (gBest.pY - pY(i, :)); % Calculate y-axis velocity


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
time=toc;
disp(['Fitness：' num2str(gBest.Fitness)]);
disp(['Netspeed：' num2str(gBest.NET)]);
disp(['Time：' num2str(time)]);
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
