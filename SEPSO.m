clc;
clear;

%% Read in image
tic;
%read data
[dsm_array,dsm_refmat] = readgeoraster("DSM_mock.tif");
[building_array,building_refmat] = readgeoraster("building_mock.tif");

building_array=building_array(:,:,1);
building_array=bwareaopen(building_array,20);%Excluding buildings with smaller areas

not_Building=0;%Pixel values not represented by buildings
radio_Toground_Dis=1; %The distance between the radio station and the ground

cell_width=dsm_refmat.CellExtentInWorldX;
cell_height=dsm_refmat.CellExtentInWorldY;
img_X=dsm_refmat.XWorldLimits;
img_Y=dsm_refmat.YWorldLimits;

%Read the upper, lower, left, and right coordinates
img_X0=img_X(1);
img_X1=img_X(2);
img_Y0=img_Y(2);
img_Y1=img_Y(1);

[img_height,img_width]=size(building_array);%Obtain image size

%% Binary processing for extracting buildings
disp('Processing buildings...')
[building_pos_array,building_Num]=bwlabel(building_array); % Annotating Connected Components in 2D Binary Images
conncomp=bwconncomp(building_array); % Search for connected components in binary images and count them
building_info=struct(); %Used to store information about building grids

point_Interval=20; %The interval between generated points
%Generate a series of test points within the image range to calculate network speed
test_Points=[];
test_Points_Z=15; %The elevation of the test point, which is relative to the model's elevation

for i = floor(img_Y0-point_Interval):-point_Interval:floor(img_Y1+point_Interval)
    for j = floor(img_X0+point_Interval):point_Interval:floor(img_X1-point_Interval)
        %The fourth column represents the maximum internet speed that a point can obtain, the fifth column represents whether the point is located inside the building, 1 represents yes, 0 represents no

        [tmp_row,tmp_col]=Point2Ras(j, i,img_X0, img_Y0, cell_width, cell_height, dsm_array);

        if building_array(tmp_row,tmp_col)==not_Building
            test_Points=[test_Points;[j,i,test_Points_Z,0,0]];
        else
            test_Points=[test_Points;[j,i,test_Points_Z,0,1]];
        end


    end
end

[test_Points_Num,~]=size(test_Points);

% NSAM：Y=41.2771-0.0149X_1-0.8211X_2-0.0003795X_1^2
b0=41.2771;
b1=-0.0149;
b2=-0.8211;
b3=-0.0003795;
min_net_speed=0; %Minimum network speed

for i=1:building_Num
    building_pos=conncomp.PixelIdxList{i}; %Obtain all pixel positions representing the building
    [h,~]=size(building_pos); %Obtain the total sum of all pixels representing the building
    building_info(i).id=i;
    building_info(i).pos=building_pos;
    tmp_max_altitude=0;
    tmp_max_altitude_pos=0;
    for j=1:h
        if dsm_array(building_pos(j))>tmp_max_altitude
            tmp_max_altitude=dsm_array(building_pos(j));
            tmp_max_altitude_pos=building_pos(j);
        end
    end
    building_info(i).max_altitude=tmp_max_altitude; %Obtain the maximum height of the building

    max_altitude_pos_X=mod(tmp_max_altitude_pos,img_height);
    if max_altitude_pos_X==0
        max_altitude_pos_X=img_height;
    end
    max_altitude_pos_Y=ceil(tmp_max_altitude_pos/img_height);

    tmp_max_altitude_pos=[max_altitude_pos_X max_altitude_pos_Y];%Obtain the pixel position of the maximum height of the building
    building_info(i).max_altitude_pos=tmp_max_altitude_pos;
    %Obtain the coordinates for placing the radio station at the highest position on the building
    [lat, lon]=Ras2Point(tmp_max_altitude_pos(1), tmp_max_altitude_pos(2), ...
        img_X0, img_Y0, cell_width, cell_height);

    %Calculate the communication network speed between each radio station and the test point
    for k =1:test_Points_Num
        tespoint=test_Points(k,:);
        %Calculate penetration distance
        [building_dis,tree_dis]=GetLineNetSpeed(lat,lon,tmp_max_altitude+radio_Toground_Dis, ...
            tespoint(1),tespoint(2),tespoint(3),img_X0,img_Y0,cell_width,cell_height,dsm_array,building_array);
        net_speed=floor(b0+b1*tree_dis+b2*building_dis+b3*tree_dis*tree_dis);%Calculate network speed
        isLastPointPredicted=0;
        last_building_dis=building_dis;
        if net_speed>min_net_speed
            test_Points(k,4)=net_speed;
        end
    end
    building_info(i).NET=sum(test_Points(:,4))/test_Points_Num; %Calculate the network speed obtained when placing a radio station at the highest position of each building
    building_info(i).Fitness=(sum(test_Points(:,4))/test_Points_Num)/(1+std(test_Points(:,4))); %Calculate the network speed obtained when placing a radio station at the highest position of each building
    test_Points(:,4)=0;

end

% Arrange candidate buildings based on the fitness values obtained by placing the radio station at the highest height
[b,index]=sort([building_info.Fitness],'descend'); %Sort the network speed values in descending order
building_info=building_info(index);
selectRatio=0.5;%Select the proportion of buildings
selected_Building_Num=ceil(building_Num*selectRatio);%Calculate the number of buildings selected as radio station placement points
selected_Building_List=building_info(1:selected_Building_Num); %List of selected buildings

%Obtain a list of selected buildings
selected_Building_Weight=zeros(1,selected_Building_Num);
for i=1:selected_Building_Num
    pos_List=selected_Building_List(i).pos;
    tmp_list=[mod(pos_List,img_height) ceil(pos_List./img_height)];
    tmp_list(tmp_list==0)=img_height;
    selected_Building_List(i).pos=tmp_list;
    selected_Building_List(i).orderNum=i;
    selected_Building_Weight(i)=selected_Building_List(i).Fitness;
end
disp('Building processing completed')

%% Parameters required before initialization iteration begins
groupNum = 30; % Number of chromosomes
serverNum = 5; % Number of base stations
pNum = 5; % Number of particles, 5 particles are optimal
userNum = 20;
radioUserRatio = serverNum / userNum;

% Initialize particle positions
pX = zeros(groupNum, pNum, serverNum);
pY = zeros(groupNum, pNum, serverNum);
randBuildingList = struct();
selected_Building_Num_List = 1:selected_Building_Num;
% For each population, first select the buildings to be placed, generate particles, and then randomly select locations to place base stations in each building
for groupIndex = 1:groupNum
    tmp_RandBuildingList = zeros(1, serverNum);
    serverIndex = 1;
    % Obtain the building numbers for each population
    while serverIndex <= serverNum
        randBuilding = randsample(selected_Building_Num_List, 1, true, selected_Building_Weight);
        if ismember(randBuilding, tmp_RandBuildingList)
            continue;
        end
        tmp_RandBuildingList(serverIndex) = randBuilding;
        serverIndex = serverIndex + 1;
    end

    randBuildingList(groupIndex).selectedBuilding = tmp_RandBuildingList;
    for i = 1:pNum
        for j = 1:serverNum
            building_Index = tmp_RandBuildingList(j);
            tmp_Building_Pos = selected_Building_List(building_Index).pos;
            [building_Area, ~] = size(tmp_Building_Pos);
            rand_Pos = randperm(building_Area, 1);
            tmp_P_Pos = tmp_Building_Pos(rand_Pos, :);

            [lat, lon] = Ras2Point(tmp_P_Pos(1), tmp_P_Pos(2), img_X0, img_Y0, cell_width, cell_height);
            pX(groupIndex, i, j) = lat;
            pY(groupIndex, i, j) = lon;
        end
    end
end

pVx = zeros(groupNum, pNum, serverNum); % Particle X velocity, initialized to 0
pVy = zeros(groupNum, pNum, serverNum); % Particle Y velocity, initialized to 0
pBests = struct('pX', zeros(1, serverNum), 'pY', zeros(1, serverNum), 'NET', 0, 'Fitness', 0); % Historical best position for each particle in the population
gBests = struct('pX', zeros(1, serverNum), 'pY', zeros(1, serverNum), 'NET', 0, 'Fitness', 0); % Historical best position for each population
best_Points = struct(); % Record data for the globally best moments

%% Calculate initial pbest and gbest
disp('Initialization started...')
for groupIndex = 1:groupNum
    max_fitness = 0; % Maximum average network speed
    pp = 1;
    for i = 1:pNum
        for j = 1:serverNum
            [tmp_row, tmp_col, pZ] = Point2Ras(pX(groupIndex, i, j), pY(groupIndex, i, j), ...
                img_X0, img_Y0, cell_width, cell_height, dsm_array);
            if building_array(tmp_row, tmp_col) ~= not_Building
                pZ = pZ + radio_Toground_Dis;
            else
                pZ = test_Points_Z + radio_Toground_Dis;
            end
            serverpoint = [pX(groupIndex, i, j), pY(groupIndex, i, j), pZ];
            for k = 1:test_Points_Num
                tespoint = test_Points(k, :);
                % Calculate penetration distance
                [building_dis, tree_dis] = GetLineNetSpeed(serverpoint(1), serverpoint(2), serverpoint(3), ...
                    tespoint(1), tespoint(2), tespoint(3), img_X0, img_Y0, cell_width, cell_height, dsm_array, building_array);
                net_speed = floor(b0 + b1 * tree_dis + b2 * building_dis + b3 * tree_dis * tree_dis); % Calculate network speed
                if net_speed > min_net_speed && radioUserRatio * net_speed > test_Points(k, 4)
                    test_Points(k, 4) = radioUserRatio * net_speed;
                end
            end
        end

        % Calculate the weighted average of the maximum and average network speeds for each point
        sum_speed = sum(test_Points(:, 4)) / test_Points_Num;
        fitness = sum_speed / (1 + std(test_Points(:, 4)));
        % Update pBests
        pBests(groupIndex, i).pX = pX(groupIndex, i, :);
        pBests(groupIndex, i).pY = pY(groupIndex, i, :);
        pBests(groupIndex, i).NET = sum_speed;
        pBests(groupIndex, i).Fitness = fitness;
        if fitness > max_fitness
            max_fitness = fitness;
            pp = i;
            best_Points(groupIndex).points = test_Points;
        end
        % Update gBest
        gBests(groupIndex).pX = pX(groupIndex, pp, :);
        gBests(groupIndex).pY = pY(groupIndex, pp, :);
        gBests(groupIndex).NET = pBests(groupIndex, pp).NET;
        gBests(groupIndex).Fitness = pBests(groupIndex, pp).Fitness;
        gBests(groupIndex).id = groupIndex;
        test_Points(:, 4) = 0;
    end
end
disp('Initialization completed')

%% Set parameters for velocity update formula
groupEvolutionNum = 50; % Number of evolution iterations
particleUpdateNum = 1; % Number of particle updates per evolution iteration

surviveGroupChangeRatio = 0.5; % This value represents the ratio of surviving populations that undergo mutation or crossover in each evolution iteration
groupDie = 0.1; % Probability of population extinction
groupDieList = ones(1, groupNum); % Record populations that have already become extinct
groupMutationP = 0.75; % Population mutation probability
groupCrossP = 1 - groupMutationP; % Population crossover probability
x = 1:(10-1)/groupEvolutionNum:10;
C1 = (-x.^5/200000 + 1); % C1 learning rate
C2 = (-exp(x.^-0.5)/1.5 + 2); % C2 learning rate
iter = 0:2.3/groupEvolutionNum:2.3;
lowFitnessGroup = [];


%% Start iteration
for groupEvolutionTimes=1:groupEvolutionNum
    groupEvolutionTimes
    for groupIndex=1:groupNum
        if groupDieList(groupIndex)==0
            continue;
        end
        for n=1:particleUpdateNum
            weight=0.4+0.5*exp(-iter(groupEvolutionTimes)*iter(groupEvolutionTimes));
            for i=1:pNum
                r1=rand(1);
                crossRand=rand(1);
                r2=rand(1);

                pVx(groupIndex,i,:)=C1(groupEvolutionTimes)*r1*(pBests(groupIndex,i).pX-pX(groupIndex,i,:))+C2(groupEvolutionTimes)*r2*(gBests(groupIndex).pX-pX(groupIndex,i,:));%Calculate the velocity in the x direction
                pVy(groupIndex,i,:)=C1(groupEvolutionTimes)*r1*(pBests(groupIndex,i).pY-pY(groupIndex,i,:))+C2(groupEvolutionTimes)*r2*(gBests(groupIndex).pY-pY(groupIndex,i,:));%Calculate the velocity in the y direction

                %Update particle positions
                dTimes=zeros(1,serverNum);%The number of times a particle is divided by 2

                %Prevent crossing image and building boundaries. If particles fly out of the building, reduce the update speed by half and recalculate
                while true
                    tmp_px=pX(groupIndex,i,:)+pVx(groupIndex,i,:);
                    tmp_py=pY(groupIndex,i,:)+pVy(groupIndex,i,:);
                    isOver=true;

                    for pxnum=1:serverNum

                        if tmp_px(pxnum)<img_X0
                            tmp_px(pxnum)=img_X0;
                        end
                        if tmp_px(pxnum)>img_X1
                            tmp_px(pxnum)=img_X1;
                        end
                        if tmp_py(pxnum)<img_Y1
                            tmp_py(pxnum)=img_Y1;
                        end
                        if tmp_py(pxnum)>img_Y0
                            tmp_py(pxnum)=img_Y0;
                        end
                        [tmp_row,tmp_col]=Point2Ras(tmp_px(pxnum), tmp_py(pxnum), ...
                            img_X0, img_Y0, cell_width, cell_height, dsm_array);
                        if building_array(tmp_row,tmp_col)==not_Building
                            pVx(groupIndex,i,pxnum)=pVx(groupIndex,i,pxnum)/2;
                            pVy(groupIndex,i,pxnum)=pVy(groupIndex,i,pxnum)/2;
                            isOver=false;
                            dTimes(pxnum)=dTimes(pxnum)+1;
                        end
                        if dTimes(pxnum)>5
                            pVx(groupIndex,i,pxnum)=0;
                            pVy(groupIndex,i,pxnum)=0;
                        end

                    end
                    if isOver
                        break
                    end
                end

                %Assign new positions to particles
                pX(groupIndex,i,:)=tmp_px;
                pY(groupIndex,i,:)=tmp_py;

                for j =1:serverNum
                    [tmp_row,tmp_col,pZ]=Point2Ras(pX(groupIndex,i,j), pY(groupIndex,i,j), img_X0, img_Y0, cell_width, cell_height, dsm_array);
                    if building_array(tmp_row,tmp_col) ~= not_Building
                        pZ=pZ+radio_Toground_Dis;
                    else
                        pZ=test_Points_Z+radio_Toground_Dis;
                    end
                    serverpoint=[pX(groupIndex,i,j),pY(groupIndex,i,j),pZ];

                    for k =1:test_Points_Num
                        tespoint=test_Points(k,:);
                        %Calculate penetration distance
                        [building_dis,tree_dis]=GetLineNetSpeed(serverpoint(1),serverpoint(2),serverpoint(3), ...
                            tespoint(1),tespoint(2),tespoint(3),img_X0,img_Y0,cell_width,cell_height,dsm_array,building_array);
                        net_speed=floor(b0+b1*tree_dis+b2*building_dis+b3*tree_dis*tree_dis);%Calculate network speed
                        isLastPointPredicted=0;
                        last_building_dis=building_dis;
                        if net_speed>min_net_speed && radioUserRatio*net_speed>test_Points(k,4)
                            test_Points(k,4)=radioUserRatio*net_speed;
                        end

                    end
                end

                sum_speed=sum(test_Points(:,4))/test_Points_Num;
                fitness=sum_speed/(1+std(test_Points(:,4)));

                %Update pbest
                if fitness>pBests(groupIndex,i).Fitness
                    pBests(groupIndex,i).NET=sum_speed;
                    pBests(groupIndex,i).Fitness=fitness;
                    pBests(groupIndex,i).pX=pX(groupIndex,i,:);
                    pBests(groupIndex,i).pY=pY(groupIndex,i,:);
                end
                %Update gbest
                if fitness>gBests(groupIndex).Fitness
                    gBests(groupIndex).NET=sum_speed;
                    gBests(groupIndex).Fitness=fitness;
                    gBests(groupIndex).pX=pX(groupIndex,i,:);
                    gBests(groupIndex).pY=pY(groupIndex,i,:);
                    best_Points(groupIndex).points=test_Points;
                end
                test_Points(:,4)=0;
            end

        end
    end
    surviveGroupChangeNum=floor(sum(groupDieList==1)*surviveGroupChangeRatio);
    [b,index]=sort([gBests.Fitness],'descend');
    ggBest=gBests(index(1));
    descenSsortedgBests=gBests(fliplr(index));%Reverse ranking
    evolutionIndex=0;
    evolutionTimes=0;

    %Population Evolution
    while evolutionTimes<surviveGroupChangeNum
        evolutionIndex=evolutionIndex+1;
        currentGroupId=descenSsortedgBests(evolutionIndex).id;
        if groupDieList(currentGroupId)==0
            continue;
        end
        %Generate random numbers to determine whether the current population is dead
        if rand(1)<groupDie
            groupDieList(currentGroupId)=0;
            evolutionTimes=evolutionTimes+1;
            continue
        end

        %Add low fitness values to the list for storage
        cur_LowFitnessGroup=randBuildingList(currentGroupId).selectedBuilding;
        lowFitnessGroup=[lowFitnessGroup;cur_LowFitnessGroup];

        if rand(1)<groupMutationP
            mutateFailTimes=0;%Allow up to three failed mutations
            while mutateFailTimes<3
                isMutated=rand(1,serverNum);%Values less than 0.5 undergo variation, while others remain unchanged
                for serverIndex=1:serverNum
                    if isMutated(serverIndex)<0.5
                        while true
                            randBuilding=randsample(selected_Building_Num_List,1,true,selected_Building_Weight);
                            if ismember(randBuilding,cur_LowFitnessGroup)
                                continue;
                            end
                            cur_LowFitnessGroup(serverIndex)=randBuilding;
                            break
                        end
                    end
                end
                if ismember(cur_LowFitnessGroup,lowFitnessGroup,'rows')
                    mutateFailTimes=mutateFailTimes+1;
                    continue;
                end
                randBuildingList(currentGroupId).selectedBuilding=cur_LowFitnessGroup;
                break
            end

            %Reinitialize the population
            for i=1:pNum
                for j=1:serverNum
                    building_Index=cur_LowFitnessGroup(j);
                    tmp_Building_Pos=selected_Building_List(building_Index).pos;
                    [building_Area,~]=size(tmp_Building_Pos);
                    rand_Pos=randperm(building_Area,1);
                    tmp_P_Pos=tmp_Building_Pos(rand_Pos,:);

                    [lat, lon] = Ras2Point(tmp_P_Pos(1), tmp_P_Pos(2), img_X0, img_Y0, cell_width, cell_height);
                    pX(currentGroupId,i,j)=lat;
                    pY(currentGroupId,i,j)=lon;
                end
            end
            pVx(currentGroupId,:,:)=0;
            pVy(currentGroupId,:,:)=0;
            max_fitness=0;
            pp=1;
            for i=1:pNum
                for j =1:serverNum
                    [tmp_row,tmp_col,pZ]=Point2Ras(pX(currentGroupId,i,j), pY(currentGroupId,i,j), ...
                        img_X0, img_Y0, cell_width, cell_height, dsm_array);
                    if building_array(tmp_row,tmp_col) ~= not_Building
                        pZ=pZ+radio_Toground_Dis;
                    else
                        pZ=test_Points_Z+radio_Toground_Dis;
                    end

                    serverpoint=[pX(currentGroupId,i,j),pY(currentGroupId,i,j),pZ];
                    for k =1:test_Points_Num
                        tespoint=test_Points(k,:);

                        [building_dis,tree_dis]=GetLineNetSpeed(serverpoint(1),serverpoint(2),serverpoint(3), ...
                            tespoint(1),tespoint(2),tespoint(3),img_X0,img_Y0,cell_width,cell_height,dsm_array,building_array);
                        net_speed=floor(b0+b1*tree_dis+b2*building_dis+b3*tree_dis*tree_dis);
                        isLastPointPredicted=0;
                        last_building_dis=building_dis;
                        if net_speed>min_net_speed && radioUserRatio*net_speed>test_Points(k,4)
                            test_Points(k,4)=radioUserRatio*net_speed;
                        end

                    end
                end

                sum_speed=sum(test_Points(:,4))/test_Points_Num;
                fitness=sum_speed/(1+std(test_Points(:,4)));

                pBests(currentGroupId,i).pX=pX(currentGroupId,i,:);
                pBests(currentGroupId,i).pY=pY(currentGroupId,i,:);
                pBests(currentGroupId,i).NET=sum_speed;
                pBests(currentGroupId,i).Fitness=fitness;
                if fitness>max_fitness
                    max_fitness=fitness;
                    pp=i;
                    best_Points(currentGroupId).points=test_Points;
                end

                gBests(currentGroupId).pX=pX(currentGroupId,pp,:);
                gBests(currentGroupId).pY=pY(currentGroupId,pp,:);
                gBests(currentGroupId).NET=pBests(currentGroupId,pp).NET;
                gBests(currentGroupId).Fitness=pBests(currentGroupId,pp).Fitness;
                gBests(currentGroupId).id=currentGroupId;
                test_Points(:,4)=0;
            end
        else
            isCrossed=rand(1,serverNum);
            selectBestRatio=0.2;%Select the surviving population from the top 20% as the hybridization target and obtain new radio station positions from them
            crossGroupNum=ceil(surviveGroupChangeNum*selectBestRatio);
            crossFailTimes=0;%A maximum of 5 hybridization failures are allowed
            while crossFailTimes<5
                isCrossed=rand(1,serverNum);%Hybridization occurs with values less than 0.5, while others remain unchanged
                %Select hybridized objects
                randGroupList=index(1,1:crossGroupNum);
                randnum=randperm(length(randGroupList)); %Randomly generate matrix positions
                randGroup=randGroupList(randnum(1)); %Randomly extract one group from the data, which is the hybridized object
                if groupDieList(randGroup)==0
                    continue;
                end
                for serverIndex=1:serverNum
                    if isCrossed(serverIndex)<0.5
                        tmp_List=randBuildingList(randGroup).selectedBuilding;
                        selectedCrossedBuilding=tmp_List(serverIndex); %Extract a building from the list of buildings for hybridization
                        if ismember(selectedCrossedBuilding,cur_LowFitnessGroup)
                            continue;
                        end
                        cur_LowFitnessGroup(serverIndex)=selectedCrossedBuilding;
                    end
                end
                if ismember(cur_LowFitnessGroup,lowFitnessGroup,'rows')
                    crossFailTimes=crossFailTimes+1;
                    continue;
                end
                randBuildingList(currentGroupId).selectedBuilding=cur_LowFitnessGroup;
                break
            end

            %Reinitialize the population
            for i=1:pNum
                for j=1:serverNum
                    building_Index=cur_LowFitnessGroup(j);
                    tmp_Building_Pos=selected_Building_List(building_Index).pos;
                    [building_Area,~]=size(tmp_Building_Pos);
                    rand_Pos=randperm(building_Area,1);
                    tmp_P_Pos=tmp_Building_Pos(rand_Pos,:);

                    [lat, lon] = Ras2Point(tmp_P_Pos(1), tmp_P_Pos(2), img_X0, img_Y0, cell_width, cell_height);
                    pX(currentGroupId,i,j)=lat;
                    pY(currentGroupId,i,j)=lon;
                end
            end
            pVx(currentGroupId,:,:)=0;
            pVy(currentGroupId,:,:)=0;
            max_fitness=0;
            pp=1;
            for i=1:pNum
                for j =1:serverNum
                    [tmp_row,tmp_col,pZ]=Point2Ras(pX(currentGroupId,i,j), pY(currentGroupId,i,j), ...
                        img_X0, img_Y0, cell_width, cell_height, dsm_array);
                    if building_array(tmp_row,tmp_col) ~= not_Building
                        pZ=pZ+radio_Toground_Dis;
                    else
                        pZ=test_Points_Z+radio_Toground_Dis;
                    end

                    serverpoint=[pX(currentGroupId,i,j),pY(currentGroupId,i,j),pZ];

                    for k =1:test_Points_Num
                        tespoint=test_Points(k,:);

                        [building_dis,tree_dis]=GetLineNetSpeed(serverpoint(1),serverpoint(2),serverpoint(3), ...
                            tespoint(1),tespoint(2),tespoint(3),img_X0,img_Y0,cell_width,cell_height,dsm_array,building_array);
                        net_speed=floor(b0+b1*tree_dis+b2*building_dis+b3*tree_dis*tree_dis);
                        isLastPointPredicted=0;
                        last_building_dis=building_dis;
                        if net_speed>min_net_speed && radioUserRatio*net_speed>test_Points(k,4)
                            test_Points(k,4)=radioUserRatio*net_speed;
                        end

                    end
                end

                sum_speed=sum(test_Points(:,4))/test_Points_Num;
                fitness=sum_speed/(1+std(test_Points(:,4)));

                pBests(currentGroupId,i).pX=pX(currentGroupId,i,:);
                pBests(currentGroupId,i).pY=pY(currentGroupId,i,:);
                pBests(currentGroupId,i).NET=sum_speed;
                pBests(currentGroupId,i).Fitness=fitness;
                if fitness>max_fitness
                    max_fitness=fitness;
                    pp=i;
                    best_Points(currentGroupId).points=test_Points;
                end

                gBests(currentGroupId).pX=pX(currentGroupId,pp,:);
                gBests(currentGroupId).pY=pY(currentGroupId,pp,:);
                gBests(currentGroupId).NET=pBests(currentGroupId,pp).NET;
                gBests(currentGroupId).Fitness=pBests(currentGroupId,pp).Fitness;
                gBests(currentGroupId).id=currentGroupId;
                test_Points(:,4)=0;
            end
        end
        evolutionTimes=evolutionTimes+1;
    end
end
time=toc;
disp(['Fitness：' num2str(ggBest.Fitness)]);
disp(['Netspeed：' num2str(ggBest.NET)]);
disp(['Time：' num2str(time)]);
%% Output deployment location
col={'X' 'Y' 'Z'};
pZ=zeros(serverNum,1);
for i= 1:serverNum
    [tmp_row,tmp_col,h]=Point2Ras(ggBest.pX(1,1,i), ggBest.pY(1,1,i), img_X0, img_Y0, cell_width, cell_height,dsm_array);
    if building_array(tmp_row,tmp_col) ~= not_Building
        h=h+radio_Toground_Dis;
    else
        h=test_Points_Z+radio_Toground_Dis;
    end
    pZ(i)=h;
end

ANET_Pos_table=table((ggBest.pX(1,:))',(ggBest.pY(1,:))',pZ,'VariableNames',col);
writetable(ANET_Pos_table, 'SEPSO_ANET_Pos.csv');
col={'X' 'Y' 'Z' 'NET'};
result_table=table(best_Points(index(1) ...
    ).points(:,1),best_Points(index(1)).points(:,2),best_Points(index(1)).points(:,3),best_Points(index(1)).points(:,4),'VariableNames',col);

writetable(result_table, 'SEPSO_testPointsNet.csv');
