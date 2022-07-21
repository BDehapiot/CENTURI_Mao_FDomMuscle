clearvars
%%% FDomMuscle_1D
%%% Benoit Dehapiot, PhD
%%% benoit.dehapiot@univ-amu.fr
%%% CENTURI (Turing Center for Living Systems)
%%% Multi-Engineering Platform
%%% Aix Marseille Université

%% Options

    %%% Inputs %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
    
    RootPath = pwd;
    DataPath = strcat(RootPath, '\data');
    OutputsPath = strcat(RootPath, '\outputs');
    pixSize = 0.0830266; % pixel size in µm
    preload = 0; % Directly open a previously saved FiberList
    
    %%% Parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    Channel = 'C1';
    minSD = 0; % Min. S.D. value between low & upBound to be considered as valid pattern
    minProm = 0.6; % % Min. peak promincence to be considered as valid pattern
    fillGap = 10; % Fill gap if un-patterned region size is < to fillGap (µm)
    minSize = 6; % Discard patterned region smaller than minSize (µm)
    
    BoundLow = 1.0; BoundUp = 3.0; % Lower and upper bound for peak analysis (µm)
    
    nBin = 3; % Number of vertical bins to analyse local auto-correlation
    wAvgSize1 = 20; % Walking average size for local auto-correlation analysis (µm)
    wAvgSize2 = 2; % Post-processing walking average size for local auto-correlation display (µm)
    winSize = 20; % Window size for displaying local auto-correlation analysis (µm)

    % winSize can't be > wAvgSize!!!
    
%% Initialize

    % DirList & Foldlist
    DirList = dir(DataPath); 
    DirList(1:4,:) = []; % Remove ghost files
    FiberList = cell(0,0);
    for i=1:size(DirList,1)
        if DirList(i).isdir == 1
            FiberList{end+1,1} = DirList(i).name;
        end
    end
    
    % Get Parameters in pixels 
    wAvgSize1Pix = floor(wAvgSize1/pixSize);
    wAvgSize2Pix = floor(wAvgSize2/pixSize);
    winSizePix = floor(winSize/pixSize);
    BoundLowPix = floor(BoundLow/pixSize);
    BoundUpPix = ceil(BoundUp/pixSize);
    fillGapPix = fillGap/pixSize;
    minSizePix = minSize/pixSize;
        
%% Open data

if preload == 1 
    % Directly open a previously saved FiberList 
    FiberList = load(strcat(DataPath,'\FiberList.mat'));
    FiberList = FiberList.FiberList;
elseif preload == 0
    % Create a new FiberList 
    for i=1:size(FiberList,1)
        tempList = dir(strcat(DataPath,'\',FiberList{i,1}));
        tempCell = cell(0,0);
        for j=1:size(tempList,1)
            if tempList(j).bytes > 10000 % Remove ghost files
                tempCell{end+1,1} = tempList(j).name;
                tempPath = strcat(DataPath,'\',FiberList{i,1},'\',tempList(j).name);
                tempCell{end,2} = uint16(loadtiff(tempPath,1));
            end
        end
        FiberList{i,2} = tempCell; clear tempCell
    end
end

%% Process data

if preload == 2 
    % Directly open a previously saved FiberList 
    FiberList = load(strcat(DataPath,'\FiberList.mat'));
    FiberList = FiberList.FiberList;
end
if preload < 2
    parfor i=1:size(FiberList,1)
        for j=1:size(FiberList{i,2},1)
            if contains(FiberList{i,2}{j,1},Channel) == 1 % Check for channel

                % Get tempFiber
                tempFiber = FiberList{i,2}{j,2};
                FiberSize = size(tempFiber,2);

                % Substract Background
                tempBlur = imgaussfilt(tempFiber,20);
                tempProcess = tempFiber-tempBlur;
                tempProcess = imgaussfilt(tempProcess,1);

                % Make 1D average and HighPass IIR
                tempAvg = mean(tempProcess,1)';
                tempIIR = customIIR(tempAvg,300,0,1,0);
             
                % Measure AutoCorr
                tempAutoCorr = customAutoCorr(tempIIR,tempIIR,'msubmean');
                [~,tempMid] = max(tempAutoCorr);
                temp1 = flipud(tempAutoCorr(1:tempMid));
                temp2 = tempAutoCorr(tempMid:end);
                tempAutoCorr = horzcat(temp1,temp2);
                tempAutoCorr = mean(tempAutoCorr,2);
                tempxAutoCorr = (0:pixSize:(size(tempAutoCorr,1)*pixSize)-pixSize)';

                % Extract vertical bins
                BinSize = size(tempProcess,1)/nBin;
                tempAvgBin = NaN(FiberSize,nBin);
                tempIIRBin = NaN(FiberSize,nBin);
                for k=1:nBin
                    tempAvgBin(:,k) = mean(tempProcess((k*BinSize)-(BinSize-1):k*BinSize,:))';
                    tempIIRBin(:,k) = customIIR(tempAvgBin(:,k),300,0,1,0);
                end
                
                % Get local AutoCorr per bin
                MergedAutoCorrLocal = NaN(winSizePix,FiberSize,nBin);
                MergedAutoCorrLocalSD = NaN(FiberSize,nBin);
                MaxSDAutoCorrLocal = NaN(winSizePix,FiberSize);
                for k=1:nBin
                    for m=1:FiberSize
                        if m > wAvgSize1Pix/2 && m < FiberSize-((wAvgSize1Pix/2)-1)
                            tempIIRLocal = tempIIRBin(m-wAvgSize1Pix/2:m+((wAvgSize1Pix/2)-1),k);
                            tempACCLocal = customAutoCorr(tempIIRLocal,tempIIRLocal,'msubmean');
                            [~,tempMid] = max(tempACCLocal);
                            temp1 = flipud(tempACCLocal(tempMid-winSizePix+1:tempMid,1));
                            temp2 = tempACCLocal(tempMid:tempMid+winSizePix-1,1);
                            temp3 = mean(horzcat(temp1,temp2),2);
                            MergedAutoCorrLocal(:,m,k) = temp3;
                            MergedAutoCorrLocalSD(m,k) = std(temp3(BoundLowPix:BoundUpPix,1)); % Measure SD between low & upBound
                        end
                    end
                end
                
                % Get local AutoCorr of MaxSD
                [MaxSD,MaxSDx] = max(MergedAutoCorrLocalSD,[],2);
                for m=1:FiberSize
                    MaxSDAutoCorrLocal(:,m) = MergedAutoCorrLocal(:,m,MaxSDx(m,1));
                end
                MaxSDAutoCorrLocal = movmean(MaxSDAutoCorrLocal,wAvgSize2Pix,2);
                
%                 % Get local AutoCorr "parameters" (pks, locs, width, prom)
%                 MergedPeaksinfo = NaN(FiberSize,4);
%                 for m=1:FiberSize
%                     [pks,locs,width,prom] = findpeaks(MaxSDAutoCorrLocal(BoundLowPix:BoundUpPix,m),'MinPeakDistance',6);
%                     if ~isempty(pks)
%                         [maxprom,idxprom] = max(prom);
%                         MergedPeaksinfo(m,1) = pks(idxprom,1);
%                         MergedPeaksinfo(m,2) = (locs(idxprom,1)*pixSize)+BoundLow; % Add back the BoundLow
%                         MergedPeaksinfo(m,3) = width(idxprom,1)*pixSize;
%                         MergedPeaksinfo(m,4) = prom(idxprom,1);
%                     else
%                         MergedPeaksinfo(m,1:4) = NaN;
%                     end
%                 end  
                
                % Get local AutoCorr "parameters" (pks, locs, width, prom)
                MergedPeaksinfo = NaN(FiberSize,4);
                for m=1:FiberSize
                    [pks,locs,width,prom] = findpeaks(MaxSDAutoCorrLocal(BoundLowPix:BoundUpPix,m),'MinPeakDistance',6);
                    if ~isempty(pks)
                        temp_locs = (locs*pixSize)+BoundLow;
                        temp_locs = abs(temp_locs-2);
                        [minlocs,idxlocs] = min(temp_locs);
                        MergedPeaksinfo(m,1) = pks(idxlocs,1);
                        MergedPeaksinfo(m,2) = (locs(idxlocs,1)*pixSize)+BoundLow; % Add back the BoundLow
                        MergedPeaksinfo(m,3) = width(idxlocs,1)*pixSize;
                        MergedPeaksinfo(m,4) = prom(idxlocs,1);
                    else
                        MergedPeaksinfo(m,1:4) = NaN;
                    end
                end  
                
                % Append Fiberlist
                FiberList{i,2}{j,3} = tempProcess;
                FiberList{i,2}{j,4}(:,1) = tempAvg;
                FiberList{i,2}{j,4}(:,2) = tempIIR;
                FiberList{i,2}{j,4}(:,3) = tempAvg-tempIIR;
                FiberList{i,2}{j,4}(:,4) = MaxSD;
                FiberList{i,2}{j,5}(:,1) = tempxAutoCorr;
                FiberList{i,2}{j,5}(:,2) = tempAutoCorr;
                FiberList{i,2}{j,6} = MaxSDAutoCorrLocal;
                FiberList{i,2}{j,7} = MergedPeaksinfo;
                
            end
        end
    end
    save(strcat(DataPath,'\FiberList.mat'),'FiberList','-v7.3');
end

%% Get Results

% Read ClassTable
ClassTable = table2cell(readtable(strcat(DataPath,'\ClassTable.csv')));

FiberResults = cell(0,0);
% Get results
for i=1:size(FiberList,1)
    for j=1:size(FiberList{i,2},1)
        if contains(FiberList{i,2}{j,1},Channel) == 1 % Check for channel

                    % Extract data from FiberList
                    FiberResults{end+1,1} = str2num(FiberList{i,1}(:,5:6)); % Day
                    FiberResults{end,2} = FiberList{i,1}; % Image
                    FiberResults{end,3} = FiberList{i,2}{j,1}; % ROI
                    FiberResults{end,4} = FiberList{i,2}{j,2}; % Raw
                    FiberResults{end,5} = FiberList{i,2}{j,3}; % Processed
                    FiberResults{end,6} = FiberList{i,2}{j,5}(1:120,:); % AutoCorr
                    FiberResults{end,7} = FiberList{i,2}{j,4}; % Data1D (part#1)
                    FiberResults{end,7}(:,5:8) = FiberList{i,2}{j,7}; % Data1D (part#2)

                    % Extract valid pattern ( >= minSD & >= minProm )
                    BorderSize = (wAvgSize1Pix+wAvgSize2Pix)/2;
                    Data1D = FiberResults{end,7}(BorderSize+1:end-BorderSize-1,:);
                    FiberSize = size(Data1D,1);            
                    for m=1:FiberSize
                        if Data1D(m,4) >= minSD && Data1D(m,8) >= minProm
                            Data1D_Valid_Idx(m,1) = 1;
                        else
                            Data1D_Valid_Idx(m,1) = 0;
                        end
                    end

                    % Fill gap if un-patterned region size is < to fillGap (µm)
                    CC = bwconncomp(imcomplement(Data1D_Valid_Idx),4);
                    S1 = regionprops(CC,'Area');              
                    for m=1:size(S1,1)
                        if S1(m).Area <= fillGapPix
                            Data1D_Valid_Idx(CC.PixelIdxList{1,m}) = 1;
                        end
                    end 
                    
                    % Discard patterned region smaller than minSize (µm)
                    CC = bwconncomp(Data1D_Valid_Idx,4);
                    S1 = regionprops(CC,'Area');              
                    for m=1:size(S1,1)
                        if S1(m).Area <= minSizePix
                            Data1D_Valid_Idx(CC.PixelIdxList{1,m}) = 0;
                        end
                    end                     

                    % Get Data1D_Valid
                    for m=1:FiberSize
                        if Data1D_Valid_Idx(m,1) == 1
                            Data1D_Valid(m,1:8) = Data1D(m,1:8);
                        else
                            Data1D_Valid(m,1:8) = NaN;
                        end
                    end
                    Valid_Count = sum(~isnan(Data1D_Valid(:,1)),1);
                    CC = bwconncomp(Data1D_Valid_Idx,4);
                    S1 = regionprops(CC,'Area');
                    tempNaN1 = NaN(BorderSize,8); tempNaN2 = NaN(1,8);
                    Data1D_Valid = vertcat(tempNaN1,Data1D_Valid,tempNaN1,tempNaN2);
                    clear tempNaN1 tempNaN2

                    % Append FiberResults
                    FiberResults{end,8} = Data1D_Valid; % Valid Data1D ( >= minSD & >= minProm )
                    FiberResults{end,9} = FiberList{i,2}{j,6}; % Local auto-correlation
                    FiberResults{end,10} = Valid_Count/FiberSize; % Valid percentage
                    FiberResults{end,11} = (CC.NumObjects/(FiberSize*pixSize))*1000; % Avg. number of valid regions (per mm)
                    FiberResults{end,12} = mean(cat(1,S1.Area))*pixSize; % Avg. length of valid regions
                    FiberResults{end,13} = nanmean(Data1D_Valid(:,6)); % Avg. locs of valid pattern
                    FiberResults{end,14} = nanstd(Data1D_Valid(:,6)); % S.D. locs of valid pattern
                    FiberResults{end,15} = nanmean(Data1D_Valid(:,7)); % Avg. width of valid pattern 
                    FiberResults{end,16} = nanmean(Data1D_Valid(:,8)); % Avg. prom of valid pattern

                    clear BorderSize FiberSize Valid_Count CC S1 Data1D Data1D_Valid Data1D_Valid_Idx  

        end
    end
end

% Add class to FiberResults 
ClassTable = table2cell(readtable(strcat(DataPath,'\ClassTable.csv')));
for i=1:size(FiberResults,1)
    tempName = FiberResults{i,3}(:,1:end-7);
    for j=1:size(ClassTable,1)
        if contains(ClassTable{j,5},tempName) == 1 
            if ClassTable{j,6} == 'N'
                FiberResults{i,17} = 1;
            elseif ClassTable{j,6} == 'I'
                FiberResults{i,17} = 2;
            elseif ClassTable{j,6} == 'T'
                FiberResults{i,17} = 3;
            elseif ClassTable{j,6} == 'M'
                FiberResults{i,17} = 4;
            end
%             FiberResults{i,18} = ClassTable{j,8};
        end
    end
end
clear tempName

%% Data Preparation (All data)

MergedValidPercentage = NaN(size(FiberResults,1),15);
MergedValidNumber = NaN(size(FiberResults,1),15);
MergedValidLength = NaN(size(FiberResults,1),15);
MergedValidLocs = NaN(size(FiberResults,1),15);
MergedValidLocsSD = NaN(size(FiberResults,1),15);
MergedValidProm = NaN(size(FiberResults,1),15);
MergedValidWidth = NaN(size(FiberResults,1),15);
MergedClass = NaN(size(FiberResults,1),15);
for i=1:size(FiberResults,1)
    tempDay = FiberResults{i,1};
    MergedValidPercentage(i,tempDay) = FiberResults{i,10};
    MergedValidNumber(i,tempDay) = FiberResults{i,11};
    MergedValidLength(i,tempDay) = FiberResults{i,12};
    MergedValidLocs(i,tempDay) = FiberResults{i,13};
    MergedValidLocsSD(i,tempDay) = FiberResults{i,14};
    MergedValidWidth(i,tempDay) = FiberResults{i,15};
    MergedValidProm(i,tempDay) = FiberResults{i,16};
    MergedClass(i,tempDay) = FiberResults{i,17};
end
clear tempDay

for i=1:size(MergedValidPercentage,2)
    MergedValidCount(i,1) = nanmean(MergedValidPercentage(:,i));
    MergedValidCount(i,2) = 1-nanmean(MergedValidPercentage(:,i));
end

tempCount = MergedClass;
tempCount(~isnan(tempCount)) = 1;
tempCount(isnan(tempCount)) = 0;
for i=1:size(MergedClass,2)
    MergedClassCount(i,1) = sum(MergedClass(:,i)==1);
    MergedClassCount(i,2) = sum(MergedClass(:,i)==2);
    MergedClassCount(i,3) = sum(MergedClass(:,i)==3);
    MergedClassCount(i,4) = sum(MergedClass(:,i)==4);
end
for i=1:size(MergedClass,2)
    MergedClassCountNorm(i,:) = MergedClassCount(i,:)/nansum(tempCount(:,i));
end
clear tempCount

%% Plot sections : "uncomment" and run ------------------------------------

%% *** Plot *** "FiberAutoCorr" per Day and Class

% % Per day
% 
%     Days = [2,3,4,5,7,9,12,15]';
%     FiberAutoCorr_Day = cell(0,0);
%     for j=1:size(Days,1)
%         tempAutoCorr = NaN(120,0);
%         for i=1:size(FiberResults,1)         
%             if FiberResults{i,1} == Days(j,1)
%                 tempAutoCorr(:,end+1) = FiberResults{i,6}(:,2);
%             end
%         end
%         FiberAutoCorr_Day{j,1} = Days(j,1);
%         FiberAutoCorr_Day{j,2} = tempAutoCorr;
%         FiberAutoCorr_Day{j,3}(:,1) = nanmean(tempAutoCorr,2);
%         FiberAutoCorr_Day{j,3}(:,2) = nanstd(tempAutoCorr,2);
%         FiberAutoCorr_Day{1,4}(:,1) = FiberResults{1,6}(:,1); % x vector
%         FiberAutoCorr_Day{1,4}(:,j+1) = nanmean(tempAutoCorr,2);
%         FiberAutoCorr_Day{1,4}(:,j+1+size(Days,1)) = nanstd(tempAutoCorr,2);
%     end
% 
%     % display
%     FiberAutoCorr_Day_Display = NaN(120,0);
%     for j=1:size(Days,1)
%         tempn = size(FiberAutoCorr_Day{j,2},2);
%         Idx = (1:tempn)';
%         Idx = Idx(randperm(size(Idx,1)));
%         Idx = Idx(1:32,1);
%         for i=1:32
%             FiberAutoCorr_Day_Display(:,end+1) = FiberAutoCorr_Day{j,2}(:,Idx(i,1));
%         end
%     end    
%     figure(1);
%     h = pcolor(FiberAutoCorr_Day_Display(1:75,:)); % plot
%     set(h,'EdgeColor','none');
%     colorbar
%     caxis([-0.6 0.6])
%     clear tempAutoCorr Idx tempn h
%     
%     % Legend 
%     FiberAutoCorr_Day = array2table(FiberAutoCorr_Day,'VariableNames',...
%     {'Day','AutoCorr','AutoCorr_Avg','AutoCorr_Avg_Merged'});
% 
% % Per category
% 
%     Cat = [1,2,3,4]';
%     CatNames = ["Non-sarcomeric"; "Immature"; "Transitional"; "Mature"];
%     FiberAutoCorr_Cat = cell(0,0);
%     for j=1:size(Cat,1)
%         tempAutoCorr = NaN(120,0);
%         for i=1:size(FiberResults,1)         
%             if FiberResults{i,17} == Cat(j,1)
%                 tempAutoCorr(:,end+1) = FiberResults{i,6}(:,2);
%             end
%         end
%         FiberAutoCorr_Cat{j,1} = CatNames(j,1);
%         FiberAutoCorr_Cat{j,2} = tempAutoCorr;
%         FiberAutoCorr_Cat{j,3}(:,1) = nanmean(tempAutoCorr,2);
%         FiberAutoCorr_Cat{j,3}(:,2) = nanstd(tempAutoCorr,2);
%         FiberAutoCorr_Cat{1,4}(:,1) = FiberResults{1,6}(:,1); % x vector
%         FiberAutoCorr_Cat{1,4}(:,j+1) = nanmean(tempAutoCorr,2);
%         FiberAutoCorr_Cat{1,4}(:,j+1+size(Cat,1)) = nanstd(tempAutoCorr,2);
%     end
%     
%     % display
%     FiberAutoCorr_Cat_Display = NaN(120,0);
%     for j=1:size(Cat,1)
%         tempn = size(FiberAutoCorr_Cat{j,2},2);
%         Idx = (1:tempn)';
%         Idx = Idx(randperm(size(Idx,1)));
%         Idx = Idx(1:36,1);
%         for i=1:36
%             FiberAutoCorr_Cat_Display(:,end+1) = FiberAutoCorr_Cat{j,2}(:,Idx(i,1));
%         end
%     end
%     figure(2);
%     h = pcolor(FiberAutoCorr_Cat_Display(1:75,:)); % plot
%     set(h,'EdgeColor','none');
%     colorbar
%     caxis([-0.6 0.6])
%     clear tempAutoCorr Idx tempn h
%     
%     % Legend 
%     FiberAutoCorr_Cat = array2table(FiberAutoCorr_Cat,'VariableNames',...
%     {'Cat','AutoCorr','AutoCorr_Avg','AutoCorr_Avg_Merged'});
% 
%     clear i j Days Cat CatNames
   
%% *** Plot *** "MergedResults" per Class

% %General variables
% nDays = 15;
% nClass = 4;
% nFibers = size(FiberResults,1);
% 
% %Set Y axis
% minY1 = 0; maxY1 = 1; nTick1 = 6;  % MergedValidPercentage
%     YLim1 = [minY1-(maxY1*0.05) maxY1+maxY1*0.05]; 
%     YTick1 = (minY1:(maxY1-minY1)/(nTick1-1):maxY1);
%     
% minY2 = 0; maxY2 = 50; nTick2 = 6;  % MergedValidNumber
%     YLim2 = [minY2-(maxY2*0.05) maxY2+maxY2*0.05]; 
%     YTick2 = (minY2:(maxY2-minY2)/(nTick2-1):maxY2);
%     
% minY3 = 0; maxY3 = 200; nTick3 = 6;  % MergedValidLength
%     YLim3 = [minY3-(maxY3*0.05) maxY3+maxY3*0.05]; 
%     YTick3 = (minY3:(maxY3-minY3)/(nTick3-1):maxY3);
%     
% minY4 = 1; maxY4 = 3; nTick4 = 5;  % MergedValidLocs
%     YLim4 = [minY4-(maxY4*0.05) maxY4+maxY4*0.05]; 
%     YTick4 = (minY4:(maxY4-minY4)/(nTick4-1):maxY4);
%     
% minY5 = 0; maxY5 = 0.6; nTick5 = 7;  % MergedValidLocsSD
%     YLim5 = [minY5-(maxY5*0.05) maxY5+maxY5*0.05]; 
%     YTick5 = (minY5:(maxY5-minY5)/(nTick5-1):maxY5);
%     
% minY6 = 0.0; maxY6 = 1.0; nTick6 = 5; % MergedValidProm
%     YLim6 = [minY6-(maxY6*0.05) maxY6+maxY6*0.05]; 
%     YTick6 = (minY6:(maxY6-minY6)/(nTick6-1):maxY6);
%     
%     
% %Set box plots appearance
% DotSize = 20;
% SquareSize = 100;
% ClassLabels = {'Myotube','Immature','Transitional','Mature'};
% DotClassColor(1,:) = [102/255 255/255 102/255];
% DotClassColor(2,:) = [51/255 204/255 51/255];
% DotClassColor(3,:) = [0/255 153/255 51/255];
% DotClassColor(4,:) = [0/255 51/255 0/255];
% DotClassColor(5,:) = [0.80 0.80 0.80];
% 
% %Spread dots on box plots
% DataX = ones(numel(MergedValidPercentage),nClass); 
% r = -0.25+(0.25+0.25).*rand(numel(MergedValidPercentage),1);
% for i=1:nClass
%     DataX(:,i) = (DataX(:,i)*i) + r;
% end
% clear r
% 
% figure('Renderer', 'painters', 'Position', [50 100 1800 1200])
% 
% % MergedValidPercentage
%     subplot(3,2,1)  
%     for i=1:nClass
%         tempX = DataX(:,i);
%         tempData = MergedValidPercentage;
%         tempData(MergedClass~=i) = NaN;
%         tempData = reshape(tempData,[numel(tempData),1]);
%         tempDataMerged(:,i) = tempData;
%         scatter(tempX,tempData,DotSize,'filled','MarkerFaceColor',DotClassColor(5,:));
%         hold on
%         scatter(i,nanmean(tempData(:)),SquareSize,'filled','s','black');
%     end
%     boxplot(tempDataMerged,'symbol','','Labels',ClassLabels); MergedValidPercentage_class = tempDataMerged;
%     set(findobj(gca,'type','line'),'Color','k','LineWidth',1,'LineStyle','-')
%     set(findobj(gca,'type','line','Tag','Median'),'Color','k','LineWidth',1.5,'LineStyle','-')
%     ax = gca; set(gca,'box','off'); ax.TickDir = 'out'; ax.FontSize = 12;
%     ax.XLim = [0.5 4.5]; ax.YLim = YLim1; ax.YTick = YTick1;
%     title('Patterned proportion (per fiber)')
%     
% % MergedValidNumber
%     subplot(3,2,2) 
%     for i=1:nClass
%         tempX = DataX(:,i);
%         tempData = MergedValidNumber;
%         tempData(MergedClass~=i) = NaN;
%         tempData = reshape(tempData,[numel(tempData),1]);
%         tempDataMerged(:,i) = tempData;
%         scatter(tempX,tempData,DotSize,'filled','MarkerFaceColor',DotClassColor(5,:));
%         hold on
%         scatter(i,nanmean(tempData(:)),SquareSize,'filled','s','black');
%     end
%     boxplot(tempDataMerged,'symbol','','Labels',ClassLabels); MergedValidNumber_class = tempDataMerged;
%     set(findobj(gca,'type','line'),'Color','k','LineWidth',1,'LineStyle','-')
%     set(findobj(gca,'type','line','Tag','Median'),'Color','k','LineWidth',1.5,'LineStyle','-')
%     ax = gca; set(gca,'box','off'); ax.TickDir = 'out'; ax.FontSize = 12;
%     ax.XLim = [0.5 4.5]; ax.YLim = YLim2; ax.YTick = YTick2;
%     title('Patterned regions number (per mm of fiber)')
%         
% % MergedValidLength
%     subplot(3,2,3) 
%     for i=1:nClass
%         tempX = DataX(:,i);
%         tempData = MergedValidLength;
%         tempData(MergedClass~=i) = NaN;
%         tempData = reshape(tempData,[numel(tempData),1]);
%         tempDataMerged(:,i) = tempData;
%         scatter(tempX,tempData,DotSize,'filled','MarkerFaceColor',DotClassColor(5,:));
%         hold on
%         scatter(i,nanmean(tempData(:)),SquareSize,'filled','s','black');
%     end
%     boxplot(tempDataMerged,'symbol','','Labels',ClassLabels); MergedValidLength_class = tempDataMerged;
%     set(findobj(gca,'type','line'),'Color','k','LineWidth',1,'LineStyle','-')
%     set(findobj(gca,'type','line','Tag','Median'),'Color','k','LineWidth',1.5,'LineStyle','-')
%     ax = gca; set(gca,'box','off'); ax.TickDir = 'out'; ax.FontSize = 12;
%     ax.XLim = [0.5 4.5]; ax.YLim = YLim3; ax.YTick = YTick3;
%     title('Length of patterned regions (µm)')
%     
% % MergedValidLocs
%     subplot(3,2,4) 
%     for i=1:nClass
%         tempX = DataX(:,i);
%         tempData = MergedValidLocs;
%         tempData(MergedClass~=i) = NaN;
%         tempData = reshape(tempData,[numel(tempData),1]);
%         tempDataMerged(:,i) = tempData;
%         scatter(tempX,tempData,DotSize,'filled','MarkerFaceColor',DotClassColor(5,:));
%         hold on
%         scatter(i,nanmean(tempData(:)),SquareSize,'filled','s','black');
%     end
%     boxplot(tempDataMerged,'symbol','','Labels',ClassLabels); MergedValidLocs_class = tempDataMerged;
%     set(findobj(gca,'type','line'),'Color','k','LineWidth',1,'LineStyle','-')
%     set(findobj(gca,'type','line','Tag','Median'),'Color','k','LineWidth',1.5,'LineStyle','-')
%     ax = gca; set(gca,'box','off'); ax.TickDir = 'out'; ax.FontSize = 12;
%     ax.XLim = [0.5 4.5]; ax.YLim = YLim4; ax.YTick = YTick4;
%     title('Main corr. peak localization (µm)')
%     
% % MergedValidLocsSD
%     subplot(3,2,5) 
%     for i=1:nClass
%         tempX = DataX(:,i);
%         tempData = MergedValidLocsSD;
%         tempData(MergedClass~=i) = NaN;
%         tempData = reshape(tempData,[numel(tempData),1]);
%         tempDataMerged(:,i) = tempData;
%         scatter(tempX,tempData,DotSize,'filled','MarkerFaceColor',DotClassColor(5,:));
%         hold on
%         scatter(i,nanmean(tempData(:)),SquareSize,'filled','s','black');
%     end
%     boxplot(tempDataMerged,'symbol','','Labels',ClassLabels); MergedValidLocsSD_class = tempDataMerged; 
%     set(findobj(gca,'type','line'),'Color','k','LineWidth',1,'LineStyle','-')
%     set(findobj(gca,'type','line','Tag','Median'),'Color','k','LineWidth',1.5,'LineStyle','-')
%     ax = gca; set(gca,'box','off'); ax.TickDir = 'out'; ax.FontSize = 12;
%     ax.XLim = [0.5 4.5]; ax.YLim = YLim5; ax.YTick = YTick5;
%     title('Main corr. peak localization S.D.')
%     
% % MergedValidProm
%     subplot(3,2,6) 
%     for i=1:nClass
%         tempX = DataX(:,i);
%         tempData = MergedValidProm;
%         tempData(MergedClass~=i) = NaN;
%         tempData = reshape(tempData,[numel(tempData),1]);
%         tempDataMerged(:,i) = tempData;
%         scatter(tempX,tempData,DotSize,'filled','MarkerFaceColor',DotClassColor(5,:));
%         hold on
%         scatter(i,nanmean(tempData(:)),SquareSize,'filled','s','black');
%     end
%     boxplot(tempDataMerged,'symbol','','Labels',ClassLabels); MergedValidProm_class = tempDataMerged; 
%     set(findobj(gca,'type','line'),'Color','k','LineWidth',1,'LineStyle','-')
%     set(findobj(gca,'type','line','Tag','Median'),'Color','k','LineWidth',1.5,'LineStyle','-')
%     ax = gca; set(gca,'box','off'); ax.TickDir = 'out'; ax.FontSize = 12;
%     ax.XLim = [0.5 4.5]; ax.YLim = YLim6; ax.YTick = YTick6;
%     title('Main corr. peak prominence')
%     
%     clear i j ax nDays nClass nFibers
%     clear -regexp temp minY maxY nTick YLim YTick
%     clear SquareSize DotSize DotClassColor DataX

%% *** Plot *** "MergedResults" per Day
  
% % General variables
% nDays = 15;
% nClass = 4;
% nFibers = size(FiberResults,1);
% 
% % Set Y axis
% minY1 = 0; maxY1 = 1; nTick1 = 6;  % MergedValidPercentage
%     YLim1 = [minY1-(maxY1*0.05) maxY1+maxY1*0.05]; 
%     YTick1 = (minY1:(maxY1-minY1)/(nTick1-1):maxY1);
%     
% minY2 = 0; maxY2 = 50; nTick2 = 6;  % MergedValidNumber
%     YLim2 = [minY2-(maxY2*0.05) maxY2+maxY2*0.05]; 
%     YTick2 = (minY2:(maxY2-minY2)/(nTick2-1):maxY2);
%     
% minY3 = 0; maxY3 = 200; nTick3 = 6;  % MergedValidLength
%     YLim3 = [minY3-(maxY3*0.05) maxY3+maxY3*0.05]; 
%     YTick3 = (minY3:(maxY3-minY3)/(nTick3-1):maxY3);
%     
% minY4 = 1; maxY4 = 3; nTick4 = 5;  % MergedValidLocs
%     YLim4 = [minY4-(maxY4*0.05) maxY4+maxY4*0.05]; 
%     YTick4 = (minY4:(maxY4-minY4)/(nTick4-1):maxY4);
%     
% minY5 = 0; maxY5 = 0.6; nTick5 = 7;  % MergedValidLocsSD
%     YLim5 = [minY5-(maxY5*0.05) maxY5+maxY5*0.05]; 
%     YTick5 = (minY5:(maxY5-minY5)/(nTick5-1):maxY5);
%     
% minY6 = 0.0; maxY6 = 1.0; nTick6 = 5; % MergedValidProm
%     YLim6 = [minY6-(maxY6*0.05) maxY6+maxY6*0.05]; 
%     YTick6 = (minY6:(maxY6-minY6)/(nTick6-1):maxY6);
% 
% % Set box plots appearance
% DotSize = 20;
% SquareSize = 100;
% DayLabels = {'d01','d02','d03','d04','d05','d06','d07','d08','d09','d10','d11','d12','d13','d14','d15'};
% DotClassColor(1,:) = [102/255 255/255 102/255];
% DotClassColor(2,:) = [51/255 204/255 51/255];
% DotClassColor(3,:) = [0/255 153/255 51/255];
% DotClassColor(4,:) = [0/255 51/255 0/255];
% DotClassColor(5,:) = [0.80 0.80 0.80];
%     
% % Spread dots on box plots
% DataX = ones(nFibers,nDays); 
% r = -0.33+(0.33+0.33).*rand(nFibers,1);
% for i=1:nDays
%     DataX(:,i) = (DataX(:,i)*i) + r;
% end
% clear r
% 
% figure('Renderer', 'painters', 'Position', [50 100 1800 1200])
% 
%     % MergedValidPercentage
%     subplot(3,2,1)  
%     tempDataMerged = MergedValidPercentage;
%     tempClassMean = NaN(4,1);
%     for i=1:nDays
%         tempX = DataX(:,i);
%         for j=1:nClass
%             tempData = tempDataMerged(:,i);
%             tempClass = MergedClass(:,i);
%             tempData(tempClass~=j) = NaN;
%             tempClassMean(j,1) = nanmean(tempData);
%             scatter(tempX,tempData,DotSize,'filled','MarkerFaceColor',DotClassColor(5,:));
%             hold on            
%         end
%         for j=1:nClass
%             scatter(i+0.5,tempClassMean(j,1),SquareSize,'*','MarkerEdgeColor',DotClassColor(j,:),'LineWidth',1.5);
%         end
%         scatter(i,nanmean(tempDataMerged(:,i)),SquareSize,'filled','s','black');
%     end
%     boxplot(tempDataMerged,'symbol','','Labels',DayLabels);
%     set(findobj(gca,'type','line'),'Color','k','LineWidth',1,'LineStyle','-')
%     set(findobj(gca,'type','line','Tag','Median'),'Color','k','LineWidth',1.5,'LineStyle','-')
%     ax = gca; set(gca,'box','off'); ax.TickDir = 'out'; ax.FontSize = 12;
%     ax.XLim = [0.5 15.5]; ax.YLim = YLim1; ax.YTick = YTick1;
%     title('Patterned proportion (per fiber)')
%     
%     % MergedValidNumber
%     subplot(3,2,2)  
%     tempDataMerged = MergedValidNumber;
%     tempClassMean = NaN(4,1);
%     for i=1:nDays
%         tempX = DataX(:,i);
%         for j=1:nClass
%             tempData = tempDataMerged(:,i);
%             tempClass = MergedClass(:,i);
%             tempData(tempClass~=j) = NaN;
%             tempClassMean(j,1) = nanmean(tempData);
%             scatter(tempX,tempData,DotSize,'filled','MarkerFaceColor',DotClassColor(5,:));
%             hold on
%         end
%         for j=1:nClass
%             scatter(i+0.5,tempClassMean(j,1),SquareSize,'*','MarkerEdgeColor',DotClassColor(j,:),'LineWidth',1.5);
%         end
%         scatter(i,nanmean(tempDataMerged(:,i)),SquareSize,'filled','s','black');
%     end
%     boxplot(tempDataMerged,'symbol','','Labels',DayLabels);
%     set(findobj(gca,'type','line'),'Color','k','LineWidth',1,'LineStyle','-')
%     set(findobj(gca,'type','line','Tag','Median'),'Color','k','LineWidth',1.5,'LineStyle','-')
%     ax = gca; set(gca,'box','off'); ax.TickDir = 'out'; ax.FontSize = 12;
%     ax.XLim = [0.5 15.5]; ax.YLim = YLim2; ax.YTick = YTick2;
%     title('Patterned regions number (per mm of fiber)')
%     
%     % MergedValidLength
%     subplot(3,2,3)  
%     tempDataMerged = MergedValidLength;
%     tempClassMean = NaN(4,1);
%     for i=1:nDays
%         tempX = DataX(:,i);
%         for j=1:nClass
%             tempData = tempDataMerged(:,i);
%             tempClass = MergedClass(:,i);
%             tempData(tempClass~=j) = NaN;
%             tempClassMean(j,1) = nanmean(tempData);
%             scatter(tempX,tempData,DotSize,'filled','MarkerFaceColor',DotClassColor(5,:));
%             hold on
%         end
%         for j=1:nClass
%             scatter(i+0.5,tempClassMean(j,1),SquareSize,'*','MarkerEdgeColor',DotClassColor(j,:),'LineWidth',1.5);
%         end
%         scatter(i,nanmean(tempDataMerged(:,i)),SquareSize,'filled','s','black');
%     end
%     boxplot(tempDataMerged,'symbol','','Labels',DayLabels);
%     set(findobj(gca,'type','line'),'Color','k','LineWidth',1,'LineStyle','-')
%     set(findobj(gca,'type','line','Tag','Median'),'Color','k','LineWidth',1.5,'LineStyle','-')
%     ax = gca; set(gca,'box','off'); ax.TickDir = 'out'; ax.FontSize = 12;
%     ax.XLim = [0.5 15.5]; ax.YLim = YLim3; ax.YTick = YTick3;
%     title('Length of patterned regions (µm)')
%     
%         
%     % MergedValidLocs
%     subplot(3,2,4)  
%     tempDataMerged = MergedValidLocs;
%     tempClassMean = NaN(4,1);
%     for i=1:nDays
%         tempX = DataX(:,i);
%         for j=1:nClass
%             tempData = tempDataMerged(:,i);
%             tempClass = MergedClass(:,i);
%             tempData(tempClass~=j) = NaN;
%             tempClassMean(j,1) = nanmean(tempData);
%             scatter(tempX,tempData,DotSize,'filled','MarkerFaceColor',DotClassColor(5,:));
%             hold on
%         end
%         for j=1:nClass
%             scatter(i+0.5,tempClassMean(j,1),SquareSize,'*','MarkerEdgeColor',DotClassColor(j,:),'LineWidth',1.5);
%         end
%         scatter(i,nanmean(tempDataMerged(:,i)),SquareSize,'filled','s','black');
%     end
%     boxplot(tempDataMerged,'symbol','','Labels',DayLabels);
%     set(findobj(gca,'type','line'),'Color','k','LineWidth',1,'LineStyle','-')
%     set(findobj(gca,'type','line','Tag','Median'),'Color','k','LineWidth',1.5,'LineStyle','-')
%     ax = gca; set(gca,'box','off'); ax.TickDir = 'out'; ax.FontSize = 12;
%     ax.XLim = [0.5 15.5]; ax.YLim = YLim4; ax.YTick = YTick4;
%     title('Main corr. peak localization (µm)')
%     
%     % MergedValidLocsSD
%     subplot(3,2,5)  
%     tempDataMerged = MergedValidLocsSD;
%     tempClassMean = NaN(4,1);
%     for i=1:nDays
%         tempX = DataX(:,i);
%         for j=1:nClass
%             tempData = tempDataMerged(:,i);
%             tempClass = MergedClass(:,i);
%             tempData(tempClass~=j) = NaN;
%             tempClassMean(j,1) = nanmean(tempData);
%             scatter(tempX,tempData,DotSize,'filled','MarkerFaceColor',DotClassColor(5,:));
%             hold on
%         end
%         for j=1:nClass
%             scatter(i+0.5,tempClassMean(j,1),SquareSize,'*','MarkerEdgeColor',DotClassColor(j,:),'LineWidth',1.5);
%         end
%         scatter(i,nanmean(tempDataMerged(:,i)),SquareSize,'filled','s','black');
%     end
%     boxplot(tempDataMerged,'symbol','','Labels',DayLabels);
%     set(findobj(gca,'type','line'),'Color','k','LineWidth',1,'LineStyle','-')
%     set(findobj(gca,'type','line','Tag','Median'),'Color','k','LineWidth',1.5,'LineStyle','-')
%     ax = gca; set(gca,'box','off'); ax.TickDir = 'out'; ax.FontSize = 12;
%     ax.XLim = [0.5 15.5]; ax.YLim = YLim5; ax.YTick = YTick5;
%     title('Main corr. peak localization S.D.')
%     
%     % MergedValidProm
%     subplot(3,2,6)  
%     tempDataMerged = MergedValidProm;
%     tempClassMean = NaN(4,1);
%     for i=1:nDays
%         tempX = DataX(:,i);
%         for j=1:nClass
%             tempData = tempDataMerged(:,i);
%             tempClass = MergedClass(:,i);
%             tempData(tempClass~=j) = NaN;
%             tempClassMean(j,1) = nanmean(tempData);
%             scatter(tempX,tempData,DotSize,'filled','MarkerFaceColor',DotClassColor(5,:));
%             hold on
%         end
%         for j=1:nClass
%             scatter(i+0.5,tempClassMean(j,1),SquareSize,'*','MarkerEdgeColor',DotClassColor(j,:),'LineWidth',1.5);
%         end
%         scatter(i,nanmean(tempDataMerged(:,i)),SquareSize,'filled','s','black');
%     end
%     boxplot(tempDataMerged,'symbol','','Labels',DayLabels);
%     set(findobj(gca,'type','line'),'Color','k','LineWidth',1,'LineStyle','-')
%     set(findobj(gca,'type','line','Tag','Median'),'Color','k','LineWidth',1.5,'LineStyle','-')
%     ax = gca; set(gca,'box','off'); ax.TickDir = 'out'; ax.FontSize = 12;
%     ax.XLim = [0.5 15.5]; ax.YLim = YLim6; ax.YTick = YTick6;
%     title('Main corr. peak prominence')
%     
%     clear i j ax nDays nClass nFibers
%     clear -regexp temp minY maxY nTick YLim YTick
%     clear SquareSize DotSize DotClassColor DataX

%% *** Plot *** "LocalAutoCorr" Select #fiber 
    
%     SelectFiber = 650;
%     SelectDay = FiberResults{SelectFiber,1};
%     SelectCat = FiberResults{SelectFiber,17};
%     
%     %Get variables
%     FiberImage = FiberResults{SelectFiber,4};
%     LocalAutoCorr = FiberResults{SelectFiber,9};
%     Data1D = FiberResults{SelectFiber,7};
%     Data1DValid = FiberResults{SelectFiber,8};
%     FiberSize = size(Data1D,1);
% 
%     %Set pcolor axis
%     X = repmat((0:FiberSize-1)*pixSize,120,1);
%     Y = repmat(((0:120-1)')*pixSize,1,FiberSize);
% 
%     %Plot
%     subplot(4,1,1)
%         FiberImage = imresize(FiberImage,[90 FiberSize]);    
%         imshow(imadjust(FiberImage))
%         set(gcf,'Position',[0 0 FiberSize/2 900])
%         title(strcat('Day#',num2str(SelectDay),' Cat#',num2str(SelectCat),' Fiber#',num2str(SelectFiber)))
% 
%     subplot(4,1,2) 
%         h = pcolor(X,Y,LocalAutoCorr);    
%         caxis([-1.0 1.0])
%         axis([0 FiberSize*pixSize 0 10])
%         set(h,'EdgeColor','none');
%         set(gca,'Layer','top');
%         ax = gca; ax.TickDir = 'out';
%     hold on
%         plot(X(1,:),Data1DValid(:,6),'-k')
% 
%     subplot(4,1,3)  
%         plot(X(1,:),Data1D(:,4),'-r')
%         yline(minSD,':r','minSD','LabelHorizontalAlignment','left');
%         axis([0 FiberSize*pixSize 0 1])
%         ax = gca; ax.TickDir = 'out';
% 
%     subplot(4,1,4)  
%         plot(X(1,:),Data1D(:,8),'-b')
%         yline(minProm,':b','minProm','LabelHorizontalAlignment','left');
%         axis([0 FiberSize*pixSize 0 2])
%         ax = gca; ax.TickDir = 'out';
%          
% %     % Save image
% %     tempName = FiberResults{SelectFiber,3};
% % 	FileName = strcat(tempName(:,1:end-4),'_LocalAutoCorr.eps');
% %     set(gcf,'renderer','Painters')
% %     saveas(h,strcat(OutputsPath,'\',FileName),'epsc');
%     
%     clear i j h ans ax X Y
%     clear SelectFiber SelectDay SelectCat
%     clear Data1D Data1DValid FiberImage FiberSize
%     clear -regexp temp 

%% *** Plot *** "LocalAutoCorr" Select #Day & #Cat

% SelectDay = 2;
% SelectCat = 2;
% nFibers = size(FiberResults,1);
% 
% randFiber = (1:nFibers)';
% randFiber = randFiber(randperm(nFibers));
% 
% for i=1:nFibers
% 
%     % Get testFiber
%     testFiber = randFiber(i,1); 
%     nDay = FiberResults{testFiber,1};
%     nClass = FiberResults{testFiber,17};
%     
%     if nDay == SelectDay && nClass == SelectCat
%     
%     % Get variables
%     FiberImage = FiberResults{testFiber,4};
%     LocalAutoCorr = FiberResults{testFiber,9};
%     Data1D = FiberResults{testFiber,7};
%     Data1DValid = FiberResults{testFiber,8};
%     FiberSize = size(Data1D,1);
% 
%     % Set pcolor axis
%     X = repmat((0:FiberSize-1)*pixSize,120,1);
%     Y = repmat(((0:120-1)')*pixSize,1,FiberSize);
% 
%     % Plot
%     subplot(4,1,1)
%         FiberImage = imresize(FiberImage,[90 FiberSize]);    
%         imshow(imadjust(FiberImage))
%         set(gcf,'Position',[0 0 FiberSize/2 900])
%         title(strcat('Day#',num2str(SelectDay),' Cat#',num2str(SelectCat),' Fiber#',num2str(randFiber(i,1))))
% 
%     subplot(4,1,2) 
%         h = pcolor(X,Y,LocalAutoCorr);    
%         caxis([-1.0 1.0])
%         axis([0 FiberSize*pixSize 0 10])
%         set(h,'EdgeColor','none');
%         set(gca,'Layer','top');
%         ax = gca; ax.TickDir = 'out';
%     hold on
%         plot(X(1,:),Data1DValid(:,6),'-k')
% 
%     subplot(4,1,3)  
%         plot(X(1,:),Data1D(:,4),'-r')
%         yline(minSD,':r','minSD','LabelHorizontalAlignment','left');
%         axis([0 FiberSize*pixSize 0 1])
%         ax = gca; ax.TickDir = 'out';
% 
%     subplot(4,1,4)  
%         plot(X(1,:),Data1D(:,8),'-b')
%         yline(minProm,':b','minProm','LabelHorizontalAlignment','left');
%         axis([0 FiberSize*pixSize 0 2])
%         ax = gca; ax.TickDir = 'out';
%         
%         drawnow
%         waitforbuttonpress;
%         close all
%         
%     end 
% end
% 
% clear i j h ans ax X Y
% clear SelectDay SelectCat randFiber
% clear nClass nDay nFibers testFiber
% clear Data1D Data1DValid FiberImage FiberSize

%% *** Plot *** "LocalAutoCorr" Save all

% nFibers = size(FiberResults,1);
% 
% for i=1:nFibers
% 
%     % Get nDay and nClass
%     nDay = FiberResults{i,1};
%     nClass = FiberResults{i,17};
%     
%     % Get variables
%     FiberImage = FiberResults{i,4};
%     LocalAutoCorr = FiberResults{i,9};
%     Data1D = FiberResults{i,7};
%     Data1DValid = FiberResults{i,8};
%     FiberSize = size(Data1D,1);
% 
%     % Set pcolor axis
%     X = repmat((0:FiberSize-1)*pixSize,120,1);
%     Y = repmat(((0:120-1)')*pixSize,1,FiberSize);
% 
%     % Plot
%     h = figure;
%     h = figure('visible','off');
%     subplot(4,1,1)
%         FiberImage = imresize(FiberImage,[90 FiberSize]);    
%         imshow(imadjust(FiberImage))
%         set(gcf,'Position',[0 0 FiberSize/2 900])
%         title(FiberResults{i,3})
% 
%     subplot(4,1,2) 
%         h = pcolor(X,Y,LocalAutoCorr);    
%         caxis([-1.0 1.0])
%         axis([0 FiberSize*pixSize 0 10])
%         set(h,'EdgeColor','none');
%         set(gca,'Layer','top');
%         ax = gca; ax.TickDir = 'out';
%     hold on
%         plot(X(1,:),Data1DValid(:,6),'-k')
% 
%     subplot(4,1,3)  
%         plot(X(1,:),Data1D(:,4),'-r')
%         yline(minSD,':r','minSD','LabelHorizontalAlignment','left');
%         axis([0 FiberSize*pixSize 0 1])
%         ax = gca; ax.TickDir = 'out';
% 
%     subplot(4,1,4)  
%         plot(X(1,:),Data1D(:,8),'-b')
%         yline(minProm,':b','minProm','LabelHorizontalAlignment','left');
%         axis([0 FiberSize*pixSize 0 2])
%         ax = gca; ax.TickDir = 'out';
%         
%     tempName = FiberResults{i,3};
% 	FileName = strcat(tempName(:,1:end-4),'_LocalAutoCorr.jpg');
%     saveas(h,strcat(OutputsPath,'\',FileName),'jpeg');
%     
%     close all
% 
% end
% 
% clear i j h ans ax X Y nDay nClass nFibers
% clear Data1D Data1DValid FiberImage FiberSize
% clear -regexp temp 

%% Extract 1D valid patterns

% nFibers = size(FiberResults,1);
% 
% % ValidPattern1D 
% 
%     % All at ones
%     ValidPattern1D = NaN(8000,nFibers);
%     for i=1:nFibers
%         temp = FiberResults{i,8}(:,5);
%         temp(isnan(temp)) = 0;
%         ValidPattern1D(1:size(temp,1),i) = logical(temp);
%     end
% 
%     % One by one
%     for i=1:nFibers
%         temp = FiberResults{i,8}(:,5); % What you gonna save 
%         temp(isnan(temp)) = 0; % 1D Valid
%         temp = logical(temp); % 1D Valid
%         tempName = FiberResults{i,3};
%         csvwrite(strcat(OutputsPath,'\', tempName(:,1:end-4),'_Valid1D.csv'),temp); % You can modify files name accordingly (blue writting)
%     end
%     
% % PromPattern1D
% 
%     % All at ones
%     PromPattern1D = NaN(8000,nFibers);
%     for i=1:nFibers
%         temp = FiberResults{i,8}(:,8);
%         PromPattern1D(1:size(temp,1),i) = temp;
%     end
% 
%     % One by one
%     for i=1:nFibers
%         temp = FiberResults{i,8}(:,8); % What you gonna save 
%         tempName = FiberResults{i,3};
%         csvwrite(strcat(OutputsPath,'\', tempName(:,1:end-4),'_Prom1D.csv'),temp); % You can modify files name accordingly (blue writting)
%     end
% 
% clear i j nFibers
% clear PromPattern1D ValidPattern1D
% clear temp tempName

%% Legend FiberResults

FiberResultsArray = FiberResults;

for i=1:size(FiberResultsArray,1)
    tempData1D = FiberResultsArray{i,7};
    tempData1D = array2table(tempData1D,'VariableNames',...
        {'Profile1D','HighPassIIR','Subtracted','SD','pks','locs','width','prom'});
    FiberResultsArray{end,7} = tempData1D;
    clear tempData1D
end
FiberResultsArray = array2table(FiberResultsArray,'VariableNames',...
    {'Day','Image','ROI','Raw','Processed','AutoCorr','Data1D','Data1DValid','LocalAutoCorr','ValidPercent','ValidNumber_mm','ValidLength','AvgLocs','AvgLocsSD','AvgWidth','AvgProm','Class'});

%% Clear
clear i j k m preload
clear nBin wAvgSize1 wAvgSize1Pix wAvgSize2 wAvgSize2Pix winSize winSizePix BoundLow BoundLowPix BoundUp BoundUpPix fillGap fillGapPix minSizePix