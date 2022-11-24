function velCourses = getExptVelocityCourses(timeList)

PVD = readtable('C:\Users\omeacock\Desktop\SeanData\processedVelocityData.csv');

ExptNames = unique(PVD{:,4});

velCourses = zeros(size(timeList,2),4,3);
repCnt = ones(4,1);

for i = 1:size(ExptNames,1)
    currDat = PVD{strcmp(PVD{:,4},ExptNames{i}),[2:3,5]};
    densInd = -log10(currDat(1,1))+1;
    velCourses(:,densInd,repCnt(densInd)) = interp1(currDat(:,3),currDat(:,2),timeList);
    repCnt(densInd) = repCnt(densInd) + 1;
end