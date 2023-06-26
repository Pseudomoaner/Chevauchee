function velCourses = getExptVelocityCourses(timeList)

noiseMedianV = 0.003; %The measured value of the median PIV signal in frames, used to correct the measured sigal

PVD = readtable('C:\Users\Olivier\OneDrive - Université de Lausanne\Simulations\Driveby\processedVelocityDataOld.csv');

ExptNames = unique(PVD{:,4});

velCourses = zeros(size(timeList,2),4,3);
repCnt = ones(4,1);

for i = 1:size(ExptNames,1)
    currDat = PVD{strcmp(PVD{:,4},ExptNames{i}),[2:3,5]};
    densInd = -log10(currDat(1,1))+1;
    velCourses(:,densInd,repCnt(densInd)) = interp1(currDat(:,3),currDat(:,2),timeList);
    repCnt(densInd) = repCnt(densInd) + 1;
end

velCourses = mean(velCourses,3);

%Apply noise correction to the measured median velocities, under the
%assumption that the velocity distribution is a Reyleigh and the noise
%is Gaussian
velCourses = sqrt(velCourses.^2 - noiseMedianV^2);