function velCourses = getExptVelocityCourses(timeList)
%GETEXPTVELOCITYCOURSES extracts the velocity timecourses from the PIV
%analyses of the experimental data.
%   
%   INPUTS:
%       -timeList: Times at which you want the velocities to be sampled.
%       Velocities will be linearly interpolated for times that were not
%       timepoints in the original data.
%
%   OUTPUTS:
%       -velCourses: Interpolated median velocities, measured by PIV. A
%       txcxr array, where t is the number of timepoints in timeList, c is
%       the number of conditions (starting densities) and r is the number
%       of replicates.
%
%   Author: Oliver J. Meacock, 2023

noiseMedianV = 0.003; %The measured value of the median PIV signal in frames, used to correct the measured sigal

PVD = readtable('processedVelocityDataOld.csv');

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