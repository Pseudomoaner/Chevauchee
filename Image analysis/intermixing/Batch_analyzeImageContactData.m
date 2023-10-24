%Script that calculates the packing fractions and intermixing indices for 
%all datasets within a specified directory.

clear all
close all

%set actual root, pull root folders from inside
rootFold = 'D:\Sean\SurfaceColonyPIV\Fluorescence_Blocks\';

d = dir(rootFold);
folds = d([d(:).isdir]) ;
folds = folds(~ismember({folds(:).name},{'.','..'}));
folds = {folds.name};

for loop = 1:size(folds,2)
    disp(folds(loop));
    Root = strcat(rootFold,char(folds(loop)));
    %from 3 channel time stacks, use imageJ macro to save in separate folder

    %'Channel_X\Frame_0000.tif'
    BFchan = 'Channel_1';
    GFPchan = 'Channel_2';
    RFPchan = 'Channel_3';

    frameName = 'Frame%04d.tif'; %imageJ output was done without _
    outName = 'Analysis.mat';

    lowDensFrames = 1:3; %Frames that will be used to do the flatfield correction for the three channels
    maxT = length(dir(fullfile(Root,BFchan))) -2; %count number of files in folder since this number varies between experiments. Minus 2 because of matlab stupidity
    pxSize = 0.227;
    
    %% Construct flatfield correction images
    for i = 1:maxT
        BF = imread(fullfile(Root,BFchan,sprintf(frameName,i-1)));
        GFP = imread(fullfile(Root,GFPchan,sprintf(frameName,i-1)));
        RFP = imread(fullfile(Root,RFPchan,sprintf(frameName,i-1)));

        BFstore(:,:,i) = double(BF);
        GFPstore(:,:,i) = double(GFP);
        RFPstore(:,:,i) = double(RFP);
    end

    BFflat = imgaussfilt(mean(BFstore(:,:,lowDensFrames),3),15);
    GFPflat = imgaussfilt(mean(GFPstore(:,:,lowDensFrames),3),50);
    RFPflat = imgaussfilt(mean(RFPstore(:,:,lowDensFrames),3),50);

    %% Perfrom brightfield and fluorescence segmentation
    BFseg = splitBFseries(BFstore);

    GFPseg = zeros(size(BFstore));
    RFPseg = zeros(size(BFstore));
    fluoThresh = zeros(maxT,1);

    packFracs = zeros(size(BFstore,3),1);
    varMix = zeros(size(BFstore,3),1);

    BFmask = BFseg(:,:,microcolonyFrame);
    BFmask = imopen(BFmask,strel('disk',10));

    %Get initial fluorescence threshold estimates, prior to smoothing
    for i = 1:maxT
        if i < microcolonyFrame
            BFseg(:,:,i) = and(BFseg(:,:,i),BFmask);
        else
            BFseg(:,:,i) = imopen(BFseg(:,:,i),strel('disk',10));
        end
        fluoThresh(i) = findFluoThresh(BFseg(:,:,i),GFPstore(:,:,i),RFPstore(:,:,i),GFPflat,RFPflat);
    end

    %Smooth fluorescence threshold temporally and apply to images to segment
    %populations
    fluoThresh = smooth(fluoThresh,5);
    for i = 1:maxT
        [GFPseg(:,:,i),RFPseg(:,:,i)] = splitFluo(BFseg(:,:,i),GFPstore(:,:,i),RFPstore(:,:,i),GFPflat,RFPflat,fluoThresh(i));

        varMix(i) = measureSensitiveKillerIntermixing(GFPseg(:,:,i),RFPseg(:,:,i),pxSize);
        packFracs(i) = sum(sum(BFseg(:,:,i)))/(size(BFseg,1)*size(BFseg,2));
    end

    save(fullfile(Root,outName),'BFseg','GFPseg','RFPseg','varMix','packFracs')
    clearvars -except rootFold folds %otherwise errors happen
end