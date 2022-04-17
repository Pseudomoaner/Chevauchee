% batch processor for analyzing contacts and coverage from fluorescent
% timelapse images

clear all
close all

%set actual root, pull root folders from inside
rootFold = 'D:\Sean\SurfaceColonyPIV\Fluorescence_Blocks\';

d = dir(rootFold);
folds = d([d(:).isdir]) ;
folds = folds(~ismember({folds(:).name},{'.','..'}));
folds = {folds.name};

for loop = 1:size(folds,2);
    disp(folds(loop));
    Root = strcat(rootFold,char(folds(loop)));
%from 3 channel time stacks, use imageJ macro to save in separate folders
%'Channel_X\Frame_0000.tif'
BFchan = 'Channel_1';
GFPchan = 'Channel_2';
RFPchan = 'Channel_3';

frameName = 'Frame%04d.tif'; %imageJ output was done without _
outName = 'Analysis_V2.mat';

lowDensFrames = 1:3; %Frames that will be used to do the flatfield correction for the three channels
microcolonyFrame = 15; %Frame when microcolonies have developed, used to discount flecks of rubbish in the brightfield
maxT = length(dir(fullfile(Root,BFchan))) -2; %count number of files in folder since this number varies between experiments. Minus 2 because of matlab stupidity
pxSize = 0.227;
neighSize = 5; %Approximate size of a single cell
textThresh = 2.5;

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

BFseg = zeros(size(BFstore));
GFPseg = zeros(size(BFstore));
RFPseg = zeros(size(BFstore));

packFracs = zeros(size(BFstore,3),1);
pKs = zeros(size(BFstore,3),1);

BFmask = splitBF(BFstore(:,:,microcolonyFrame),neighSize,textThresh);
BFmask = imopen(BFmask,strel('disk',10));

for i = 1:maxT
    if i < microcolonyFrame
        BFseg(:,:,i) = and(splitBF(BFstore(:,:,i),neighSize,textThresh),BFmask);
    else
        BFseg(:,:,i) = imopen(splitBF(BFstore(:,:,i),neighSize,textThresh),strel('disk',10));
    end
    [GFPseg(:,:,i),RFPseg(:,:,i)] = splitFluo(BFseg(:,:,i),GFPstore(:,:,i),RFPstore(:,:,i),GFPflat,RFPflat);

    pKs(i) = measureSensitiveKillerContactProb(GFPseg(:,:,i),RFPseg(:,:,i),pxSize); 
    packFracs(i) = sum(sum(BFseg(:,:,i)))/(size(BFseg,1)*size(BFseg,2));

    fprintf('Image %i of %i done.\n',i,maxT)
end

save(fullfile(Root,outName),'BFseg','GFPseg','RFPseg','pKs','packFracs')
clearvars -except rootFold folds %otherwise errors happen
end