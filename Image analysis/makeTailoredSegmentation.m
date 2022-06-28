clear all
close all

Root = 'C:\Users\olijm\Desktop\SeanAna\Sample20x\OD0_01';

BFchan = 'Channel_1';
GFPchan = 'Channel_2';
RFPchan = 'Channel_3';

frameName = 'Frame_%04d.tif';
outName = 'Frame_%04d_Segmentation.tif';

lowDensFrames = 1:3; %Frames that will be used to do the flatfield correction for the three channels
maxT = 31;
pxSize = 0.227;

%Values you might want to adjust
microcolonyFrame = 18; %Frame by which microcolonies have developed, used to discount flecks of rubbish in the brightfield
maskMicrocolonies = false; %Whether to actually apply the microcolony mask (can actually reduce accuracy of segmentation if no debris in image)
tgtFrame = 15; %Frame you want to segment

textThresh = 2.5;
neighSize = 3; %Approximate size of a single cell
fluoThresh = 0.03; %Increase to increase relative fraction of green pixels

visBounds = [500,1000,500,1000]; %Boundaries of the region that will be displayed

%% Construct flatfield correction images
for i = lowDensFrames
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

%% Generate BF segmentation of microcolony mask
BFmask = splitBF(imread(fullfile(Root,BFchan,sprintf(frameName,microcolonyFrame))),neighSize,textThresh);
BFmask = imopen(BFmask,strel('disk',5));

%% Generate BF segmentation of target frame
BFraw = imread(fullfile(Root,BFchan,sprintf(frameName,tgtFrame)));
BFseg = splitBF(BFraw,neighSize,textThresh);

if maskMicrocolonies
    BFseg = and(BFmask,BFseg);
end

tiledlayout(1,3,'TileSpacing','tight')
nexttile()
imshow(cat(3,double(BFseg),double(BFraw)/double(max(BFraw(:))),double(BFraw)/double(max(BFraw(:)))),[])
title('BF segmentation')
axis(visBounds)

%% Generate fluorescence segmentation
GFPraw = double(imread(fullfile(Root,GFPchan,sprintf(frameName,tgtFrame))));
RFPraw = double(imread(fullfile(Root,RFPchan,sprintf(frameName,tgtFrame))));

[GFPseg,RFPseg] = splitFluo(BFseg,GFPraw,RFPraw,GFPflat,RFPflat,fluoThresh);

%% Display population segmentation
nexttile()
imshow(double(cat(3,RFPseg,GFPseg,zeros(size(RFPseg)))));
title('Fluo segmentation')
axis(visBounds)

nexttile()
RFPraw = imgaussfilt(RFPraw,2);
GFPraw = imgaussfilt(GFPraw,2);
RFPresc = (RFPraw-min(RFPraw(:)))/(max(RFPraw(:))-min(RFPraw(:)));
GFPresc = (GFPraw-min(GFPraw(:)))/(max(GFPraw(:))-min(GFPraw(:)));
imshow(cat(3,RFPresc,GFPresc,zeros(size(RFPraw))))
title('Original fluo')
axis(visBounds)

%% Run this code to save the image when you're happy with the segmentation
%imwrite(cat(3,RFPseg,GFPseg,zeros(size(RFPseg))),fullfile(Root,sprintf(outName,tgtFrame)))