clear all
close all

Root = 'C:\Users\olijm\Desktop\SeanAna\Sample20x';

BFchan = 'Channel_1';
GFPchan = 'Channel_2';
RFPchan = 'Channel_3';

frameName = 'Frame_%04d.tif';

lowDensFrames = 1:10; %Frames that will be used to do the flatfield correction for the three channels
maxT = 31;
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

pKs = zeros(size(BFstore,3),1);

for i = 1:maxT
    BFseg(:,:,i) = splitBF(BFstore(:,:,i),neighSize,textThresh);
    [GFPseg(:,:,i),RFPseg(:,:,i)] = splitFluo(BFseg(:,:,i),GFPstore(:,:,i),RFPstore(:,:,i),GFPflat,RFPflat);

    pKs(i) = measureSensitiveKillerContactProb(GFPseg(:,:,i),RFPseg(:,:,i),pxSize); 

    fprintf('Image %i of %i done.\n',i,maxT)
end