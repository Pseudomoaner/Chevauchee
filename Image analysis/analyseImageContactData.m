clear all
close all

Root = 'C:\Users\olijm\Desktop\SeanAna\Sample20x\OD0_01';

BFchan = 'Channel_1';
GFPchan = 'Channel_2';
RFPchan = 'Channel_3';

frameName = 'Frame_%04d.tif';
outName = 'Analysis.mat';

lowDensFrames = 1:3; %Frames that will be used to do the flatfield correction for the three channels
microcolonyFrame = 15; %Frame by which microcolonies have developed, used to discount flecks of rubbish in the brightfield
maxT = 31;
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

packFracs = zeros(size(BFstore,3),1);
pKs = zeros(size(BFstore,3),1);

BFmask = BFseg(:,:,microcolonyFrame);
BFmask = imopen(BFmask,strel('disk',10));

for i = 1:maxT
    if i < microcolonyFrame
        BFseg(:,:,i) = and(BFseg(:,:,i),BFmask);
    else
        BFseg(:,:,i) = imopen(BFseg(:,:,i),strel('disk',10));
    end
    [GFPseg(:,:,i),RFPseg(:,:,i)] = splitFluo(BFseg(:,:,i),GFPstore(:,:,i),RFPstore(:,:,i),GFPflat,RFPflat);

    pKs(i) = measureSensitiveKillerContactProb(GFPseg(:,:,i),RFPseg(:,:,i),pxSize); 
    packFracs(i) = sum(sum(BFseg(:,:,i)))/(size(BFseg,1)*size(BFseg,2));

    fprintf('Image %i of %i done.\n',i,maxT)
end

save(fullfile(Root,outName),'BFseg','GFPseg','RFPseg','pKs','packFracs')