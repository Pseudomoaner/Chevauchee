clear all
close all

Root = 'C:\Users\olijm\Desktop\SeanAna\Sample20x\OD0_01';

CMchan = 'Channel_0';
GFPchan = 'Channel_2';
RFPchan = 'Channel_3';

frameName = 'Frame%04d.tif';
outName = 'Analysis.mat';

lowDensFrames = 1:3; %Frames that will be used to do the flatfield correction for the three channels
microcolonyFrame = 15; %Frame by which microcolonies have developed, used to discount flecks of rubbish in the brightfield
maxT = length(dir(fullfile(Root,CMchan))) -2; 
pxSize = 0.227;

%% Construct flatfield correction images
for i = 1:maxT
    BF = imread(fullfile(Root,CMchan,sprintf(frameName,i-1)));
    GFP = imread(fullfile(Root,GFPchan,sprintf(frameName,i-1)));
    RFP = imread(fullfile(Root,RFPchan,sprintf(frameName,i-1)));

    BFstore(:,:,i) = double(BF);
    GFPstore(:,:,i) = double(GFP);
    RFPstore(:,:,i) = double(RFP);
end

GFPflat = imgaussfilt(mean(GFPstore(:,:,lowDensFrames),3),50);
RFPflat = imgaussfilt(mean(RFPstore(:,:,lowDensFrames),3),50);

%% Perfrom brightfield and fluorescence segmentation
BFseg = BFstore;

GFPseg = zeros(size(BFstore));
RFPseg = zeros(size(BFstore));

packFracs = zeros(size(BFstore,3),1);
pKs = zeros(size(BFstore,3),1);

for i = 1:maxT
    [GFPseg(:,:,i),RFPseg(:,:,i)] = splitFluo_correlationMap(BFseg(:,:,i),GFPstore(:,:,i),RFPstore(:,:,i),GFPflat,RFPflat);

    pKs(i) = measureSensitiveKillerContactProb(GFPseg(:,:,i),RFPseg(:,:,i),pxSize); 
    packFracs(i) = sum(sum(BFseg(:,:,i)))/(size(BFseg,1)*size(BFseg,2));

    fprintf('Image %i of %i done.\n',i,maxT)
end

save(fullfile(Root,outName),'BFseg','GFPseg','RFPseg','pKs','packFracs')