function [GFPseg,RFPseg] = splitFluo(BFseg,GFPraw,RFPraw,GFPflat,RFPflat,fluoThresh)

%Once you have the fluorescence threshold from the smoothed data,
%return to less smoothed data for higher resolution results
fluoRat = log(imgaussfilt(RFPraw./RFPflat,2)./imgaussfilt(GFPraw./GFPflat,2));

GFPsegTmp = and(fluoRat <= fluoThresh, BFseg);
RFPsegTmp = and(fluoRat > fluoThresh, BFseg);

se = strel('disk',2);

GFPseg = imopen(GFPsegTmp,se);
RFPseg = imopen(RFPsegTmp,se);