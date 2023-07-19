function [GFPseg,RFPseg] = splitFluo(BFseg,GFPraw,RFPraw,GFPflat,RFPflat,fluoThresh)
%SPLITFLURO applies the specified fluorescence threshold to the given input
%data by taking the ratio between the fluorescence intensities and
%log-transforming the result.
%
%   INPUTS:
%       -BFseg: The segmentation of the brightfield image
%       -GFPraw, RFPraw: The raw fluorescence images for the two genotypes
%       -GFPflat, RFPflat: The background versions of the above two images,
%       used for flatfield correction.
%       -fluoThresh: The threshold to be applied to distinguish the two
%       populations of fluorescent cells.
%
%   OUTPUTS:
%       -GFPseg, RFPseg: The segmented versions of the two populations.
%
%   Author: Oliver J. Meacock, 2023

%Once you have the fluorescence threshold from the smoothed data,
%return to less smoothed data for higher resolution results
fluoRat = log(imgaussfilt(RFPraw./RFPflat,2)./imgaussfilt(GFPraw./GFPflat,2));

GFPsegTmp = and(fluoRat <= fluoThresh, BFseg);
RFPsegTmp = and(fluoRat > fluoThresh, BFseg);

se = strel('disk',2);

GFPseg = imopen(GFPsegTmp,se);
RFPseg = imopen(RFPsegTmp,se);