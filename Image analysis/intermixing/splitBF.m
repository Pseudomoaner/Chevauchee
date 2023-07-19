function BFseg = splitBF(BF,neighSize,textThresh)
%SPLITBF splits a given brightfield image into cell-containing and cell
%non-containing regions using a texture-based approach.
%
%   INPUTS:
%       -BF: The original brightfield image. A single frame
%       -neighSize: The size of the neighbourhood within which the texture
%       metric should be evaluated. Generally of the same size (in pixels)
%       as a single cell.
%       -textThresh: The threshold used to split the resulting texture
%       metric into foreground and background.
%
%   OUTPUTS:
%       -BFseg: The segmented brightfield frame.
%
%   Author: Oliver J. Meacock, 2023

%Apply texture analysis
stdImg = stdfilt(BF,ones(neighSize));
kernel = ones(neighSize) / neighSize^2; % Mean kernel
meanImg = conv2(BF, kernel, 'same'); % Convolve keeping size of I
seImg = stdImg./(meanImg.^0.5); %Logic here is that dividing by meanImg gives COV. However, the COV itself scales as one over the square root of the image intensity (shot noise), so multiply by that number to get just the contribution from the cells.
seImg = imgaussfilt(seImg,2);
BFseg = seImg > textThresh;