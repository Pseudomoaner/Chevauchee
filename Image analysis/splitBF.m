function BFseg = splitBF(BF,neighSize,textThresh)

%Apply texture analysis
stdImg = stdfilt(BF,ones(neighSize));
kernel = ones(neighSize) / neighSize^2; % Mean kernel
meanImg = conv2(BF, kernel, 'same'); % Convolve keeping size of I
seImg = stdImg./(meanImg.^0.5); %Logic here is that dividing by meanImg gives COV. However, the COV itself scales as one over the square root of the image intensity (shot noise), so multiply by that number to get just the contribution from the cells.
seImg = imgaussfilt(seImg,2);
BFseg = seImg > textThresh;