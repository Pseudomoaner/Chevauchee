function BFseg = splitBFseries(BFstore)

textThresh = 2.0;
neighSize = 5; %Approximate size of a single cell

BFseg = zeros(size(BFstore));

%Sweep 1: A rough sweep, using the default parameter values
for i = 1:size(BFstore,3)
    BFseg(:,:,i) = splitBF(BFstore(:,:,i),neighSize,textThresh);
end

BFfrac = squeeze(sum(sum(BFseg,1),2))/(size(BFseg,1)*size(BFseg,2));

%Sweep 2: A more refined sweep, varying the texture threshold on outlying
%timepoints so the packing fraction trace is smoothed out
margin = 0.05;

interMean = (BFfrac(1:end-2) + BFfrac(3:end))/2; % Prediction of packing fraction from neighbouring timepoints

highInds = find(BFfrac(2:end-1) > interMean * (1+margin)) + 1; %Plus one to counteract unit offset in calculation of interMean
lowInds = find(BFfrac(2:end-1) < interMean * (1-margin)) + 1;

for i = highInds'
    currFrame = BFseg(:,:,i);
    currTextThresh = textThresh;
    while sum(currFrame(:))/(size(BFseg,1)*size(BFseg,2)) > interMean(i-1) * (1+margin)
        currTextThresh = currTextThresh*1.05;
        currFrame = splitBF(BFstore(:,:,i),neighSize,currTextThresh);
    end
    BFseg(:,:,i) = currFrame;
end

for i = lowInds'
    currFrame = BFseg(:,:,i);
    currTextThresh = textThresh;
    while sum(currFrame(:))/(size(BFseg,1)*size(BFseg,2)) < interMean(i-1) * (1-margin)
        currTextThresh = currTextThresh*0.95;
        currFrame = splitBF(BFstore(:,:,i),neighSize,currTextThresh);
    end
    BFseg(:,:,i) = currFrame;
end