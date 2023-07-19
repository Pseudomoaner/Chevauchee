% export contact and packFrac data for analysis in R

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
    load(fullfile(Root,"Analysis.mat"),'pKs','packFracs')
    
    output = table(pKs,packFracs);
    writetable(output, fullfile(Root,"output.csv"));
end