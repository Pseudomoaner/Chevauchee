function [] = plotStdAsArea(tVals,data,col,ax)

dataStd = nanstd(data,0,2);
dataStd(isnan(dataStd)) = [];
dataMean = nanmean(data,2);

areaY = [dataMean-dataStd,2*dataStd];

h = area(ax,tVals,areaY,'EdgeColor','none');

h(1).FaceColor = col;
h(2).FaceColor = col;
h(1).FaceAlpha = 0;
h(2).FaceAlpha = 0.25;