function [] = plotPopFracBars(y,contNo,hitBinNo,ax)

y = reshape(y,[hitBinNo,contNo+1])';

b = bar(ax,y,'stacked','LineWidth',1.5);
b(1).FaceColor = [255, 190, 11]/255;

for i = 2:hitBinNo
    colFac = i/(hitBinNo+0.5);
    b(i).FaceColor = [0.9,colFac*0.9,0.4+colFac*0.5];
end

ax.LineWidth = 1.5;
axis([0.4,contNo + 1.6,0,1])