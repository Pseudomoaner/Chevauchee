function [] = plotPopFracBars(y,ax)

noBins = size(y,2)/5;

y = reshape(y,[noBins,5])';

b = bar(ax,y,'stacked','LineWidth',1.5);
b(1).FaceColor = [255, 190, 11]/255;

for i = 2:noBins
    colFac = i/(noBins+0.5);
    b(i).FaceColor = [0.9,colFac*0.9,0.4+colFac*0.5];
end

ax.LineWidth = 1.5;
axis([0.4,5.6,0,100])

xticklabels({'oooo','ooo.','oo..','o...','....'})