function [] = plotPopFracBars(y,ax)

noBins = size(y,2)/2;

y = [y(1:noBins);y(noBins+1:end)];

b = bar(ax,y,'stacked','LineWidth',1.5);
b(1).FaceColor = [0,0.4,1];

for i = 2:noBins
    b(i).FaceColor = [1,1-(i-2)/(noBins-1),0];
end

ax.LineWidth = 1.5;
axis([0.4,2.6,0,100])

xticklabels({'Protected','Attacked'})