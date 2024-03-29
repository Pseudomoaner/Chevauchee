function [] = plotPopFracBars(y,contNo,hitBinNo,ax)
%PLOTPOPFRACBARS generates a stacked bar chart representation of the
%populations present within a given lattice site.
%
%   INPUTS:
%       -y: List of populations associated with a given lattice site.
%       -contNo: Number of contact bins in simulation
%       -hitBinNo: Number of hit bins in simulation
%       -ax: Active axes for visualisation to be inserted into.
%
%   Author: Oliver J. Meacock, 2023

y = reshape(y,[hitBinNo,contNo+1])';

b = bar(ax,y,'stacked','LineWidth',1.5);
b(1).FaceColor = [255, 120, 10]/255; %Original yellow = [255, 190, 10];

for i = 2:hitBinNo
    colFac = i/(hitBinNo+1.5);
    b(i).FaceColor = [1,0.5 + colFac*0.5,0.05+colFac*0.95]; %Original magenta = [0.9,colFac*0.9,0.4+colFac*0.5]
end

ax.LineWidth = 1.5;
axis(ax,[0.4,contNo + 1.6,0,1])

tickLabs = {};
for i = 1:contNo + 1
    tickLabs = [tickLabs;num2str(i-1)];
end

ax.XTickLabel = tickLabs;

xlabel(ax,'Number of sensitive contacts to attackers')