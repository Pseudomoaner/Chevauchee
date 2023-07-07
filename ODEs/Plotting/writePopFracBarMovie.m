function [] = writePopFracBarMovie(y,t,noConts,noHitBins)
%WRITEPOPFRACBARMOVIE creates a movie showing the changes in the population
%distribution of a given bin of a continuum simulation.
%
%   INPUTS:
%       - y: The population distribution in the specified bin over time
%       - t: The list of times at which the system was sampled.
%       - noConts: The number of contact bins used (typically 5)
%       - noHitBins: The number of hit accumulation bins used.
%
%   Author: Oliver J. Meacock, 2023

interpy = interp1(t,y(:,2:end),2:2:max(t),'pchip');

figure()
fig = gcf;
fig.Units = 'Normalized';
fig.Position=[0.2,0.2,0.3,0.3];

tiledlayout(1,5,'TileSpacing','compact')
axHist1A = nexttile;
axHist1 = nexttile([1,4]);
ax = gca;

writer = VideoWriter('C:\Users\olijm\Desktop\SeanAna\test.avi');
writer.FrameRate = 60;
open(writer)

for i = 1:size(interpy,1)
    bar(axHist1A,y(1,1),'FaceColor',[0.2,0.6,1],'LineWidth',1.5)
    axis(axHist1A,[0.5,1.5,0,1]);
    axHist1A.LineWidth = 1.5;
    axHist1A.Box = 'off';
    ylabel(axHist1A, 'Fraction of total population')
    xlabel(axHist1A, 'Attacker')
    axHist1A.FontSize = 10;
    axHist1A.XTick = [];

    plotPopFracBars(interpy(i,:),noConts,noHitBins,axHist1)
    axHist1.YColor = [1,1,1];
    axHist1.YTick = [];
    axHist1.Box = 'off';
    axHist1.FontSize = 10;

    frame = getframe(gcf);
    writeVideo(writer,frame)

    cla
end

close(writer)