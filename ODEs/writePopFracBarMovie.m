function [] = writePopFracBarMovie(y,t,noConts,noHitBins)

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