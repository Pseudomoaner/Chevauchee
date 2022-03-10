function [] = writePopFracBarMovie(y,t)

interpy = interp1(t,y,2:2:max(t),'pchip');

figure()
fig = gcf;
fig.Units = 'Normalized';
fig.Position=[0.2,0.2,0.3,0.2];
ax = gca;

writer = VideoWriter('C:\Users\olijm\Desktop\SeanAna\test.avi');
writer.FrameRate = 60;
open(writer)

for i = 1:size(interpy,1)
    plotPopFracBars(interpy(i,:),ax)
    frame = getframe(gcf);
    writeVideo(writer,frame)

    cla
end

close(writer)