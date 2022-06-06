clear all
close all

%Uses the outputs of diffusiveVariableVModel.m

OrCyCmap = [linspace(0.3,1,50)',linspace(0.6,0.5,50)',linspace(0.9,0.05,50)'];
load('C:\Users\olijm\Desktop\SeanAna\popTcourse_4.mat')
OD0pt001TC = squeeze(popsTcourse(:,:,1,:));
load('C:\Users\olijm\Desktop\SeanAna\popTcourse_3.mat')
OD0pt01TC = squeeze(popsTcourse(:,:,1,:));
load('C:\Users\olijm\Desktop\SeanAna\popTcourse_2.mat')
OD0pt1TC = squeeze(popsTcourse(:,:,1,:));
load('C:\Users\olijm\Desktop\SeanAna\popTcourse_1.mat')
OD1TC = squeeze(popsTcourse(:,:,1,:));

f = figure;
f.Units = 'Normalized';
f.Position = [0.1,0.1,0.6,0.25];

tiledlayout(1,4,'TileSpacing','none')
ax1 = nexttile;
ax2 = nexttile;
ax3 = nexttile;
ax4 = nexttile;

ax1.Box = 'on';
ax1.LineWidth = 1.5;
ax2.Box = 'on';
ax2.LineWidth = 1.5;
ax3.Box = 'on';
ax3.LineWidth = 1.5;
ax4.Box = 'on';
ax4.LineWidth = 1.5;

writer = VideoWriter('C:\Users\olijm\Desktop\SeanAna\diffExpComp.avi');
writer.FrameRate = 24;
open(writer)

for i = 1:433
    imagesc(ax1,OD1TC(:,:,i))
    imagesc(ax2,OD0pt1TC(:,:,i))
    imagesc(ax3,OD0pt01TC(:,:,i))
    imagesc(ax4,OD0pt001TC(:,:,i))

    ax1.XTick = [];
    ax1.YTick = [];
    ax2.XTick = [];
    ax2.YTick = [];
    ax3.XTick = [];
    ax3.YTick = [];
    ax4.XTick = [];
    ax4.YTick = [];

    axis(ax1,'equal')
    axis(ax1,'tight')
    colormap(ax1,OrCyCmap)
    title(ax1,'OD 1')
    caxis(ax1,[0,1])
    axis(ax2,'equal')
    axis(ax2,'tight')
    colormap(ax2,OrCyCmap)
    title(ax2,'OD 0.1')
    caxis(ax2,[0,1])
    axis(ax3,'equal')
    axis(ax3,'tight')
    colormap(ax3,OrCyCmap)
    title(ax3,'OD 0.01')
    caxis(ax3,[0,1])
    axis(ax4,'equal')
    axis(ax4,'tight')
    colormap(ax4,OrCyCmap)
    title(ax4,'OD 0.001')
    caxis(ax4,[0,1])

    cb = colorbar(ax4);
    cb.Label.String = 'Local attacker fraction';
    cb.FontSize = 12;

    text(ax1,5,5,sprintf('t = %ih %im', floor((i-1)/18), round(rem((i-1),18)*3.33333)),'FontSize',12)

    frame = getframe(gcf);
    writeVideo(writer,frame)
end

close(writer)