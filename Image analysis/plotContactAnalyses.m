clear all
close all

Root = 'C:\Users\olijm\Desktop\SeanAna\Sample20x';
Branches = {'OD1','OD0_1','OD0_01','OD0_001'};
fileName = 'Analysis.mat';

colours = [1,0.5,0;0.8,0.65,0;0.6,0.8,0;0.4,1,0];
alignFrames = [7,16,22,28];
dt = 0.5;

figure(1)
subplot(2,1,1)
hold on
ax1 = gca;

subplot(2,1,2)
hold on
ax2 = gca;

for b = 1:size(Branches,2)
    load(fullfile(Root,Branches{b},fileName),'pKs','packFracs')
    
    times = (0:dt:(size(packFracs,1)-1)*dt) - (alignFrames(b)-1)*dt;
    plot(ax1,times,packFracs,'Color',colours(b,:))
    plot(ax2,times,pKs,'Color',colours(b,:))
end

ax1.Box = 'on';
ax1.LineWidth = 1.5;
ylabel(ax1,'Packing fraction')
lgd = legend('OD 1','OD 0.1','OD 0.01','OD 0.001');

ax2.Box = 'on';
ax2.LineWidth = 1.5;
axis(ax2)
ylabel(ax2,'p(K)')
xlabel(ax2,'time relative to confluency (hr)')