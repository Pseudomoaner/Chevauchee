clear all
close all

v = 0.1;
lam = 0.01;
atFrac = 1/10;
rho0 = 0.001;

dx = 10; %Granularity of coarse-grained lattice
xWidth = 300;
yHeight = 300;
noX = xWidth/dx;
noY = yHeight/dx;

tMax = 1000;
diffDt = 5; %Timestep between diffusion timesteps
noDiffTsteps = tMax/diffDt;
alphaD = 5.638; %Proportionality constant that converts velocity into cell diffusion rate
D = v*alphaD;

noHitBins = 6; %Number of different hit categories to take into consideration, from 0 to noHitBins-1 plus
noConts = 5; %Number of contacts made by each cell

plotting = false;

%Initial conditions
[startA,startS] = initialisePatchyField(dx,xWidth,yHeight,rho0,atFrac);
pops = cat(3,startA,startS,zeros(noY,noX,noHitBins*(noConts+1)-1)); %Create population array, attackers in first layer, unhit sensitives in second)

%Allow the system to equilibrate contact compartments without killing or
%diffusion
[t,pops] = ode45(@(t,y)diffusiveODEs(t,y,1,0,noX,noY,noConts,noHitBins),[0,100],pops(:));
pops = pops(end,:);
pops = reshape(pops,noY,noX,noHitBins*(noConts+1) + 1);
popsTcourse = pops;

if plotting
    f = figure;
    f.Units = 'Normalized';
    f.Position = [0.025,0.2,0.95,0.45];

    tiledlayout(1,12,'TileSpacing','compact')
    axHist1A = nexttile;
    axHist1 = nexttile([1,3]);
    axUH = nexttile([1,4]);
    axHist2A = nexttile;
    axHist2 = nexttile([1,3]);

    writer = VideoWriter('C:\Users\olijm\Desktop\SeanAna\diffTest.avi');
    writer.FrameRate = 24;
    open(writer)

    OrCyCmap = [linspace(0.3,1,50)',linspace(0.6,0.5,50)',linspace(0.9,0.05,50)'];
end

%Outer loop - macroscopic mixing (diffusion)
for t = 1:noDiffTsteps
    %Run diffusion of each of the populations separately
    for i = 1:size(pops,3)
        pops(:,:,i) = diffTimestepCN(pops(:,:,i),diffDt,dx,0,true);
    end

    %Inner loop - microscopic mixing (contact swapping)
    [t,pops] = ode45(@(t,y)diffusiveODEs(t,y,v,lam,noX,noY,noConts,noHitBins),[0,diffDt],pops(:));
    pops = pops(end,:);
    pops = reshape(pops,noY,noX,noHitBins*(noConts+1) + 1);
    popsTcourse = cat(4,popsTcourse,pops);

    if plotting
        imagesc(axUH,sum(pops(:,:,2:end),3)) %Shows all the unhit sensitives
        colormap(OrCyCmap)
        caxis(axUH,[0,1])
        axis(axUH,'equal')
        axis(axUH,'tight')
        axUH.Box = 'on';
        axUH.LineWidth = 1.5;
        axUH.XTick = [];
        axUH.YTick = [];
        cb = colorbar(axUH,'northOutside','LineWidth',1.5);
        cb.Label.String = 'Local fraction of sensitives';
        cb.FontSize = 10;

        bar(axHist1A,pops(6,5,1),'FaceColor',[0.2,0.6,1],'LineWidth',1.5)
        axHist1A.XColor = [0.75,0.25,0.75]; axHist1A.YColor = [0.75,0.25,0.75];
        axis(axHist1A,[0.5,1.5,0,1]);
        axHist1A.LineWidth = 1.5;
        axHist1A.Box = 'off';
        ylabel(axHist1A, 'Fraction of local population')
        xlabel(axHist1A, 'Attacker')
        axHist1A.FontSize = 10;

        plotPopFracBars(squeeze(pops(6,5,2:end)),noConts,noHitBins,axHist1)
        axHist1.XColor = [0.75,0.25,0.75]; axHist1.YColor = [1,1,1];
        axHist1.YTick = [];
        axHist1.Box = 'off';
        axHist1.FontSize = 10;

        bar(axHist2A,pops(12,13,1),'FaceColor',[0.2,0.6,1],'LineWidth',1.5)
        axHist2A.XColor = [0.25,0.75,0.75]; axHist2A.YColor = [0.25,0.75,0.75];
        axis(axHist2A,[0.5,1.5,0,1]);
        axHist2A.LineWidth = 1.5;
        axHist2A.Box = 'off';
        ylabel(axHist2A, 'Fraction of local population')
        xlabel(axHist2A, 'Attacker')
        axHist2A.FontSize = 10;

        plotPopFracBars(squeeze(pops(12,13,2:end)),noConts,noHitBins,axHist2)
        axHist2.XColor = [0.25,0.75,0.75]; axHist2.YColor = [1,1,1];
        axHist2.YTick = [];
        axHist2.Box = 'off';
        axHist2.FontSize = 10;

        rectangle(axUH,'Position',[4.5,5.5,1,1],'EdgeColor',[0.75,0.25,0.75],'LineWidth',1.5)
        rectangle(axUH,'Position',[12.5,11.5,1,1],'EdgeColor',[0.25,0.75,0.75],'LineWidth',1.5)

        pause(0.1)

        frame = getframe(gcf);
        writeVideo(writer,frame)

        cla(axUH)
    end
end

if plotting
    close(writer)
end
