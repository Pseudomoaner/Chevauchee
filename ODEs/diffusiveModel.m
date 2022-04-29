clear all
close all

v = 0;
lam = 0.01;
atFrac = 1/10;

dx = 5; %Granularity of coarse-grained lattice
xWidth = 200;
yHeight = 200;
noX = xWidth/dx;
noY = yHeight/dx;

tMax = 1000;
diffDt = 5; %Timestep between diffusion timesteps
noDiffTstesp = tMax/diffDt;
D = v*4; %4 is an approximate scaling factor between the single-cell velocity and the diffusion constant

noHitBins = 6; %Number of different hit categories to take into consideration, from 0 to noHitBins-1 plus
noConts = 5; %Number of contacts made by each cell 

%Initial conditions
[startA,startS] = initialisePatchyField(dx,xWidth,yHeight,0.016,atFrac);
pops = cat(3,startA,startS,zeros(noX,noY,noHitBins*(noConts+1)-1)); %Create population array, attackers in first layer, unhit sensitives in second)

%Allow the system to equilibrate contact compartments without killing or
%diffusion
[t,pops] = ode45(@(t,y)diffusiveODEs(t,y,1,0,noX,noY,noConts,noHitBins),[0,100],pops(:));
pops = pops(end,:);
pops = reshape(pops,noY,noX,noHitBins*(noConts+1) + 1);
popsTcourse = pops;

%Outer loop - macroscopic mixing (diffusion)
for t = 1:noDiffTstesp
    %Run diffusion of each of the populations separately (can we do this?
    %Will they be guaranteed to add up to one everywhere if they're run
    %separately?)
    for i = 1:size(pops,3)
        pops(:,:,i) = diffTimestepCN(pops(:,:,i),diffDt,dx,D,true);
    end

    %Inner loop - microscopic mixing (contact swapping)
    [t,pops] = ode45(@(t,y)diffusiveODEs(t,y,v,lam,noX,noY,noConts,noHitBins),[0,diffDt],pops(:));
    pops = pops(end,:);
    pops = reshape(pops,noY,noX,noHitBins*(noConts+1) + 1);
    popsTcourse = cat(4,popsTcourse,pops);

    imagesc(1-pops(:,:,2)-pops(:,:,1))
    caxis([0,1])
    axis equal
    pause(0.1)
end