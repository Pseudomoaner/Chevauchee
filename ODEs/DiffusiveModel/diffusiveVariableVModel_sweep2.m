clear all
close all

vmaxes = [0.0124,0.0207,0.0457,0.0630]/0.8;
confTs = [2,2,2,2]*3600;%[4.5,9.5,11.5,14]*3600;
rho0s = [1,0.1,0.01,0.001]*0.064;
colours = [10,58,92;31,126,193;76,178,250;156,212,252]/255;

lam = 0.003; %CDI firing rate (s^-1) - 0.001 corresponds to 3.6 firings / hr
hitEfficiency = 0.1; %The impact of each hit on the targeted cell's ongoing growth rate
atFrac = 5.5/11; %Attacker fraction
vrate_original = 3000; %Width of Gaussian velocity profile (s^-1)

%Dimensional parameters
dx = 10; %Granularity of coarse-grained lattice
xWidth = 600; %780 = 624 um / 0.8 um
yHeight = 600; %630 = 501 um / 0.8 um
noX = xWidth/dx;
noY = yHeight/dx;
w = 1; %Cell width
tMax = 4*3600;
diffDt = 200; %Timestep between diffusion timesteps
noDiffTsteps = tMax/diffDt;
tList = linspace(0,tMax,noDiffTsteps+1);

%Rate-related parameters
alphaD = 5.638; %Proportionality constant that converts velocity into cell diffusion rate
noHitBins = 11; %Number of different hit categories to take into consideration, from 0 to noHitBins-1 plus
noConts = 5; %Number of contacts made by each cell

data = zeros(1,16);

%hit thershold
hitThresh = 1;
hits = [1,2,3,4,5,6,7,8,9,10];

%ranges of values to check
ranges = [0.75,1,1.5]; %width of gaussian velocity profile
Vranges = [0.75,1.0,1.5]; %range of estimates for max velocity and confluency time

%firing rates
fires = [0.00125, 0.0025, 0.005, 0.01];

%firing rate loop
for fire = 1:size(fires,2)
lam = fires(fire);
%gaussian width loop
    for w = 1:size(ranges,2)
    vrate = vrate_original*ranges(w);
        



            for i = 1:4
                rho0 = rho0s(i); %Seeding density (cells cellWidth^-2)
                confTime = confTs(i); %Time to confluency
            for v = 1:size(Vranges,2)
            vmax = vmaxes(i)*Vranges(v); %Maximum velocity (cewllWidths s^-1)

                [startA,startS] = initialisePatchyField(dx,xWidth,yHeight,rho0,atFrac);
                pops = cat(3,startA,startS,zeros(noY,noX,noHitBins*(noConts+1)-1)); %Create population array, attackers in first layer, unhit sensitives in second)

                %Allow the system to equilibrate contact compartments without killing or
                %diffusion
                [~,pops] = ode45(@(t,y)diffusiveODEs(t,y,1,0,noX,noY,noConts,noHitBins),[0,3600],pops(:));
                pops = pops(end,:);
                pops = reshape(pops,noY,noX,noHitBins*(noConts+1) + 1);
                popsTcourse = pops;

                for t = 1:noDiffTsteps
                    currV = vmax*exp(-((t*diffDt-confTime)/vrate)^2);
                    currD = currV*alphaD;

                    %Run diffusion of each of the populations separately
                    for j = 1:size(pops,3)
                        pops(:,:,j) = diffTimestepCN(pops(:,:,j),diffDt,dx,currD,true);
                    end

                    %Inner loop - microscopic mixing (contact swapping)
                    [~,pops] = ode45(@(t,y)diffusiveODEs(t,y,currV,lam,noX,noY,noConts,noHitBins),[0,diffDt],pops(:));
                    pops = pops(end,:);
                    pops = reshape(pops,noY,noX,noHitBins*(noConts+1) + 1);
                    popsTcourse = cat(4,popsTcourse,pops);
                end

                atPopSize = sum(sum(popsTcourse(:,:,1,1),1),2);
                sensPopSize = sum(sum(sum(popsTcourse(:,:,2:end,1),1),2),3);

                J = squeeze(popsTcourse(:,:,1,:));
                val1 = max(squeeze(mean(mean(J.*(1-J),1),2))./(sensPopSize/(atPopSize+sensPopSize)));

                val2 = zeros(10,1);
                unhitTcourse = squeeze(sum(sum(sum(popsTcourse(:,:,2:noHitBins:end,:),1),2),3));
                liveSus = unhitTcourse;
                val2(1) = max(atPopSize./liveSus)./(atPopSize/sensPopSize);

                for threshloop = 2:10
                   liveSus = unhitTcourse + squeeze(sum(sum(sum(popsTcourse(:,:,(threshloop+1):noHitBins:end,:),1),2),3)) * max(0,(1-threshloop*hitEfficiency));
                   val2(threshloop) = max(atPopSize./liveSus)./(atPopSize/sensPopSize);
                end

                fprintf('firing rate %i\n', fires(fire))
                fprintf('v width %i\n', vrate)
                fprintf('density %i\n', rho0)
                fprintf('vmax %i\n', vmax)
                fprintf('confluence time %i\n', confTime)
                data = vertcat(data,[lam, rho0, vrate, vmax, confTime, val1, val2(1),val2(2),val2(3),val2(4), val2(5),val2(6),val2(7),val2(8),val2(9),val2(10)]);
                fprintf('next run\n\n')


            %     save(sprintf('C:\\Users\\olijm\\Desktop\\SeanAna\\popTcourse_%i.mat',i),'popsTcourse','vmax','rho0','confTime')
            end
        end
    end
end


writematrix(data,'new_sweep3.csv');