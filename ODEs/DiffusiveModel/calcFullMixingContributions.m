function [base,contEx,homo] = calcFullMixingContributions(dx,startA,startS,atFrac,noConts,noHitBins,noX,noY,lam,vs,diffDt,tMax,hitEfficiency)
%CALCFULLMIXINGCONTRIBUTIONS calculates the contributions of target
%switching and homogenisation to the intoxification dynamics of a continuum
%simulation with the specified settings.
%
%   INPUTS:
%       - dx: spacing of continuum lattice
%       - startA: Initial distribution of the attacker population
%       - startS: Initial distribution of the sensitive population
%       - atFrac: Overall fraction of attacker population
%       - noConts: Number of contact bins to be included in simulation.
%       Typically 5
%       - noHitBins: Number of bins for labelling different numbers of
%       hits.
%       - noX, noY: Size of simulation lattice
%       - lam: Intoxification rate
%       - vs: List of velocities, used to define the mixing rates
%       - diffDt: Timestep size between diffusion steps
%       - tMax: Maximal simulation timepoint
%       - hitEfficiency: Fractional decrease in the growth rate of
%       sensitive cell upon application of a contact-dependent hit.
%
%   OUTPUTS:
%       - base: The reduction in the growth of the sensitive
%       population with no motility compared to if there were no attackers.
%       Relative to 1-atFrac
%       - contEx: The reduction in growth of the sensitive population
%       relative to the base condition when contact exchange alone is active.
%       - homo: The reduction in the growth of the sensitive population
%       relative to the contact-exchange-only condition when both contact
%       exchange and genotypic mixing are active.
%
%   Author: Oliver J. Meacock, 2023

alphaD = 5.638; %Proportionality constant that converts velocity into cell diffusion rate
noDiffTsteps = tMax/diffDt + 1;

%Initial conditions
pops = cat(3,startA,startS,zeros(noY,noX,noHitBins*(noConts+1)-1)); %Create population array, attackers in first layer, unhit sensitives in second)

%Allow the system to equilibrate contact compartments without killing or
%diffusion
[t,pops] = ode45(@(t,y)diffusiveODEs(t,y,1,0,noX,noY,noConts,noHitBins),[0,100],pops(:));
pops = pops(end,:);
popsInit = reshape(pops,noY,noX,noHitBins*(noConts+1) + 1);

%% Part 1: No mixing
[~,popsNM] = ode45(@(t,y)diffusiveODEs(t,y,0,lam,noX,noY,noConts,noHitBins),[0,tMax],popsInit(:));
popsNM = popsNM(end,:);
popsNM = reshape(popsNM,noY,noX,noHitBins*(noConts+1) + 1);

%% Part 2: Contact exchange only
popsCE = popsInit;
for tS = 1:noDiffTsteps
    [~,popsCE] = ode45(@(t,y)diffusiveODEs(t,y,vs(tS),lam,noX,noY,noConts,noHitBins),[0,diffDt],popsCE(:));
    popsCE = popsCE(end,:);
    popsCE = reshape(popsCE,noY,noX,noHitBins*(noConts+1) + 1);
end

%% Part 3: Both mixing types
popsB = popsInit;

for tS = 1:noDiffTsteps
    %Run diffusion of each of the populations separately
    D = alphaD * vs(tS);
    for i = 1:size(popsB,3)
        popsB(:,:,i) = diffTimestepCN(popsB(:,:,i),diffDt,dx,D,true);
    end

    %Inner loop - microscopic mixing (contact swapping)
    [~,popsB] = ode45(@(t,y)diffusiveODEs(t,y,vs(tS),lam,noX,noY,noConts,noHitBins),[0,diffDt],popsB(:));

    popsB = popsB(end,:);
    popsB = reshape(popsB,noY,noX,noHitBins*(noConts+1) + 1);
end
 
baseSens = 0;
contExSens = 0;
bothSens = 0;
for j = 0:noConts
    baseSens = baseSens + squeeze(sum(sum(sum(popsNM(:,:,2+j:noHitBins:end),1),2),3)) * max(0,(1-j*hitEfficiency));
    contExSens = contExSens + squeeze(sum(sum(sum(popsCE(:,:,2+j:noHitBins:end),1),2),3)) * max(0,(1-j*hitEfficiency));
    bothSens = bothSens + squeeze(sum(sum(sum(popsB(:,:,2+j:noHitBins:end),1),2),3)) * max(0,(1-j*hitEfficiency));
end

base = 1-atFrac-baseSens/(noX*noY);
contEx = 1-atFrac-contExSens/(noX*noY) - base;
homo = 1-atFrac-bothSens/(noX*noY) - base - contEx;