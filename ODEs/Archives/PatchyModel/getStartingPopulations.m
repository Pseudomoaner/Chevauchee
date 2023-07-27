function startPops = getStartingPopulations(J,popSize,noHitBins,noConts)

probs = binopdf(0:noConts,noConts,J)*(popSize-J); %Assume all of these different contact populations are in the first hit bin (zero hits)
listRec = [probs;zeros(noHitBins-1,size(probs,2))];
startPops = listRec(:);