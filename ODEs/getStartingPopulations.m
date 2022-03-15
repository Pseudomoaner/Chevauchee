function startPops = getStartingPopulations(J,popSize,noBins)

probs = binopdf(0:4,4,J)*popSize; %Assume all of these different contact populations are in the first hit bin (zero hits)
listRec = [probs;zeros(noBins-1,size(probs,2))];
startPops = listRec(:);