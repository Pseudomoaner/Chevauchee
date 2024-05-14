load('C:\Users\omeacock\Downloads\SprResults_f1_0.000000_f2_1.000000_fire1_0.020000_fire2_0.020000.mat')

zetas = 0:0.001:0.1;
CAs = zeros(size(zetas,2));

sampTime = 100;

noAt = sum(trackableData.Population{1} == 't');
noSens = sum(trackableData.Population{1} == 's');

for zetIndT = 1:size(zetas,2)
    zetT = zetas(zetIndT); %Potency of the 't' population's toxin (motile pop.)
    for zetIndS = 1:size(zetas,2)
        zetS = zetas(zetIndS); %Potency of the 's' population's toxin (non-motile pop.)

        finThits = trackableData.Hit{sampTime}(trackableData.Population{sampTime}=='t');
        finShits = trackableData.Hit{sampTime}(trackableData.Population{sampTime}=='s');

        finTimpact = finThits*zetS; %Final impact on the 't' population from the 's' population
        finTimpact(finTimpact>1) = 1;
        finSimpact = finShits*zetT;
        finSimpact(finSimpact>1) = 1;

        finSens = sum(1-finSimpact);
        finAt = sum(1-finTimpact);

        CAs(zetIndT,zetIndS) = (finAt/finSens)/(noAt/noSens);
    end
end

imagesc(zetas,zetas,log10(CAs))

rbcUp = [interp1(rbcmap(:,1),1:0.1:11)',interp1(rbcmap(:,2),1:0.1:11)',interp1(rbcmap(:,3),1:0.1:11)'];
colormap(rbcUp)

xlabel('Toxicity of strain 2')
ylabel('Toxicity of strain 1')

caxis([-1,1])

ax = gca;
ax.YDir = 'normal';
ax.Box = 'on';
ax.LineWidth = 1.5;

cb = colorbar;
cb.Label.String = 'log_{10}(Competitive advantage of strain 1)';
cb.Label.FontSize = 12;
cb.LineWidth = 1.5;

title('Strain 1 mobile, strain 2 static')