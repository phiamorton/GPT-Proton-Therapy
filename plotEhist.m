% Read the data from the text file using readtable
data = readtable('phia_run_energyMod_spreadBragg.txt');
data90 = readtable('phia_simulationsEnergyMod_phi1.57hist.txt');
data180 = readtable('phia_run_energyMod_spreadBragg_pihist.txt');
data270 = readtable('phia_simulationsEnergyMod_phi2.35hist.txt');
data360 = readtable('phia_simulationsEnergyMod_phi6.28hist.txt');
% 
% % Extract the columns from the table
G = data.G;
G90 = data90.G;
G180 = data180.G;
G270 = data270.G;
G360 = data360.G;

E=938.272*(G-1); %MeV
E90=938.272*(G90-1); %MeV
E180=938.272*(G180-1); %MeV
E270=938.272*(G270-1); %MeV
E360=938.272*(G360-1); %MeV
% 
% % Create a plot
figure
histogram(E)
%ylim([0,12])
xlabel('Energy [MeV]');
title('No Phase shift');
legend show;
% grid on;

figure
histogram(E90)
%ylim([0,12])
title('\pi/2 phase shift');
xlabel('Energy [MeV]');
legend show;


figure
histogram(E180)
%ylim([0,12])
title('\pi phase shift');
xlabel('Energy [MeV]');
legend show;


figure
histogram(E270)
%ylim([0,12])
title('3/4 \pi phase shift');
xlabel('Energy [MeV]');
legend show;


figure
histogram(E360)
%ylim([0,12])
title('2 \pi phase shift');
xlabel('Energy [MeV]');
legend show;