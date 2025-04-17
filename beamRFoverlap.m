phase = 0.00;  % Define phase value for grabbing the photo
energyspreadpercent= 0.03
energy0=228.5 %MeV
uniform=true
a = .005; %0.5 cm
overlap=0.5; %amount of overlap with RF pulse, 0= no overlap, 1=100% overlap
%%files with same parameters but RF=0 in the reference, must already be
%%created, use NoRF=true in phia_test_EMod_spreadBragg.m and run
%%sim_auto_uniform_noRF.bat with correct E, Esp, phase params
masterfilename = sprintf('output_EnergyMod_phi%.2f_E%.2f_Esp%.2f_uniform', phioffsetE, energy0, energyspreadpercent);
reference_noRFfilename= sprintf('output_noRF_EnergyMod_phi%.2f_E%.2f_Esp%.2f_uniform', phioffsetE, energy0, energyspreadpercent);

%truncate that the overlap mark and add the data 

dataRF = readtable(sprintf('%s.txt',masterfilename));
dataNoRF = readtable(sprintf('%s.txt',reference_noRFfilename));


    %% Extract the columns from the table
GRF = dataRF.G;
ERF=938.272*(GRF-1); %MeV
    %G=G(~isnan(G));
xRF=dataRF.x;
yRF=dataRF.y;
zRF=dataRF.z; 
lengthRF=length(zRF);

GNoRF = dataNoRF.G;
ENoRF=938.272*(GNoRF-1); %MeV
    %G=G(~isnan(G));
xNoRF=dataNoRF.x;
yNoRF=dataNoRF.y;
zNoRF=dataNoRF.z; 
lengthNoRF=length(zNoRF)

    