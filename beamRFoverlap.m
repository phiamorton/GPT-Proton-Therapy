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
%zRF_sorted=sort(zRF)
lengthRF=length(zRF);
k = ceil(lengthRF * overlap);

% --- Step 2: Find indices of top 'k' values in zRF
[~, sorted_idx] = sort(zRF, 'descend');
top_idx = sorted_idx(1:k);

% --- Step 3: Use these indices to grab corresponding data
x_top = xRF(top_idx);
y_top = yRF(top_idx);
z_top = zRF(top_idx);
E_top = ERF(top_idx);

z_threshold = min(z_top) %where to splice the NoRF data

GNoRF = dataNoRF.G;
ENoRF=938.272*(GNoRF-1); %MeV
    %G=G(~isnan(G));
xNoRF=dataNoRF.x;
yNoRF=dataNoRF.y;
zNoRF=dataNoRF.z; 
lengthNoRF=length(zNoRF);

% --- Step 2: Filter zNoRF values >= threshold
idx_keep = zNoRF <= z_threshold;

% --- Step 3: Extract corresponding values
zNoRF_kept = zNoRF(idx_keep);
xNoRF_kept = xNoRF(idx_keep);
yNoRF_kept = yNoRF(idx_keep);
ENoRF_kept = ENoRF(idx_keep);
%length(ENoRF_kept)
%%plot 
figure;
hold on;
scatter(z_top, E_top, 'DisplayName','RF');

scatter(zNoRF_kept,ENoRF_kept, 'DisplayName', sprintf('No RF'));
hold off;
xlabel('z [cm]');
ylabel('Energy [MeV]');
legend();