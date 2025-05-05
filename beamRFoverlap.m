phase = 0.00;  % Define phase value for grabbing the photo
phioffsetE=phase
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

%combine
z_combined=[zNoRF_kept;z_top];
E_combined=[ENoRF_kept;E_top];
y_combined=[yNoRF_kept;y_top]*100;
x_combined=[xNoRF_kept;x_top]*100;


%%plot 
% figure;
% hold on;
% histogram([ENoRF_kept;E_top], 100, 'FaceColor', "#7E2F8E", 'DisplayName', sprintf('combined')); %concatenate
% 
% histogram(E_top, 100, 'FaceColor', "#0072BD", 'DisplayName','RF');
% 
% histogram(ENoRF_kept, 50, 'FaceColor', "#77AC30", 'DisplayName', sprintf('No RF'));
% 
% hold off;
% xlabel('Energy [Mev]');
% ylabel('counts [simulated particles]');
% legend();
% 
% figure;
% hold on;
% %scatter([zNoRF_kept;z_top],[ENoRF_kept;E_top], 'Color', "#7E2F8E", 'DisplayName', sprintf('combined')); %concatenate
% 
% scatter(z_top, E_top,'Color', "#0072BD", 'DisplayName','RF');
% 
% scatter(zNoRF_kept,ENoRF_kept, 'Color', "#0072BD", 'DisplayName', sprintf('No RF'));
% 
% hold off;
% xlabel('z [cm]');
% ylabel('Energy [Mev]');
% legend();
% 
% 
% figure;
% hold on;
% %scatter([zNoRF_kept;z_top],[ENoRF_kept;E_top], 'Color', "#7E2F8E", 'DisplayName', sprintf('combined')); %concatenate
% 
% %scatter(y_top, E_top, 'DisplayName','RF','Color',"red",'Marker','*');
% %scatter(yNoRF_kept,ENoRF_kept, 'DisplayName', sprintf('No RF'),'Color',"red",'Marker','+');
% 
% scatter([yNoRF_kept;y_top],[ENoRF_kept;E_top])
% hold off;
% xlabel('y [cm]');
% ylabel('Energy [Mev]');
% legend();
%need to feed this into making an image to do reconstruction with

%%make the output image

newWidth = 1292;
newHeight = 964;
%real camera calibration ruler photo 1292x964 pixels, 18.5cm-6.8cm
realwidth=18.5-6.8; %in width
realpixelsize=realwidth/newWidth;
realyrange=realpixelsize*newHeight;

mevion_25nA=false;
mevion_1nA=true;

if mevion_1nA==true
    xrms0 = 3.495/1000 ;%m 4.9mm
    %based on mevion numbers this will be ~3-4mm at 1nA or 5-6mm at 25 nA
    yrms0 = 4.007/1000; %m 6mm
    %divergence of beam
    divangx0 = (3.794-3.496)/120; %.58; %change in x [mm] over 12 cm
    divangy0 = (4.299-4.007)/120; % .67;
end  
if mevion_25nA==true
    xrms0 = 4.906/1000 ;%m 4.9mm
    %based on mevion numbers this will be ~3-4mm at 1nA or 5-6mm at 25 nA
    yrms0 = 6.039/1000; %m 6mm
    %divergence of beam
    divangx0 = (5.198-4.906)/120; %.58; %change in x [mm] over 12 cm
    divangy0 = (6.289-6.039)/120; % .67;
end 
yrms0=yrms0*100; %cm
xrms0=xrms0*100; %cm
A_beam=pi*xrms0*yrms0; %cm^2
rho=1;

num_particles = length(z_combined);

% --- Setup dose grid ---
z_max = 50;                  
num_z = newWidth;                
z_vals = linspace(0, z_max, num_z); 
dx=z_max/num_z;

% Increase y_bins range if needed
realyrange=8.7297;
y_bins = linspace(-realyrange/2, realyrange/2, newHeight);  % Wide enough to show thin beam
dose_map = zeros(num_z, length(y_bins));    % z vs y (dose map)

% --- Accumulate dose per particle ---
for i = 1:num_particles
    % Compute dose profile for this particle
    dose_profile = computeDoseBetheBloch(E_combined(i),dx,num_z,rho,A_beam);
    
    % Find the nearest y-bin index for this particle's y position
    [~, y_idx] = min(abs(y_bins - y_combined(i)));
    
    % Accumulate the dose for the corresponding y-bin
    dose_map(:, y_idx) = dose_map(:, y_idx) + dose_profile;
    i;
end

% --- Plot result ---
figure;
dose_map;
imagesc(z_vals,y_bins, transpose(dose_map));
% xlabel('Z Position (cm)');
% ylabel('Y Position (cm)');
% colorbar;


ylabel('y [cm]');
xlabel('depth [cm]');
xmin=28;
xlim([xmin,xmin+11.7]);
ylim([-realyrange/2,realyrange/2]); 
a=colorbar;
a.Label.String = 'Dose [a.u.]';
title(sprintf('Bragg Curve for Proton in Water, phase offset = %.2f rad, energy=%.2f', phioffsetE,energy0), 'FontSize', 14);
%set(gca,'XTick',[], 'YTick', [],'Visible', 'off', 'Color','none');
saveas(gcf,sprintf('%sBraggIm_unknownoverlap.png', masterfilename))

%%Functions
function dose = computeDoseBetheBloch(E0, dx,num_z,rho,A_beam)
    dose_vals = zeros(num_z(:));
    E_new=E0;
    for ii = 2:num_z  
        dEdX = calc_stoppingpower(E_new);
        E_loss = dEdX * dx;
        E_new = E_new - E_loss;
        dEdX_new=calc_stoppingpower(E_new);
        dose_vals(ii)= 1/rho*dEdX_new*1/A_beam;
    end
    dose=transpose(dose_vals(1:num_z));
    size(dose);
end


% Function to calculate beta^2
function beta2 = calcbeta2(E, A, e_0)
    beta2 = 1 - (e_0^2 / (e_0 + E/A)^2);
    % Ensure beta^2 is valid
    if beta2 < 0
        beta2 = 0; % Set to 0 if beta^2 is negative
    end
end
 
% % Function to calculate W_max (maximum energy transfer)
function W_max = calcW_max(E, A, mec2, e_0)
    beta_val2 = calcbeta2(E, A, e_0);
    % Prevent invalid values for W_max
    if beta_val2 >= 1
        W_max = 0; % Set to 0 if invalid
    else
        W_max = 2 * (mec2) * beta_val2 / (1 - beta_val2); % MeV
    end
end

% % Function to calculate stopping power (dE/dx)
function dEdX = calc_stoppingpower(E_val)
    I = 75 *10^(-6); % MeV or 80.8+-0.3
    re = 2.817e-13; % classical electron radius in centimeters
    mec2 = 0.511; % MeV (rest mass energy of electron)
    q = 1; % charge of proton
    e_0 = 931.5; % MeV
    A = 1; % A for proton=1
    n= 3.34*10^23; %electron density of water in 1/cc
    const = 4*pi*n*re^2*mec2*q^2;
    W_max_new = calcW_max(E_val, A, mec2, e_0);
    if W_max_new <= I
        dEdX = 0; % Set to 0 if W_max is too small or invalid
    else
        beta2 = calcbeta2(E_val, A, e_0);
        % Stopping power formula with an increased dEdX near Bragg peak
        dEdX = const / (beta2) * (log(W_max_new /(I)) - beta2);
    end
end