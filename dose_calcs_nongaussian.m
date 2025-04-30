
phioffsets = [0.00]; %[0  0.33        0.66        0.99        1.32        1.65        1.98        2.31        2.65        2.98      3.14   3.31        3.64        3.97         4.30        4.63        4.96        5.29        5.62        5.95        6.28]; %3.4; %in rad, 0-2pi
energyspreadpercent= 0.03;
energy0=228.5; %MeV
phioffsetE=phioffsets;
masterfilename = sprintf('output_EnergyMod_phi%.2f_E%.2f_Esp%.2f', phioffsetE, energy0, energyspreadpercent);
data = readtable(sprintf('%s.txt',masterfilename));
rho=1;

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

%% Extract the columns from the table
G = data.G(1:100);
x_beam=data.x(1:100)*100;
y_beam=data.y(1:100)*100;
E0 =938.272*(G-1); %MeV, starting energy
num_particles = length(G);

% --- Setup dose grid ---
z_max = 50;                  
num_z = 1000;                
z_vals = linspace(0, z_max, num_z); 
dx=z_max/num_z

% Increase y_bins range if needed
realyrange=8.7297;
y_bins = linspace(-realyrange/2, realyrange/2, 200);  % Wide enough to show thin beam
dose_map = zeros(num_z, length(y_bins));    % z vs y (dose map)

% --- Accumulate dose per particle ---
for i = 1:num_particles
    % Compute dose profile for this particle
    dose_profile = computeDoseBetheBloch(E0(i),dx,num_z,rho,A_beam);
    
    % Find the nearest y-bin index for this particle's y position
    [~, y_idx] = min(abs(y_bins - y_beam(i)));
    
    % Accumulate the dose for the corresponding y-bin
    dose_map(:, y_idx) = dose_map(:, y_idx) + dose_profile;
    i
end

% --- Plot result ---
figure;
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
saveas(gcf,sprintf('%sBraggIm_nongaussian.png', masterfilename))

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