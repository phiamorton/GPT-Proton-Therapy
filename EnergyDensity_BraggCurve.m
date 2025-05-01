% Read the data from the text file using readtable
phioffsets = [0.00]; %[0  0.33        0.66        0.99        1.32        1.65        1.98        2.31        2.65        2.98      3.14   3.31        3.64        3.97         4.30        4.63        4.96        5.29        5.62        5.95        6.28]; %3.4; %in rad, 0-2pi
 %[0.33]; %3.4; %in rad, 0-2pi
%[0  0.33        0.66        0.99        1.32        1.65        1.98        2.31        2.65        2.98        3.31        3.64        3.97         4.3        4.63        4.96        5.29        5.62        5.95        6.28]; %3.4; %in rad, 0-2pi
energyspreadpercent= 0.03; % in %
energy0=228.5; %MeV
phioffsetE=phioffsets;
masterfilename = sprintf('output_EnergyMod_phi%.2f_E%.2f_Esp%.2f', phioffsetE, energy0, energyspreadpercent);
%%Material and constants
I = 75 *10^(-6); % MeV or 80.8+-0.3
rho = 1; % g/cc for water
TZ = 10; % target Z
TA = 18; % target A
NA = 6.023e23; % Avogadro's number
re = 2.817e-13; % classical electron radius in centimeters
mec2 = 0.511; % MeV (rest mass energy of electron)
Z = 1; % proton
q = 1; % charge of proton
e_0 = 931.5; % MeV
A = 1; % A for proton=1
n= 3.34*10^23; %electron density of water in 1/cc
const = 4*pi*n*re^2*mec2*q^2;  %MeV/cm %4 * pi * NA * re^2 * mec2  * Z * q^2 / TA %MeV/cm


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
%translate simulated particles into real # of particles 
    Qtot0 = 4.2e-15; %from phia_test_Emod_spreadBragg.m, going to assume this is in C
    Qproton =1.6e-19; %C
    numrealprotons= Qtot0/Qproton; %this is total number of real protons
    numsimpart=2000; %from phia_test_EMod_spreadBragg.m
    sim_particles_scaling=numrealprotons/numsimpart; %converting between simulated and real particles
    
%%Functions and arrays
% Function to calculate beta^2
function beta2 = calcbeta2(E, A, e_0)
    beta2 = 1 - (e_0^2 / (e_0 + E/A)^2);
    % Ensure beta^2 is valid
    if beta2 < 0
        beta2 = 0; % Set to 0 if beta^2 is negative
    end
end


% Function to calculate W_max (maximum energy transfer)
function W_max = calcW_max(E, A, mec2, e_0)
    beta_val2 = calcbeta2(E, A, e_0);
    % Prevent invalid values for W_max
    if beta_val2 >= 1
        W_max = 0; % Set to 0 if invalid
    else
        W_max = 2 * (mec2) * beta_val2 / (1 - beta_val2); % MeV
    end
end

% Function to calculate stopping power (dE/dx)
function dEdX = calc_stoppingpower(E_val, A, I, e_0, rho, const, Z, mec2, TA)
    W_max_new = calcW_max(E_val, A, mec2, e_0);
    if W_max_new <= I
        dEdX = 0; % Set to 0 if W_max is too small or invalid
    else
        beta2 = calcbeta2(E_val, A, e_0);
        % Stopping power formula with an increased dEdX near Bragg peak
        dEdX = const / (beta2) * (log(W_max_new /(I)) - beta2); %MeV/cm
    end
end

% Material and simulation parameters
material_length = 80; % cm 
numsteps = 10000; % Increase the number of steps for better resolution
dx = material_length / numsteps;
x_values = linspace(0, material_length, numsteps); %cm
%plotbrowser

%%the no RF case
comparison_noE = readtable(sprintf('output_noRF_EnergyMod_phi%.2f_E%.2f_Esp%.2f_uniform.txt', phioffsetE, energy0, energyspreadpercent)); %sprintf('comparison_phi0.00_E228_Esp0.03.txt'));
G_comp =comparison_noE.G;
G_comp=G_comp(~isnan(G_comp));
meandEdX_comp=zeros(numsteps);
dose_comp=zeros(numsteps);
npart_comp=length(G_comp);
%length(G_comp)
for proton = 1:length(G_comp)
    E_0 =938.272*(G_comp(proton)-1); %MeV, starting energy
    % Loop to calculate energy loss and stopping power
    E_new=E_0;
    if isnan(E_new)==1
        break;
    end
    for ii = 2:numsteps
        if E_new <= 0 
            break;
        end
        if isnan(E_new)==0 && ii==2
            E_new;
        end
        % If energy becomes too low, break out of the loop

        % Calculate stopping power at the current energy
        dEdX_comp = calc_stoppingpower(E_new, A, I, e_0, rho, const, Z, mec2, TA);

        % Energy loss at each step
        E_loss = dEdX_comp * dx; %MeV

        % Update proton's energy (E decreases)
        E_new = E_new - E_loss;

        % Store the new stopping power value at the current step
        dEdX_new_comp=calc_stoppingpower(E_new, A, I, e_0, rho, const, Z, mec2, TA);
        %mass_stop_power_nominal = 1/rho * meandEdX_comp %MeVcm^2/g
        %dEdX_values(proton, ii) = dEdX_new;

        meandEdX_comp(ii) = meandEdX_comp(ii)+ dEdX_new_comp;
        dose_comp(ii) = dose_comp(ii) +1/rho*dEdX_new_comp*npart_comp/A_beam; %MeV/cm *cm^3/g *1/cm^2 = MeV/g, 
        %rho is g/cm^3, 1g/cc for water, A is beam lateral size in cm^2

    end
end

%%the RF case for each phase
max_Egain=zeros(1, length(phioffsets));

for pp = 1:length(phioffsets)
    phase = phioffsets(pp);
    datafile=sprintf('%s.txt',masterfilename);
    data = readtable(datafile);

    %% Extract the columns from the table
    G = data.G;
    %G= phiasimulationsEnergyModphi0.VarName6;
    G=G(~isnan(G));
    max_Egain(pp)= 938.272*(max(G)-1);
    dEdX_values = zeros(length(G), numsteps);
    dose_vals = zeros(numsteps);
    meandEdX=zeros(numsteps);
    size(dEdX_values);
    npart0=length(G);
    for proton = 1:length(G)
        E_0 =938.272*(G(proton)-1); %MeV, starting energy
        % Loop to calculate energy loss and stopping power
        E_new=E_0;
        if isnan(E_new)==1
            break;
        end
        for ii = 2:numsteps
            if E_new <= 0 
                break;
            end
            if isnan(E_new)==0 && ii==2
                E_new;
            end
            % If energy becomes too low, break out of the loop

            % Calculate stopping power at the current energy
            dEdX = calc_stoppingpower(E_new, A, I, e_0, rho, const, Z, mec2, TA);

            % Energy loss at each step
            E_loss = dEdX * dx; %MeV

            % Update proton's energy (E decreases)
            E_new = E_new - E_loss;

            % Store the new stopping power value at the current step
            dEdX_new=calc_stoppingpower(E_new, A, I, e_0, rho, const, Z, mec2, TA); %MeV/cm
            %dEdX_values(proton, ii) = dEdX_new;

            meandEdX(ii) = meandEdX(ii)+ dEdX_new; 
            dose_vals(ii)= dose_vals(ii) + 1/rho*dEdX_new*npart0/A_beam; 


        end
    end
    %meandEdX=meandEdX/length(G);
    dose_comp= dose_comp * 1.602e-10;  %MeV/g to J/kg [Gy]
    dose_vals= dose_vals * 1.602e-10 ; %MeV/g to J/kg [Gy]
    %% Create a plot
    %figure(gcf)
    figure('Visible','off')
    %scale up from sim particles to real particles
    scaled_dose_comp= dose_comp * sim_particles_scaling;  
    scaled_dose_vals= dose_vals * sim_particles_scaling; 

    subplot(2,1,1)
    %plot(x_values(1:numsteps), meandEdX(1:numsteps), 'b', 'LineWidth', 2,'Color',"#0072BD", 'DisplayName', 'With RF, 0.03% energy spread');
    plot(x_values(1:numsteps), scaled_dose_vals(1:numsteps), 'b', 'LineWidth', 2,'Color',"#0072BD", 'DisplayName', sprintf('With RF, E=%.2f MeV %.2f%% energy spread',energy0, energyspreadpercent));
    hold on;
    %plot(x_values(1:numsteps), meandEdX_comp(1:numsteps), 'b', 'LineWidth', 2,'LineStyle',':', 'Color',"#7E2F8E",'DisplayName','No RF');
    plot(x_values(1:numsteps), scaled_dose_comp(1:numsteps), 'b', 'LineWidth', 2,'LineStyle',':', 'Color',"#7E2F8E",'DisplayName','No RF, 0.03% energy spread');
    
    hold off;
   
    xlim([29,35]);
    xlabel('Position [cm]', 'FontSize', 12);
    %ylabel('Stopping Power (dE/dx) [MeV/cm]', 'FontSize', 12);
    ylabel('Dose [Gy]', 'FontSize', 12);
    legend('Location','northwest')
    title(sprintf('Bragg Curve for Proton in Water, phase offset = %.2f rad', phase), 'FontSize', 14);
    grid on;
    % title(sprintf('phase offset is %.2f radians', phase));
    % grid on;
    %figure('WindowStyle','docked', 'Name', sprintf('Phase %.2f', phase), 'NumberTitle', 'off')
    subplot(2,1,2)
    E_spec=938.272*(G-1);
    E_spec_comp=938.272*(G_comp-1);
    histogram(E_spec_comp, 100, 'FaceColor',"#7E2F8E", 'DisplayName','No RF, 0.03% energy spread');
    hold on;
    hh=histogram(E_spec, 100, 'FaceColor',"#0072BD", 'DisplayName', sprintf('With RF, E=%.2f MeV %.2f%% energy spread',energy0, energyspreadpercent));
    xlabel('Energy [MeV]');
    ylabel('Simulated Particles');
    bincounts= hh.BinCounts;
    hold on;
    ylim([0,max(bincounts)])
    yyaxis right 
    ylabel('Particles')
    bincounts_scaled=bincounts*sim_particles_scaling;
    ylim([0,max(bincounts_scaled)])
    title(sprintf('phase offset is %.2f radians', phase));
    legend('Location','northwest')
    saveas(gcf,sprintf('%s.png',masterfilename));
    %shg
end

% figure
% plot(phioffsets(1:length(phioffsets)),max_Egain(1:length(phioffsets)), 'b', 'LineWidth', 2)
% xlabel('Phase')
% ylabel('Max energy reached [MeV]')
% xlim([0,6.28])
% 
max_Egain
max_E_diff = -(energy0 - max_Egain) %MeV
length_cell = 0.0236;  %m
% %use with input 0 energy spread
MaxGradSeen = max_E_diff/(2*length_cell); %MV/m %2 cells %peak surface E field in this case is 68 MV/m
% Compute the average energy gain
AvgEGain = (mean(E_0 - energy0))*10^6;
AvgGrad=AvgEGain/(2*length_cell);
shuntimpedance= 54.8e6 %Ohm/m
fprintf('Max gradient %.2f MV/m \n', MaxGradSeen)
%p_disp=max_E_diff^2/(2*shuntimpedance*2*length_cell)
p_disp=(MaxGradSeen*1000000)^2/(2*shuntimpedance); %W
fprintf('Power calculated from max gradient %.2f MW \n', p_disp/1e6)
fprintf('Power calculated from max gradient %.2f MW per cavity \n', p_disp/1e6/2)
power_avg= AvgGrad^2/(2*shuntimpedance); %MV^2/m*MOhm/m = V^2/Ohm= W 
%power loss of klystron/power into cavities ?
fprintf('Power calculated from average gradient %.2f MW \n', power_avg/1e6)
fprintf('Power calculated from average gradient %.2f MW per cavity \n', power_avg/1e6/2)
%disp(power_avg)
%power_max = max_E_diff^2/(2*length_cell*shuntimpedance)  %MV^2/m*MOhm/m = V^2/Ohm= W %power loss of klystron/power into cavities
%https://cds.cern.ch/record/1005047/files/p145.pdf 
% %Cavity run with 2.5 MW into each cell
% cell_Power=2.5; %MW
% RF_power = 400; %kW
% 
% %Shunt_Imp = MaxGradSeen^2*length_cell/cell_Power %MeV^2/(m* MW)
% %https://accelconf.web.cern.ch/l06/papers/TUP044.pdf
Avg_grad_calc = 15*sqrt(power_avg/100000);  %MV/m, assumes P is in W
fprintf('Average gradient calculated from power %.2f MV/m \n', Avg_grad_calc)
%disp(AvgGrad-Avg_grad_calc)
% MaxGradSeen/Avg_grad; %2.26

