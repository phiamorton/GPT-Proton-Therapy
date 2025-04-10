% Read the data from the text file using readtable
phioffsets = [0.00] %[0  0.33        0.66        0.99        1.32        1.65        1.98        2.31        2.65        2.98      3.14   3.31        3.64        3.97         4.30        4.63        4.96        5.29        5.62        5.95        6.28]; %3.4; %in rad, 0-2pi
 %[0.33]; %3.4; %in rad, 0-2pi
%[0  0.33        0.66        0.99        1.32        1.65        1.98        2.31        2.65        2.98        3.31        3.64        3.97         4.3        4.63        4.96        5.29        5.62        5.95        6.28]; %3.4; %in rad, 0-2pi

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
rad_beam = 0.5/2; %cm - 5mm beam diameter
A_beam= pi*rad_beam^2; %cm^2

%scaling between simlated and actual number of particles
Qtot0 = 4.2e-15; %from phia_test_Emod_spreadBragg.m, going to assume this is in C
Qproton =1.6e-19; %C
sim_particles_scaling= Qtot0/Qproton;
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
material_length = 80; % cm i hope
numsteps = 10000; % Increase the number of steps for better resolution
dx = material_length / numsteps;
x_values = linspace(0, material_length, numsteps); %cm
%plotbrowser

%the no RF case
comparison_noE = readtable(sprintf('phia_simulationsEnergyMod_phi0.00_0.03Espread_nominalhist.txt'));
G_comp =comparison_noE.G;
G_comp=G_comp(~isnan(G_comp));
meandEdX_comp=zeros(numsteps);
dose_comp=zeros(numsteps);
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
        dose_comp(ii) = dose_comp(ii) +1/rho*dEdX_new_comp*1/A_beam; %MeV/cm *cm^3/g *1/cm^2 = MeV/g, 
        %rho is g/cm^3, 1g/cc for water, A is beam lateral size in cm^2

    end
end

%the RF case for each phase
max_Egain=zeros(1, length(phioffsets));
for pp = 1:length(phioffsets)
    phase = phioffsets(pp);
    data = readtable(sprintf('phia_simulationsEnergyMod_phi%.2fhist.txt',phase));

    %% Extract the columns from the table
    G = data.G;
    %G= phiasimulationsEnergyModphi0.VarName6;
    G=G(~isnan(G));
    max_Egain(pp)= 938.272*(max(G)-1);
    dEdX_values = zeros(length(G), numsteps);
    dose_vals = zeros(numsteps);
    meandEdX=zeros(numsteps);
    size(dEdX_values);
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
            dose_vals(ii)= dose_vals(ii) + 1/rho*dEdX_new*1/A_beam; 


        end
    end
    %meandEdX=meandEdX/length(G);
    dose_comp= dose_comp * 1.602*10^(-13)*1000;  %MeV/g to J/kg [Gy]
    dose_vals= dose_vals * 1.602*10^(-13)*1000 ; %MeV/g to J/kg [Gy]
    %% Create a plot
    figure(gcf)
    
    %scale up from sim particles to real particles
    scaled_dose_comp= dose_comp * sim_particles_scaling;  
    scaled_dose_vals= dose_vals * sim_particles_scaling ; 

    subplot(2,1,1)
    %plot(x_values(1:numsteps), meandEdX(1:numsteps), 'b', 'LineWidth', 2,'Color',"#0072BD", 'DisplayName', 'With RF, 0.03% energy spread');
    plot(x_values(1:numsteps), scaled_dose_vals(1:numsteps), 'b', 'LineWidth', 2,'Color',"#0072BD", 'DisplayName', 'With RF, 0.03% energy spread');
    hold on;
    %plot(x_values(1:numsteps), meandEdX_comp(1:numsteps), 'b', 'LineWidth', 2,'LineStyle',':', 'Color',"#7E2F8E",'DisplayName','No RF');
    plot(x_values(1:numsteps), scaled_dose_comp(1:numsteps), 'b', 'LineWidth', 2,'LineStyle',':', 'Color',"#7E2F8E",'DisplayName','No RF, 0.03% energy spread');
    
    hold off;
    % hold on;
    %plot(x_values(1:ii), dEdX_values(1:ii), 'b', 'LineWidth', 2, 'Color', 'magenta');
    %xline(32.75, '--r', '228.5 MeV Proton', 'LineWidth', 2);
    %xline(38, '--r', '250 MeV Proton', 'LineWidth', 2);
    %legend(sprintf('%.2f phase', phase));
    %xlim([0 material_length]);
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
    hh=histogram(E, 100);
    xlabel('Energy [MeV]');
    ylabel('Simulated Particles')
    bincounts= hh.BinCounts
    ylim([0,max(bincounts)])
    yyaxis right 
    ylabel('Particles')
    bincounts_scaled=hh.BinCounts*sim_particles_scaling
    ylim([0,max(bincounts_scaled)])
    title(sprintf('phase offset is %.2f radians', phase));
    shg
end
% 
% figure
% plot(phioffsets(1:length(phioffsets)),max_Egain(1:length(phioffsets)), 'b', 'LineWidth', 2)
% xlabel('Phase')
% ylabel('Max energy reached [MeV]')
% xlim([0,6.28])
% 
% 
% max_E_diff = -(228.5 - max(max_Egain)) %MeV
% length_cell = 0.0236;  %m
% %use with input 0 energy spread
% MaxGradSeen = max_E_diff/length_cell; %peak surface E field in this case is 68 MV/m
% 
% %Cavity run with 2.5 MW into each cell
% cell_Power=2.5; %MW
% RF_power = 400; %kW
% 
% %Shunt_Imp = MaxGradSeen^2*length_cell/cell_Power %MeV^2/(m* MW)
% %https://accelconf.web.cern.ch/l06/papers/TUP044.pdf
% Avg_grad = 15*sqrt(cell_Power/0.1);  %MV/m
% MaxGradSeen/Avg_grad; %2.26
