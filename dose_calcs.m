%https://gray.mgh.harvard.edu/media/com_dpattachments/attachments/com_content.article/Techniques-of-Proton-Radiotherapy-04-Basics.pdf
% Read the data from the text file using readtable
phioffsets = [0.00] ;%[0  0.33        0.66        0.99        1.32        1.65        1.98        2.31        2.65        2.98      3.14   3.31        3.64        3.97         4.30        4.63        4.96        5.29        5.62        5.95        6.28]; %3.4; %in rad, 0-2pi
 %[0.33]; %3.4; %in rad, 0-2pi
%[0  0.33        0.66        0.99        1.32        1.65        1.98        2.31        2.65        2.98        3.31        3.64        3.97         4.3        4.63        4.96        5.29        5.62        5.95        6.28]; %3.4; %in rad, 0-2pi
energyspreadpercent= 0.03
energy0=228.5; %MeV
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
%rad_beam = 0.5/2; %cm - 5mm beam diameter
%A_beam= pi*rad_beam^2; %cm^2

mevion_25nA=false;
mevion_1nA=true;

if mevion_1nA==true
    divangx0 = (3.794-3.496)/120; %.58; %change in x [mm] over 12 cm
    divangy0 = (4.299-4.007)/120; % .67;
    xiso=3.604/1000 ;%m 4.9mm
    yiso = 4.129/1000; %m 6mm
    divangx0*isocenter;
    divangy0*isocenter;
    xrms0 = xiso-divangx0*isocenter;
    yrms0= yiso-divangy0*isocenter;
    %gives y=4.15 at iso and x=3.61
    
end  
if mevion_25nA==true
    divangx0 = (5.198-4.906)/120; %.58; %change in x [mm] over 12 cm
    divangy0 = (6.289-6.039)/120; % .67;
    xiso=5.05/1000; %m 4.9mm
    yiso = 6.23/1000 ;%m 6mm
    divangx0*isocenter;
    divangy0*isocenter;
    xrms0 = xiso-divangx0*isocenter;
    yrms0= yiso-divangy0*isocenter;
    %remember, div in GPT MUST be rad/m ie div/sizerms
    %this actually gives isocenter yrms0=6.28 and xrms0=5.09
    %based on mevion numbers this will be ~3-4mm at 1nA or 5-6mm at 25 nA
    
end 

yrms0=yrms0*100; %cm
xrms0=xrms0*100; %cm

A_beam=pi*xrms0*yrms0; %cm^2
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
        dEdX = const / (beta2) * (log(W_max_new /(I)) - beta2);
    end
end

% Material and simulation parameters
material_length = 80; % cm i hope
numsteps = 10000; % Increase the number of steps for better resolution
dx = material_length / numsteps;
x_values = linspace(0, material_length, numsteps); %cm
%plotbrowser

comparison_noE = readtable(sprintf('comparison_phi0.00_E228_Esp0.03.txt'));
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
        E_loss = dEdX_comp * dx;

        % Update proton's energy (E decreases)
        E_new = E_new - E_loss;

        % Store the new stopping power value at the current step
        dEdX_new_comp=calc_stoppingpower(E_new, A, I, e_0, rho, const, Z, mec2, TA);
        %mass_stop_power_nominal = 1/rho * meandEdX_comp %MeVcm^2/g
        %dEdX_values(proton, ii) = dEdX_new;

        meandEdX_comp(ii) = meandEdX_comp(ii)+ dEdX_new_comp;
        dose_comp(ii) = dose_comp(ii) +1/rho*dEdX_new_comp*1/A_beam; %MeV/cm *cm^3/g *1/cm^2 = MeV/g

    end
end


max_Egain=zeros(1, length(phioffsets));
for pp = 1:length(phioffsets)
    phase = phioffsets(pp)
    data = readtable(sprintf('%s.txt',masterfilename));

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
            E_loss = dEdX * dx;

            % Update proton's energy (E decreases)
            E_new = E_new - E_loss;

            % Store the new stopping power value at the current step
            dEdX_new=calc_stoppingpower(E_new, A, I, e_0, rho, const, Z, mec2, TA);
            %dEdX_values(proton, ii) = dEdX_new;

            meandEdX(ii) = meandEdX(ii)+ dEdX_new;
            dose_vals(ii)= dose_vals(ii) + 1/rho*dEdX_new*1/A_beam; 


        end
    end
    figure;
    %meandEdX=meandEdX/length(G);
    dose_comp= dose_comp * 1.602*10^(-13)*1000;  %MeV/g to J/kg [Gy]
    dose_vals= dose_vals * 1.602*10^(-13)*1000 ; %MeV/g to J/kg [Gy]
    y=linspace(-yrms0*15,yrms0*15,numsteps);
    [Z,Y]=meshgrid(dose_vals(1:numsteps),y);
    size(dose_vals(1:numsteps));
    D=Z .* 1/(yrms0*sqrt(2*pi)).*exp(-Y.^2/(2*yrms0^2));
    %figure(gcf)
    imagesc(x_values,y,D);
    ylabel('Transverse size [cm]');
    xlabel('Depth in material [cm]');
    xmin=28;
    realyrange=8.7297;
    xlim([xmin,xmin+11.7]); %for camera res
    
    ylim([-realyrange/2,realyrange/2]);
    %set(gca,'XTick',[], 'YTick', [],'Visible', 'off', 'Color','none');
    fontsize(35,"points");
    a=colorbar;
    a.Label.String = 'Dose [a.u.]';
    %title(sprintf('Bragg Curve for Proton in Water, phase offset = %.2f rad', phase), 'FontSize', 14);
    saveas(gcf,sprintf('%sBraggIm.png', masterfilename))
    
end

