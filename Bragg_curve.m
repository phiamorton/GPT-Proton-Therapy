I = 75 *10^(-6); % MeV or 80.8+-0.3
rho = 1/18*(100)^2; % mol/(m^2cm) for water
TZ = 10; % target Z
TA = 18; % target A
NA = 6.023e23; % Avogadro's number
re = 2.817e-13; % classical electron radius in centimeters
mec2 = 0.511; % MeV (rest mass energy of electron)
Z = 1; % proton
q = 1; % charge of proton
e_0 = 931.5; % MeV
E_0 = 228.5; % energy of particle in MeV
A = 1; % A for proton=1
n= 3.34*10^23 %electron density of water in 1/cc
const = 4*pi*n*re^2*mec2*q^2  %MeV/cm %4 * pi * NA * re^2 * mec2  * Z * q^2 / TA %MeV/cm

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
material_length = 40; % cm i hope
numsteps = 1000000; % Increase the number of steps for better resolution
dx = material_length / numsteps;
x_values = linspace(0, material_length, numsteps); %cm
dEdX_values = zeros(1, numsteps);
dEdX_values(1) = calc_stoppingpower(E_0, A, I, e_0, rho, const, Z, mec2, TA);
E_new = E_0;

% Loop to calculate energy loss and stopping power
for ii = 2:numsteps
    % Calculate stopping power at the current energy
    dEdX = calc_stoppingpower(E_new, A, I, e_0, rho, const, Z, mec2, TA);
    
    % Energy loss at each step
    E_loss = dEdX * dx;
    
    % Update proton's energy (E decreases)
    E_new = E_new - E_loss;
    
    % Store the new stopping power value at the current step
    dEdX_values(ii) = calc_stoppingpower(E_new, A, I, e_0, rho, const, Z, mec2, TA);
    
    % If energy becomes too low, break out of the loop
    if E_new <= 0
        break;
    end
end

% Plot the Bragg curve as a function of position
figure;
plot(x_values(1:ii), dEdX_values(1:ii), 'b', 'LineWidth', 2);
xline(38, '--r', 'Expected Range (~38 cm)', 'LineWidth', 2);
legend(sprintf('%.2f MeV proton', E_0));
xlabel('Position [cm]', 'FontSize', 12);
ylabel('Stopping Power (dE/dx) [MeV/cm]', 'FontSize', 12);
title('Bragg Curve for Proton in Water', 'FontSize', 14);
xlim([30,35])
grid on;
%A 250 MeV proton has a penetration depth of approximately 38 centimeters
%in water 