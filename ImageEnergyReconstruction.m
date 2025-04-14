phase = 0.00;  % Define phase value for grabbing the photo
energyspreadpercent= 0.03
filename = sprintf('BraggIm%.2f_E%.2f.png', phase, energyspreadpercent);  % Input file name
inpict = imread(filename);  % Read image

% Convert to grayscale to identify the white border
grayImg = rgb2gray(inpict);

% Thresholding: Identify non-white pixels (assuming white is near max intensity)
threshold = 250; % Adjust if needed (255 is pure white)
mask = grayImg < threshold; % Logical mask of non-white pixels

% Find non-white pixel bounds
[row, col] = find(mask); % Get indices of non-white pixels
top = min(row);
bottom = max(row);
left = min(col);
right = max(col);

% Crop the image using the detected bounding box
croppedImg = inpict(top:bottom, left:right, :);

% Resize the image to a fixed size (ie 300x300 pixels)
newWidth = 1292;
newHeight = 964;
%real camera calibration ruler photo 1292x964 pixels, 18.5cm-6.8cm
realwidth=18.5-6.8; %in width
realpixelsize=realwidth/newWidth;
realyrange=realpixelsize*newHeight;
resizedImg = imresize(croppedImg, [newHeight, newWidth]);

figure;
% Save the resized image
pixelated_filename = sprintf('BraggIm%.2f_E%.2fPixel.png', phase, energyspreadpercent);
imshow(resizedImg);

saveas(gcf, pixelated_filename);
inpixsize= size(imread(filename));
pix_size=size(imread(pixelated_filename));
xmin=28;
xmax=xmin+11.7;
ymin=-realyrange/2;
ymax=realyrange/2;
x_range= xmax-xmin; %cm
y_range= ymax-ymin; %cm

% Convert to grayscale for intensity analysis
grayResizedImg = rgb2gray(resizedImg);
graysize=size(grayResizedImg)
% Integrate along the y-dimension (sum over rows)
intensityProfile = sum(grayResizedImg, 1);  % Sum each column along y
intensityProfile= intensityProfile-min(intensityProfile); %rescale(intensityProfile,0,max(intensityProfile));
max_intensity=max(intensityProfile);
imagesize=size(intensityProfile);
pixelsizex= x_range/graysize(1) %cm
pixelsizey= y_range/graysize(2); %cm %change based on image scale
% Plot the integrated intensity vs. x-position
x = 1:newWidth;
x_cm= linspace(xmin, xmax ,newWidth);  %xlim([15,38]) steal from image boundaries
figure('Visible','off');
plot(x_cm, intensityProfile, 'k', 'LineWidth', 2);
title(sprintf('Reconstructed Dose vs Depth at Relative Phase %.2f',phase))
xlabel('Depth [cm]');
ylabel('Dose [a.u.]');
%title(sprintf('Phase = %.2f', phase));
grid on;

% Save the intensity plot
intensity_filename = sprintf('BraggIm%.2f_E%.2fIntensity.png', phase, energyspreadpercent);
saveas(gcf, intensity_filename);

% %map depth to original particle energy--> do this by taking an array of
% %initial energies and using dE/dX calculate stopping position, use this as
% %a map to take position(i)-->energy(i)
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

material_length = 45; % cm 
numsteps = 10000; % Increase the number of steps for better resolution
dx = material_length / numsteps;
x_values = linspace(0, material_length, numsteps); %cm
G=linspace(100,270, 10000);
stop_pos = zeros(1,length(G));

for proton = 1:length(G)
    E_0 =G(proton); %MeV, starting energy
    % Loop to calculate energy loss and stopping power
    E_new=E_0;
    if isnan(E_new)==1
        break;
    end
    for ii = 2:numsteps
        if E_new <= 0 
            stop_pos(proton)=x_values(ii);
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

    end
end

function dose_resized = calculate_energy_loss(energy, numsteps, x_values, A, I, e_0, rho, const, Z, mec2, TA, dx,xmin,xmax, A_beam,x_cm)
    % Function to calculate energy loss of a proton in a medium
    % 
    % Inputs:
    %   G        - Function handle to get initial energy
    %   proton   - Proton identifier/index
    %   numsteps - Number of steps in the simulation
    %   x_values - Array of position values
    %   A, I, e_0, rho, const, Z, mec2, TA - Physical parameters for stopping power
    %   dx       - Step size
    %
    % Output:
    %   stop_pos - Position where the proton stops

    % Get initial energy
    E_0 = energy; % MeV, starting energy
    %initial array for dEdX
    dose=zeros(numsteps);
    E_new = E_0; % Initialize energy
    stop_pos = 0; % Default to 0 if proton does not stop within range

    % Check if the initial energy is NaN
    if isnan(E_new)
        return; % Exit the function
    end

    % Loop to calculate energy loss and stopping power
    for ii = 2:numsteps
        % If energy is depleted, record stopping position and exit loop
        if E_new <= 0
            stop_pos = x_values(ii);
            break;
        end
        
        % Calculate stopping power at current energy
        dEdX = calc_stoppingpower(E_new, A, I, e_0, rho, const, Z, mec2, TA);

        % Energy loss at this step
        E_loss = dEdX * dx;
        
        % Update energy
        E_new = E_new - E_loss;
        dEdX_new=calc_stoppingpower(E_new, A, I, e_0, rho, const, Z, mec2, TA);
        dose(ii)= dose(ii) + 1/rho*dEdX_new*1/A_beam;
    end
    dose = dose* 1.602*10^(-13)*1000;  %MeV/g to J/kg [Gy]
    %dose here is across the x range 0-80
    %need to make a cut to only the same range as the photo range
    % Cut x_values to keep only values between xmin and xmax
    % Reshape to match x_cm using interpolation
    %x_cm = linspace(xmin, xmax, newWidth); % New evenly spaced points
    mask = (x_values >= xmin & x_values <= xmax); % Logical mask
    x_cut = x_values(mask);
    dose_cut = dose(mask);
    dose_resized = interp1(x_cut, dose_cut, x_cm, 'linear');
end


%take an x value in the camera and find its energy stopping range value and
%multiply by the corresponding intensity value

energy_spectrum=zeros(length(x_cm));
energy_intensity =zeros(length(x_cm));

figure('Visible','off')
for x_pos = length(x_cm):-1:1 %start at high x values, furthest out
    intensity=intensityProfile(x_pos);
    if intensity*pixelsizex> 0.0001* max_intensity
        intensity;
        
        %find energy from x position map
        x_val=x_cm(x_pos);
        
        [~, index] = min(abs(stop_pos - x_val));
        index;
        
        energy=G(index);
        
        %Next remove the intensity spectrum from that particle,
        %Otherwise it will cound the low energy buildup as coming from lots of
        %very low energy particles
        dose_resized = calculate_energy_loss(energy, numsteps, x_values, A, I, e_0, rho, const, Z, mec2, TA, dx,xmin,xmax, A_beam, x_cm);
        
        dose_resized=rescale(dose_resized,0, 1)*intensity;
        
        
    
        intensityProfile=intensityProfile-dose_resized; %need to rescale differently

        plot(x_cm,intensityProfile,'DisplayName', num2str(x_val));
            
        legend
        title('intensity after subtraction');
        hold on
        
        %Create spectrum for particle of energy G(proton)- lets take the time
        %to make this a function
        energy_spectrum(x_pos)= energy ;
        energy_intensity(x_pos)=intensity;
        %remove from the energy spectrum
    end
    % if intensity<0
    %     break;
    % end
end

hold off
energy_spectrum(1:length(x_cm));

figure
subplot(2,1,2); 
data = readtable(sprintf('phia_simulationsEnergyMod_phi%.2f_E%.2fhist.txt',phase, energyspreadpercent));

%% Extract the columns from the table
G = data.G;
E=938.272*(G-1); %MeV
%% Create a plot
%figure('WindowStyle','docked', 'Name', sprintf('Energy Spectrum at Phase %.2f', phase), 'NumberTitle', 'off')
histogram(E, 100);
xlabel('Energy [MeV]');
ylabel('Simulated Particles');
title(sprintf('Phase is %.2f radians, energy spread %.2f %', phase, energyspreadpercent));
% % grid on; 
rangemin=min(E)-0.3;
rangemax=max(E)+0.3;
xlim([rangemin, rangemax]);
%energy_spectrum_filtered=energy_spectrum(energy_spectrum>0);
%histogram(energy_spectrum(1:length(x_cm)))
subplot(2,1,1);
for x_pos = 1:length(x_cm)
    if energy_spectrum(x_pos)>0
        bar(energy_spectrum(x_pos),energy_intensity(x_pos),'BarWidth',pixelsizex, 'FaceColor','m');
        hold on;
    end
end

title(sprintf('Reconstructed Energy Spectrum, Phase is %.2f and pixel size is %.4f cm for %.0f x %.0f Pixels', phase, pixelsizex,graysize(1),graysize(2)));
%xlim([227,230])
xlabel('Energy [MeV]');
ylabel('Intensity [a.u.]');
xlim([rangemin, rangemax]);
hold off;
%saveas(f,'ReconstructedVSSimulatedEnergy','fig')