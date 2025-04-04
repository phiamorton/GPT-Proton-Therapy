% Read the data from the text file using readtable
phioffsets = [2.65]; %[0  0.33        0.66        0.99        1.32        1.65        1.98        2.31        2.65        2.98        3.31        3.64        3.97         4.3        4.63        4.96        5.29        5.62        5.95        6.28]; %3.4; %in rad, 0-2pi

%plotbrowser

for pp = 1:length(phioffsets)
    phase = phioffsets(pp)
    data = readtable(sprintf('phia_simulationsEnergyMod_phi%.2fhist.txt',phase));
    
    %% Extract the columns from the table
    G = data.G;
    E=938.272*(G-1); %MeV
    %% Create a plot
    figure('WindowStyle','docked', 'Name', sprintf('Energy Spectrum at Phase %.2f', phase), 'NumberTitle', 'off')
    histogram(E, 100)
    xlabel('Energy [MeV]');
    ylabel('Simulated Particles')
    title(sprintf('phase offset is %.2f radians', phase));
    % % grid on;
end