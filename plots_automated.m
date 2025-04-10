% Read the data from the text file using readtable
phioffsets = [0.00]; %[0  0.33        0.66        0.99        1.32        1.65        1.98        2.31        2.65        2.98        3.31        3.64        3.97         4.3        4.63        4.96        5.29        5.62        5.95        6.28]; %3.4; %in rad, 0-2pi

%plotbrowser

for pp = 1:length(phioffsets)
    phase = phioffsets(pp)
    data = readtable(sprintf('phia_simulationsEnergyMod_phi%.2fhist.txt',phase));
    
    %% Extract the columns from the table
    G = data.G;
    E=938.272*(G-1); %MeV

    %translate simulated particles into real # of particles 
    Qtot0 = 4.2e-15; %from phia_test_Emod_spreadBragg.m, going to assume this is in C
    Qproton =1.6e-19; %C
    numrealprotons= Qtot0/Qproton; %this is total number of real protons
    numsimpart=2000; %from phia_test_EMod_spreadBragg.m
    sim_particles_scaling=numrealprotons/numsimpart %converting between simulated and real particles
    

    

    %% Create a plot
    figure('WindowStyle','docked', 'Name', sprintf('Energy Spectrum at Phase %.2f', phase), 'NumberTitle', 'off')
    hh=histogram(E, 100);
    xlabel('Energy [MeV]');
    ylabel('Simulated Particles')
    bincounts= hh.BinCounts;
    ylim([0,max(bincounts)]);
    yyaxis right 
    ylabel('Particles')
    bincounts_scaled=hh.BinCounts*sim_particles_scaling;
    ylim([0,max(bincounts_scaled)]);
    title(sprintf('phase offset is %.2f radians', phase));
    % % grid on;
end