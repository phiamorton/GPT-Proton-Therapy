phase = 0.00;  % Define phase value for grabbing the photo
energyspreadpercent= 0.03
energy0=228.5 %MeV
uniform=true
a = .005; %0.5 cm
if uniform ==true
    masterfilename = sprintf('output_EnergyMod_phi%.2f_E%.2f_Esp%.2f_uniform', phioffsetE, energy0, energyspreadpercent);
else
    masterfilename = sprintf('output_EnergyMod_phi%.2f_E%.2f_Esp%.2f', phioffsetE, energy0, energyspreadpercent);
end

for pp = 1:length(phioffsets)
    phase = phioffsets(pp)
    data = readtable(sprintf('%s.txt',masterfilename));

    %% Extract the columns from the table
    G = data.G;
    E=938.272*(G-1); %MeV
    %G=G(~isnan(G));
    x=data.x;
    y=data.y;
    z=data.z;    
    
    %plot
    figure
    % Create tiled layout
    %tiledlayout(1,3)
    
    s1=subplot(1,2,1);
    %nexttile
    scatter(x,y,25,E, 'filled')
    xlabel('x [m]')
    ylabel('y [m]')
    xlim([-a-a/10,a+a/10])
    ylim([-a-a/10,a+a/10])

    s2=subplot(1,2,2);
    %nexttile
    scatter(z,x,25,E, 'filled')
    xlabel('z [m]')
    ylabel('x [m]')
    %cb = colorbar;
    %cb.Label.String = 'Energy [MeV]';
    %hold on;

    cb = colorbar;
    cb.Label.String = 'Energy [MeV]';
    
    s1.Position(1) = s1.Position(1);
    s2.Position(1) = s2.Position(1) - 0.05;
    
    if uniform==true
        saveas(gcf,sprintf('%sBeamDist_uniform.png', masterfilename))
    else
        saveas(gcf,sprintf('%sBeamDist.png', masterfilename))
    end

end