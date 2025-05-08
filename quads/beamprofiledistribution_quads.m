phase = 0.00;  % Define phase value for grabbing the photo
phioffsetE=phase;
energyspreadpercent= 0.03;
energy0=228.5; %MeV
uniform=false;
a = .005; %0.5 cm
zlen0= 3*c/freq*beta0;  %in m
zposE0 = zlen0/1.8; %.104
quadpos=zposE0;
zinit=zpos-0.1;
tinit=(zinit)/beta0/c;
zfin=zpos+0.6;
tfin=(zfin)/beta0/c;
gamma0 = (energy0+938.27)/938.27; % 1.2435;
c = 2.998e8; %m/s
beta0= sqrt(1-1/(gamma0^2));
if uniform ==true
    masterfilename = sprintf('output_EnergyMod_phi%.2f_E%.2f_Esp%.2f_uniform_quads', phioffsetE, energy0, energyspreadpercent);
else
    masterfilename = sprintf('output_EnergyMod_phi%.2f_E%.2f_Esp%.2f_quads', phioffsetE, energy0, energyspreadpercent);
end

NoRF=false;
if NoRF==true
    masterfilename= sprintf('output_noRF_EnergyMod_phi%.2f_E%.2f_Esp%.2f_uniform_quads', phioffsetE, energy0, energyspreadpercent);
    ffacE=0;
end

simavg = readtable(sprintf('avgfull_%s.txt',masterfilename));
avg = table2struct(simavg,'ToScalar',true);
times=avg.time;
stdx=avg.stdx;
stdy=avg.stdy;
avgz=avg.avgz;

figure
scatter(avgz,stdx, 'Color', "#0072BD", 'DisplayName', 'average x')
hold on
scatter(avgz,stdy, 'Color', "red", 'DisplayName', 'average y')
xline(quadpos,'--','DisplayName', 'quad position', 'LineWidth',2)
legend();
xlabel('Average Z [m]');
ylabel('Transverse Profile [m]');
% 
% G = data.G;
% E=938.272*(G-1); %MeV
% G=G(~isnan(G));
% x=data.x;
% y=data.y;
% z=data.z;    
% for pp = 1:length(phioffsets)
%     phase = phioffsets(pp)
%     data = readtable(sprintf('%s.txt',masterfilename));
% 
%     %% Extract the columns from the table
%     G = data.G;
%     E=938.272*(G-1); %MeV
%     %G=G(~isnan(G));
%     x=data.x;
%     y=data.y;
%     z=data.z;    
% 
%     %plot
%     figure
%     % Create tiled layout
%     %tiledlayout(1,3)
% 
%     s1=subplot(1,2,1);
%     %nexttile
%     scatter(x,y,25,E, 'filled')
%     xlabel('x [m]')
%     ylabel('y [m]')
%     xlim([-a-a/10,a+a/10])
%     ylim([-a-a/10,a+a/10])
% 
%     s2=subplot(1,2,2);
%     %nexttile
%     scatter(z,x,25,E, 'filled')
%     xlabel('z [m]')
%     ylabel('x [m]')
%     %cb = colorbar;
%     %cb.Label.String = 'Energy [MeV]';
%     %hold on;
% 
%     cb = colorbar;
%     cb.Label.String = 'Energy [MeV]';
% 
%     s1.Position(1) = s1.Position(1);
%     s2.Position(1) = s2.Position(1) - 0.05;
% 
%     if uniform==true
%         saveas(gcf,sprintf('%sBeamDist_uniform.png', masterfilename))
%     else
%         saveas(gcf,sprintf('%sBeamDist.png', masterfilename))
%     end
% 
% end

% Bx=data.fBx;
% By=data.fBy;
% B_tot=sqrt(By.^2+Bx.^2);
% figure;
% quiver(x,y,Bx,By);
% cb = colorbar;
% cb.Label.String = 'Field Strength [T]';
% xlabel('x [m]');
% ylabel('y [m]');
% saveas(gcf,sprintf('%sQuadField.png', masterfilename))

