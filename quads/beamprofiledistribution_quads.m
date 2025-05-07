phase = 0.00;  % Define phase value for grabbing the photo
phioffsetE=phase;
energyspreadpercent= 0.03;
energy0=228.5; %MeV
uniform=false;
a = .005; %0.5 cm
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
%%
toutlist = unique(avg.time);
tmargin = diff(toutlist(1:2))/2;
tgrab = toutlist(round(linspace(1,length(toutlist),35))); %[toutlist(1) toutlist(end)];
for kk = 1:length(tgrab)
    %system(['"' GPTpathname 'gdfselect.exe" -o input_dist.gdf resfull_' date '.gdf time ' num2str(tgrab(kk)-tmargin) ' ' num2str(tgrab(kk)+tmargin)]); 
    %system(['"' GPTpathname 'gdf2a.exe" -w 16 -o res.txt input_dist.gdf time x y z Bx By Bz G ID fEy GPTLICENSE=1329126328']);
    simres = readtable('res.txt');
    simres.time = simres.ID.*0+tgrab(kk);
    if kk == 1
        writetable(simres,['resfull_' date '.txt'],'WriteRowNames',true);
    else
        writetable(simres,['resfull_' date '.txt'],'WriteMode','Append','WriteVariableNames',false,'WriteRowNames',true);
    end
    
end

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

