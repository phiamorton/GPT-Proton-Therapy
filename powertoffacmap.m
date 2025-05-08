phioffsetE = [0.00]; %[0  0.33        0.66        0.99        1.32        1.65        1.98        2.31        2.65        2.98      3.14   3.31        3.64        3.97         4.30        4.63        4.96        5.29        5.62        5.95        6.28]; %3.4; %in rad, 0-2pi
 %[0.33]; %3.4; %in rad, 0-2pi
%[0  0.33        0.66        0.99        1.32        1.65        1.98        2.31        2.65        2.98        3.31        3.64        3.97         4.3        4.63        4.96        5.29        5.62        5.95        6.28]; %3.4; %in rad, 0-2pi
energyspreadpercent= 0.03; % in %
energy0=228.5; %MeV
ffacs= [4.5,5.5,6.5]
l=length(ffacs(:));
powers=zeros([l,1]);
for pp = 1:length(ffacs)
    ffac = ffacs(pp);
    disp(ffac)
    masterfilename = sprintf('output_EnergyMod_phi%.2f_E%.2f_Esp%.2f_ffac%.2f', phioffsetE, energy0, energyspreadpercent,ffac);
    %masterfilename = sprintf('output_EnergyMod_phi%.2f_E%.2f_Esp%.2f', phioffsetE, energy0, energyspreadpercent);
    datafile=sprintf('%s.txt',masterfilename);
    masterfilename
    data = readtable(datafile);
    %% Extract the columns from the table
    G = data.G;
    %G= phiasimulationsEnergyModphi0.VarName6;
    G=G(~isnan(G));
    max_Egain= 938.272*(max(G)-1);
    max_E_diff = -(energy0 - max_Egain) %MeV
    length_cell = 0.0236;  %m
    % %use with input 0 energy spread
    MaxGradSeen = max_E_diff/(2*length_cell); %MV/m %2 cells %peak surface E field in this case is 68 MV/m
    % Compute the average energy gain
    %AvgEGain = (mean(E_0 - energy0))*10^6;
    %AvgGrad=AvgEGain/(2*length_cell);
    shuntimpedance= 54.8e6; %Ohm/m
    fprintf('Max gradient %.2f MV/m \n', MaxGradSeen)
    %p_disp=max_E_diff^2/(2*shuntimpedance*2*length_cell)
    p_disp=(MaxGradSeen*1000000)^2/(2*shuntimpedance); %W %%should there be a factor of 2??
    fprintf('Power calculated from max gradient %.2f MW \n', p_disp/1e6)
    fprintf('Power calculated from max gradient %.2f MW per cavity \n', p_disp/1e6/2)
    Avg_grad_calc = 15*sqrt(p_disp/100000);  %MV/m, assumes P is in W
    fprintf('Average gradient calculated from power %.2f MV/m \n', Avg_grad_calc)
    powers(pp)=p_disp;
    %
end

figure;
line(ffacs,powers/1e6)
hold on
scatter(ffacs,powers/1e6)
xlabel('ffac value')
ylabel('power into both cavities [MW]')
saveas(gcf,sprintf('%sPowerToFFAC.png', masterfilename))

%power_avg= AvgGrad^2/(2*shuntimpedance); %MV^2/m*MOhm/m = V^2/Ohm= W 
    %power loss of klystron/power into cavities ?
    %fprintf('Power calculated from average gradient %.2f MW \n', power_avg/1e6)
    %fprintf('Power calculated from average gradient %.2f MW per cavity \n', power_avg/1e6/2)
    %disp(power_avg)
    %power_max = max_E_diff^2/(2*length_cell*shuntimpedance)  %MV^2/m*MOhm/m = V^2/Ohm= W %power loss of klystron/power into cavities
    %https://cds.cern.ch/record/1005047/files/p145.pdf 
    % %Cavity run with 2.5 MW into each cell
    % %Shunt_Imp = MaxGradSeen^2*length_cell/cell_Power %MeV^2/(m* MW)
    % %https://accelconf.web.cern.ch/l06/papers/TUP044.pdf
    