% runs through end of chopper slit - started 2/28/2022
% will include optional phase offset for section of linac
close all
clearvars

%%
%phioffsets = [0.00 1/4*3.14 3.14/2 3/4*3.14 3.14 3/2*3.14 6.28]; %3.4; %in rad, 0-2pi
phioffsets =  [0.00] %[0.00  0.33        0.66        0.99        1.32        1.65        1.98        2.31        2.65        2.98      3.14   3.31        3.64        3.97         4.30        4.63        4.96        5.29        5.62        5.95        6.28];  %linspace(0, 2*pi,30)
energyspreadpercent= 0.03
energy0 = 228.5 %alter energy into cavities
uniform=false
NoRF=false
ffac=false
rounded = round(phioffsets,2);
%format bank
num2str(rounded);
for pp = 1:length(phioffsets)
    phioffsetE = phioffsets(pp);
    inputfilepath = 'output_';
    fieldpathname = '""';
    GPTpathname = 'C:\bin\'; 
    ffacE = 5.5;%*10; %5.5 ;%-482; %7.5; %5.1;
    
    %Cavity run with 2.5 MW into each cell 
    %an average gradient of 30 MV/m, 400 kW of input power is required to be fed into this cell. 
    % %The peak surface E field in this case is 68 MV/m, and the peak H field is 99 kA/m. (https://doi.org/10.1063/5.0035331) 
    %E and H files should be in V/m already, need to scale to get to
    %average gradient of 30 MV/m

    if uniform ==true
        masterfilename = sprintf('EnergyMod_phi%.2f_E%.2f_Esp%.2f_uniform_quads.in', phioffsetE, energy0, energyspreadpercent);
    elseif uniform==false
        masterfilename = sprintf('EnergyMod_phi%.2f_E%.2f_Esp%.2f_quads.in', phioffsetE, energy0, energyspreadpercent);
    end

    if NoRF==true
        masterfilename= sprintf('noRF_EnergyMod_phi%.2f_E%.2f_Esp%.2f_uniform_quads.in', phioffsetE, energy0, energyspreadpercent);
        ffacE=0
    end
    if ffac==true
        masterfilename= sprintf('EnergyMod_phi%.2f_E%.2f_Esp%.2f_ffac%.2f_quads.in', phioffsetE, energy0, energyspreadpercent,ffacE);
    end

    %masterfilename;
    %sprintf('EnergyMod_phi0.00_0.03Espread_nominal.in') %for ffac=0
    %mrfilename = 'mr_test.mr';
    %date = '6_27_2024';
    
    load energyMod_phase_1_24_2022_40cells30MeV
    philistE = philist;
    %disp(philistE)
    %% Define linac parameters
    
    freq = 2.856e9;
    dcellE = 14.7*0.0254; %distance between the cells, 14.7 inches, takes input as m
    a = .005; %0.5 cm
    ncellsE = 0; %length(philistE); %changed to 2 (only 2 cell cavities)
    phasebreakE = 2;
    subplotnum = 5;
    drift = .67; %m I think
    
    %% Define beam parameters
    
    npart0 = 2000;
    
    gamma0 = (energy0+938.27)/938.27; % 1.2435;
    
    dgamma0 = (energy0*energyspreadpercent/100+938.27)/938.27-1; % .03% energy spread
    
    mevion_25nA=false;
    mevion_1nA=true;
   
    if mevion_1nA==true
        xrms0 = 3.495/1000 ;%m 4.9mm
        %based on mevion numbers this will be ~3-4mm at 1nA or 5-6mm at 25 nA
        yrms0 = 4.007/1000; %m 6mm
        %divergence of beam
        divangx0 = (3.794-3.496)/120; %.58; %change in x [mm] over 12 cm
        divangy0 = (4.299-4.007)/120; % .67;
    end  
    if mevion_25nA==true
        xrms0 = 4.906/1000 ;%m 4.9mm
        %based on mevion numbers this will be ~3-4mm at 1nA or 5-6mm at 25 nA
        yrms0 = 6.039/1000; %m 6mm
        %divergence of beam
        divangx0 = (5.198-4.906)/120; %.58; %change in x [mm] over 12 cm
        divangy0 = (6.289-6.039)/120; % .67;
    end 
    %beta0 = .5944; %v/c? 
    c = 2.998e8; %m/s
    beta0= sqrt(1-1/(gamma0^2));
    
    t_bunch= 2*10^(-6); %2 us
    %zlen0 = t_bunch*c*beta0 %in m
    %t_4rf=4/freq %4RF cycles is ~1.4e-9 seconds 
    zlen0= 3*c/freq*beta0;  %in m % will set to 3-4 RF cycles for now, actual bunch length will be 2??? us long
    %based on mevion numbers, 230 MeV beam will be 12 cm long at isocenter
    %(midway between cavities)
    emit0 = .01e-6; % 3 pi mm-mrad emittance
    Qtot0 = 4.2e-15; %in C % assumes 6 uA pulsed average current
    %current for 2us period the pulse is there
    %mevion gave average current
    zposE0 = zlen0/1.8 %.03; %what is this doing
    sc = 0;
    xoffset=0; %m
    yoffset= 0; %m
    %tdiff = .001/beta0/c;
    
    %% Initialize particle distribution entering treatment room
    if uniform == true
        buildparticles = {
        'accuracy(6);';
        ['npart = ' num2str(npart0) ';'];
        ['sc = ' num2str(sc) ';'];
        'if(npart==1){';
        ['setstartpar("beam",0,0,0,0,0,' num2str(gamma0*beta0) ',mp,-qe,' num2str(Qtot0) ');'];
        '}';
        'if(npart > 1){';
        ['setparticles("beam",' num2str(npart0) ',mp,-qe,' num2str(Qtot0) ');'];
        ['setrxydist("beam","u",' num2str(a/2) ',' num2str(a) ');'];
        ['setphidist("beam","u",0,2*pi);'];
        ['setzdist("beam","u", 0, ' num2str(zlen0) ');'];
        'setGBxdist("beam","g",0,1e-3,3,3); #primarily setting distribution shape, will be rescaled';
        'setGBydist("beam","g",0,1e-3,3,3);';
        ['setGBxemittance("beam",' num2str(emit0) ');'];
        ['setGByemittance("beam",' num2str(emit0) ');'];
        ['setGdist("beam","g",' num2str(gamma0) ',' num2str(dgamma0) ',3,3); '];
        ['setoffset("beam",' num2str(xoffset) ',' num2str(yoffset) ',0,0,0,0);'];
        ['addxdiv("beam",0,' num2str(divangx0) ');'];
        ['addydiv("beam",0,' num2str(divangy0) ');'];
        '}';
        'if(sc==1){';
        'spacecharge3dmesh();';
        '}';
        };
        
    else 
        buildparticles = {
        'accuracy(6);';
        ['npart = ' num2str(npart0) ';'];
        ['sc = ' num2str(sc) ';'];
        'if(npart==1){';
        ['setstartpar("beam",0,0,0,0,0,' num2str(gamma0*beta0) ',mp,-qe,' num2str(Qtot0) ');'];
        '}';
        'if(npart > 1){';
        ['setparticles("beam",' num2str(npart0) ',mp,-qe,' num2str(Qtot0) ');'];
        ['setxdist("beam","g",0,' num2str(xrms0) ',3,3);'];
        ['setydist("beam","g",0,' num2str(yrms0) ',3,3);'];
        ['setzdist("beam","u", 0, ' num2str(zlen0) ');'];
        'setGBxdist("beam","g",0,1e-3,3,3); #primarily setting distribution shape, will be rescaled';
        'setGBydist("beam","g",0,1e-3,3,3);';
        ['setGBxemittance("beam",' num2str(emit0) ');'];
        ['setGByemittance("beam",' num2str(emit0) ');'];
        ['setGdist("beam","g",' num2str(gamma0) ',' num2str(dgamma0) ',3,3); '];
        ['setoffset("beam",' num2str(xoffset) ',' num2str(yoffset) ',0,0,0,0);'];
        ['addxdiv("beam",0,' num2str(divangx0) ');'];
        ['addydiv("beam",0,' num2str(divangy0) ');'];
        '}';
        'if(sc==1){';
        'spacecharge3dmesh();';
        '}';
        };
    end
    %% Initialize linac iris aperture
    maxl = a*5;
    stepl = a/10;
    [X,Y,Z] = meshgrid(-maxl:stepl:maxl,-maxl:stepl:maxl,0:1);
    map.x = X(:);
    map.y = Y(:);
    map.z = Z(:);
    map.R = map.x*0;
    map.R(sqrt(map.x.^2+map.y.^2)>=a) = 1;
    
    fileID = fopen([inputfilepath 'linac_iris.txt'],'w');
    fprintf(fileID,'%s\t%s\t%s\t%s\r\n','x','y','z','R');
    for ii = 1:length(map.x)
        fprintf(fileID,'%1.5f\t%1.5f\t%1.5f\t%i\r\n',map.x(ii),map.y(ii),map.z(ii),map.R(ii));
    end
    fclose(fileID);
    
    %system(['"' 'asci2gdf.exe" -o linac_iris.gdf linac_iris.txt x 1 y 1 z ' num2str(dcellE*ncellsE) ' R 1']);
    
    %% Initialize loop over cells
    tic
    linactext = cell(2*ncellsE,1);
    zpos = zposE0;
    % for jj = 1:ncellsE
    %     if jj<phasebreakE
    %         %disp(num2str(philistE(jj)))
    %         linactext(2*(jj-1)+1:2*(jj-1)+2) = {
    %                 ['map3D_Hcomplex("wcs","z",' num2str(zpos) ', ' fieldpathname '+"Hfield_01_13_2021.gdf", "x","y","z","HxRe","HyRe","HzRe","HxIm","HyIm","HzIm", ' num2str(ffacE) ', ' num2str(philistE(jj)) ', ' num2str(2*pi*freq) ');'];
    %                 ['map3D_Ecomplex("wcs","z",' num2str(zpos) ', ' fieldpathname '+"Efield_01_13_2021.gdf", "x","y","z","ExRe","EyRe","EzRe","ExIm","EyIm","EzIm", ' num2str(ffacE) ', ' num2str(philistE(jj)) ', ' num2str(2*pi*freq) ');'];
    %                 };
    %     else
    %         linactext(2*(jj-1)+1:2*(jj-1)+2) = {
    %                 ['map3D_Hcomplex("wcs","z",' num2str(zpos) ', ' fieldpathname '+"Hfield_01_13_2021.gdf", "x","y","z","HxRe","HyRe","HzRe","HxIm","HyIm","HzIm", ' num2str(ffacE) ', ' num2str(philistE(jj)+phioffsetE) ', ' num2str(2*pi*freq) ');'];
    %                 ['map3D_Ecomplex("wcs","z",' num2str(zpos) ', ' fieldpathname '+"Efield_01_13_2021.gdf", "x","y","z","ExRe","EyRe","EzRe","ExIm","EyIm","EzIm", ' num2str(ffacE) ', ' num2str(philistE(jj)+phioffsetE) ', ' num2str(2*pi*freq) ');'];
    %                 };
    %     end
    %     zpos = zpos+dcellE; 
    % end


    %quadrupole strength in the unit of T/m~~~ dimension is IMPORTANT
    length_quad = 0.2062;
    quadpos=[0.2,0.5];
    gq1 = -2; %~36kG/m + focuses in x and - focuses in y
    %gq2 = -0.0001;
    %gq3 = 0.0001;

    
    inputfiletext=[{ ['quadrupole( "wcs","z",' num2str(zpos) ',' num2str(length_quad) ',' num2str(gq1) ');'] }];
    %quadrupole( "wcs","z", 2.84, length_quad, gq2) ;
    %quadrupole( "wcs","z", 3.68, length_quad, gq3) ;
    %{ ['quadrupole("wcs","z",' num2str(quadpos(3)) ',' num2str(length_quad) ',' num2str(gq1) ');'] };
    inputfiletext = [buildparticles; 
        { ['quadrupole("wcs","z",' num2str(quadpos(1)) ',' num2str(length_quad) ',' num2str(gq1) ');'] }; 
        { ['quadrupole("wcs","z",' num2str(quadpos(2)) ',' num2str(length_quad) ',' num2str(-gq1) ');'] };{
        ['map3D_remove("wcs","z",' num2str(zposE0-dcellE/2) ', ' fieldpathname '+"linac_iris.gdf", "x","y","z","R") ;'];
        }; linactext; {
        ['tout(' num2str((zpos-0.1)/beta0/c) ',' num2str((zpos+2)/beta0/c) ',' num2str((0.01)/beta0/c)  ');']; 
        };];
    %'tout(' num2str((2*drift)/beta0/c) ');'
    
    % Write input file
    fileID = fopen([inputfilepath masterfilename],'wt');
    for ii = 1:length(inputfiletext)
    fprintf(fileID,'%s \n',inputfiletext{ii});
    end
    fclose(fileID); 
end

%run the GPT script
system('bash "sim_auto.bat"');

%do the plotting
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

gq1abs=abs(gq1);
figure
scatter(avgz,stdx*1000, 'Color', "#0072BD", 'DisplayName', 'average x')
hold on
scatter(avgz,stdy*1000, 'Color', "red", 'DisplayName', 'average y')
xline(quadpos(1),'-','DisplayName', 'quad position 1', 'LineWidth',2)
xline(quadpos(2),'-','DisplayName', 'quad position 2', 'LineWidth',2)
%xline(quadpos(3),'-','DisplayName', 'quad position 3', 'LineWidth',2)
legend();
xlabel('Average Z [m]');
ylabel('Transverse Profile [mm]');
hold off
title(sprintf('Transverse profile, strength= %.2f T/m, positions are %.2f & %.2f m', gq1abs, quadpos(1),quadpos(2)), 'FontSize', 8);
saveas(gcf,sprintf('%sFODO.png', masterfilename))
