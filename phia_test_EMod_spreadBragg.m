% runs through end of chopper slit - started 2/28/2022
% will include optional phase offset for section of linac
close all
clearvars

%%
%phioffsets = [0.00 1/4*3.14 3.14/2 3/4*3.14 3.14 3/2*3.14 6.28]; %3.4; %in rad, 0-2pi
phioffsets =  [0.00] %[0.00  0.33        0.66        0.99        1.32        1.65        1.98        2.31        2.65        2.98      3.14   3.31        3.64        3.97         4.30        4.63        4.96        5.29        5.62        5.95        6.28];  %linspace(0, 2*pi,30)
energyspreadpercent= 0.03
energy0 = 228.5 %alter energy into cavities
uniform=true

rounded = round(phioffsets,2);
%format bank
num2str(rounded);
for pp = 1:length(phioffsets)
    phioffsetE = phioffsets(pp);
    inputfilepath = 'output_';
    fieldpathname = '""';
    GPTpathname = 'C:\bin\'; 
    if uniform ==true
        masterfilename = sprintf('EnergyMod_phi%.2f_E%.2f_Esp%.2f_uniform.in', phioffsetE, energy0, energyspreadpercent);
    else
        masterfilename = sprintf('EnergyMod_phi%.2f_E%.2f_Esp%.2f.in', phioffsetE, energy0, energyspreadpercent);
    end
    %sprintf('EnergyMod_phi0.00_0.03Espread_nominal.in') %for ffac=0
    mrfilename = 'mr_test.mr';
    date = '6_27_2024';
    
    load energyMod_phase_1_24_2022_40cells30MeV
    philistE = philist;
    %disp(philistE)
    %% Define linac parameters
    
    freq = 2.856e9;
    dcellE = 14.7*0.0254; %distance between the cells, 14.7 inches, takes input as m
    a = .005; %0.5 cm
    ffacE = 5.5%5.5 ;%-482; %7.5; %5.1;
    %Cavity run with 2.5 MW into each cell 
    %an average gradient of 30 MV/m, 400 kW of input power is required to be fed into this cell. 
    % %The peak surface E field in this case is 68 MV/m, and the peak H field is 99 kA/m. (https://doi.org/10.1063/5.0035331) 
    %E and H files should be in V/m already, need to scale to get to
    %average gradient of 30 MV/m
    ncellsE = 2; %length(philistE); %changed to 2 (only 2 cell cavities)
    phasebreakE = 2;
    subplotnum = 5;
    drift = .67; %m I think
    
    %% Define beam parameters
    
    npart0 = 2000;
    
    gamma0 = (energy0+938.27)/938.27; % 1.2435;
    
    dgamma0 = (energy0*energyspreadpercent/100+938.27)/938.27-1; % .03% energy spread
    
    rad_beam = 0.5/2; %cm - 5mm beam diameter
    %beta0 = .5944; %v/c? 
    c = 2.998e8; %m/s
    beta0= sqrt(1-1/(gamma0^2));
    xrms0 = rad_beam/100; %m I think;
    yrms0 = rad_beam/100; %.0035;
    t_bunch= 2*10^(-6); %2 us
    %zlen0 = t_bunch*c*beta0 %in m
    %t_4rf=4/freq %4RF cycles is ~1.4e-9 seconds 
    zlen0= 1*c/freq*beta0  %in m % will set to 3-4 RF cycles for now, actual bunch length will be 2??? us long
    divangx0 = 0; %.58;
    divangy0 = 0; % .67;
    emit0 = .01e-6; % 3 pi mm-mrad emittance
    Qtot0 = 4.2e-15; %in C % assumes 6 uA pulsed average current
    zposE0 = zlen0/1.8; %.03; %what is this doing
    sc = 0;
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
    for jj = 1:ncellsE
        if jj<phasebreakE
            %disp(num2str(philistE(jj)))
            linactext(2*(jj-1)+1:2*(jj-1)+2) = {
                    ['map3D_Hcomplex("wcs","z",' num2str(zpos) ', ' fieldpathname '+"Hfield_01_13_2021.gdf", "x","y","z","HxRe","HyRe","HzRe","HxIm","HyIm","HzIm", ' num2str(ffacE) ', ' num2str(philistE(jj)) ', ' num2str(2*pi*freq) ');'];
                    ['map3D_Ecomplex("wcs","z",' num2str(zpos) ', ' fieldpathname '+"Efield_01_13_2021.gdf", "x","y","z","ExRe","EyRe","EzRe","ExIm","EyIm","EzIm", ' num2str(ffacE) ', ' num2str(philistE(jj)) ', ' num2str(2*pi*freq) ');'];
                    };
        else
            linactext(2*(jj-1)+1:2*(jj-1)+2) = {
                    ['map3D_Hcomplex("wcs","z",' num2str(zpos) ', ' fieldpathname '+"Hfield_01_13_2021.gdf", "x","y","z","HxRe","HyRe","HzRe","HxIm","HyIm","HzIm", ' num2str(ffacE) ', ' num2str(philistE(jj)+phioffsetE) ', ' num2str(2*pi*freq) ');'];
                    ['map3D_Ecomplex("wcs","z",' num2str(zpos) ', ' fieldpathname '+"Efield_01_13_2021.gdf", "x","y","z","ExRe","EyRe","EzRe","ExIm","EyIm","EzIm", ' num2str(ffacE) ', ' num2str(philistE(jj)+phioffsetE) ', ' num2str(2*pi*freq) ');'];
                    };
        end
        zpos = zpos+dcellE; 
    end
    inputfiletext = [buildparticles; {
        ['map3D_remove("wcs","z",' num2str(zposE0-dcellE/2) ', ' fieldpathname '+"linac_iris.gdf", "x","y","z","R") ;'];
        }; linactext; {
        ['tout(' num2str((2*drift)/beta0/c) ');']; 
        };];
    
    
    % Write input file
    fileID = fopen([inputfilepath masterfilename],'wt');
    for ii = 1:length(inputfiletext)
    fprintf(fileID,'%s \n',inputfiletext{ii});
    end
    fclose(fileID); 
end
