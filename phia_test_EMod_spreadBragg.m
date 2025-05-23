% runs through end of chopper slit - started 2/28/2022
% will include optional phase offset for section of linac
close all
clearvars

%%
%phioffsets = [0.00 1/4*3.14 3.14/2 3/4*3.14 3.14 3/2*3.14 6.28]; %3.4; %in rad, 0-2pi
phioffsets =  [0.00] ;%[0.00  0.33        0.66        0.99        1.32        1.65        1.98        2.31        2.65        2.98      3.14   3.31        3.64        3.97         4.30        4.63        4.96        5.29        5.62        5.95        6.28];  %linspace(0, 2*pi,30)
energyspreadpercent= 0.03;
energy0 = 228.5 ;%alter energy into cavities
uniform=false;
NoRF=false;
ffac=false;
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
        masterfilename = sprintf('EnergyMod_phi%.2f_E%.2f_Esp%.2f_uniform', phioffsetE, energy0, energyspreadpercent);
    elseif uniform==false
        masterfilename = sprintf('EnergyMod_phi%.2f_E%.2f_Esp%.2f', phioffsetE, energy0, energyspreadpercent);
    end

    if NoRF==true
        masterfilename= sprintf('noRF_EnergyMod_phi%.2f_E%.2f_Esp%.2f', phioffsetE, energy0, energyspreadpercent);
        ffacE=0;
    end
    if ffac==true
        masterfilename= sprintf('EnergyMod_phi%.2f_E%.2f_Esp%.2f_ffac%.2f', phioffsetE, energy0, energyspreadpercent,ffacE);
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
    a = 0.05; %5 cm, this is the iris radius in m
    ncellsE = 2; %length(philistE); %changed to 2 (only 2 cell cavities)
    phasebreakE = 2;
    subplotnum = 5;
    drift = .67; %m I think
    
    %% Define beam parameters
    
    npart0 = 2000; %number of particles in simulation
    
    gamma0 = (energy0+938.27)/938.27; % 1.2435;
    
    dgamma0 = (energy0*energyspreadpercent/100+938.27)/938.27-1; % .03% energy spread
    
    %beta0 = .5944; %v/c? 
    c = 2.998e8; %m/s
    beta0= sqrt(1-1/(gamma0^2));
    
    t_bunch= 2*10^(-6); %2 us
    %zlen0 = t_bunch*c*beta0 %in m
    %t_4rf=4/freq %4RF cycles is ~1.4e-9 seconds 
    zlen0= 3*c/freq*beta0 ; %in m % will set to 3-4 RF cycles for now, actual bunch length will be 2??? us long
    %based on mevion numbers, 230 MeV beam will be 12 cm long at isocenter
    %(midway between cavities)
    emit0 = .01e-6; % 3 pi mm-mrad emittance
    Qtot0 = 4.2e-15; %in C % assumes 6 uA pulsed average current
    %current for 2us period the pulse is there
    %mevion gave average current
    zposE0 = zlen0/1.8; %.03; %position of first cavity
    sc = 0;
    xoffset=0; %m
    yoffset= 0; %m
    %tdiff = .001/beta0/c;
    isocenter = zposE0+dcellE/2 ;%m

    mevion_25nA=false;
    mevion_1nA=true;
   
    if mevion_1nA==true
        divangx0 = (3.794-3.496)/120; %.58; %change in x [mm] over 12 cm
        divangy0 = (4.299-4.007)/120; % .67;
        xiso=3.604/1000; %m 4.9mm
        yiso = 4.129/1000; %m 6mm
        divangx0*isocenter;
        divangy0*isocenter;
        xrms0 = xiso-divangx0*isocenter;
        yrms0= yiso-divangy0*isocenter;
        %gives y=4.15 at iso and x=3.61
        
    end  
    if mevion_25nA==true
        divangx0 = (5.198-4.906)/120; %.58; %change in x [mm] over 12 cm
        divangy0 = (6.289-6.039)/120; % .67;
        xiso=5.05/1000 ;%m 4.9mm
        yiso = 6.23/1000; %m 6mm
        divangx0*isocenter;
        divangy0*isocenter;
        xrms0 = xiso-divangx0*isocenter;
        yrms0= yiso-divangy0*isocenter;
        %remember, div in GPT MUST be rad/m ie div/sizerms
        %this actually gives isocenter yrms0=6.28 and xrms0=5.09
       
        %based on mevion numbers this will be ~3-4mm at 1nA or 5-6mm at 25 nA
        
        %divergence of beam
        
    end 
    
    
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
        ['addxdiv("beam",0,' num2str(divangx0/xiso) ');'];
        ['addydiv("beam",0,' num2str(divangy0/yiso) ');'];
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
    outputlengthfrom2ndcav=6; %in
    inputfiletext = [buildparticles; %{
        %['map3D_remove("wcs","z",' num2str(zposE0-dcellE/2) ', ' fieldpathname '+"linac_iris.gdf", "x","y","z","R") ;'];
        %}; %this is the iris
        linactext; {
        ['tout(' num2str((outputlengthfrom2ndcav*.0254+zposE0+dcellE)/beta0/c) ');']; 
        };];
    %tout is .152m (6in) + center position of 2nd cavity
    
    % Write input file
    fileID = fopen([inputfilepath sprintf('%s.in', masterfilename)],'wt');
    for ii = 1:length(inputfiletext)
    fprintf(fileID,'%s \n',inputfiletext{ii});
    end
    fclose(fileID); 
end

%run gpt and create .txt file output
%system('bash "sim_auto.bat"'); %make sure filepath here matches filepath used
system(sprintf('gpt -o output_%s.gdf output_%s.in', masterfilename, masterfilename))
%system(sprintf('gdf2his -o output_%shist.gdf output_%s.gdf G 0.00001',masterfilename, masterfilename))
system(sprintf('gdf2a -w 16 -o output_%s.txt output_%s.gdf', masterfilename, masterfilename))
system(sprintf('ex -s -c "1d3|x" output_%s.txt', masterfilename))

datafile=sprintf('output_%s.txt',masterfilename);
data = readtable(datafile);
%% Extract the columns from the table
Gout = data.G;
xout = data.x;
yout = data.y;
zout = data.z;
Bxout = data.Bx; %velocity beta (speed as % speed of light c)
Byout = data.By;
Bzout = data.Bz;
IDout = data.ID;

% Assuming all variables are column vectors or can be reshaped into column vectors
T = table(IDout(:), xout(:), yout(:), zout(:), Bxout(:), Byout(:), Bzout(:), Gout(:), ...
    'VariableNames', {'ID', 'x position [m] ', 'y position [m]', 'z position [m]', 'Bx, speed in x [%c]', 'By, speed in y [%c]', 'Bz, speed in z [%c]', 'G, relativistic gamma factor'});

% Write the table to a CSV file
writetable(T, sprintf('output_%.2finfrom2ndcavity_%s.csv',outputlengthfrom2ndcav, masterfilename));