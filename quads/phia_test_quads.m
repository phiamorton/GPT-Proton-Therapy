% runs through end of chopper slit - started 2/28/2022
% will include optional phase offset for section of linac
close all
clearvars

%%
%phioffsets = [0.00 1/4*3.14 3.14/2 3/4*3.14 3.14 3/2*3.14 6.28]; %3.4; %in rad, 0-2pi
phioffsets =  [0.00]; %[0.00  0.33        0.66        0.99        1.32        1.65        1.98        2.31        2.65        2.98      3.14   3.31        3.64        3.97         4.30        4.63        4.96        5.29        5.62        5.95        6.28];  %linspace(0, 2*pi,30)
energyspreadpercent= 0.03;
energy0 = 228.5; %alter energy into cavities
uniform=false;
NoRF=false;
ffac=false;
rounded = round(phioffsets,2);
%format bank
num2str(rounded);
%for pp = 1:length(phioffsets)
length_quad = 0.2062;
npos=20;
position2s=linspace(0.4, 0.7,npos);
quadstrengths=linspace(0.1,15,20);
     %quadrupole strength in the unit of T/m~~~ dimension is IMPORTANT
beamonitorpos=1; %m
tolerance = 0.05;  % Accept values within ±0.05 for the pos monitor
beammonitorarea=zeros([npos,length(quadstrengths)]);
divanglesy=zeros([npos,length(quadstrengths)]);
divanglesx=zeros([npos,length(quadstrengths)]);
for qps1=1:npos
    for quadstrength=1:length(quadstrengths)
        qps1;
        length_quad = 0.2062;
        quadpos=[0.15,position2s(qps1)];
        gq1 = -quadstrengths(quadstrength); %~36kG/m + focuses in x and - focuses in y
        gq2 = -gq1;
            %gq3 = 0.0001;
        phioffsetE = phioffsets;
        inputfilepath = 'output_';
        fieldpathname = '""';
        GPTpathname = 'C:\bin\'; 
        ffacE = 5.5;%*10; %5.5 ;%-482; %7.5; %5.1;
     
        if uniform ==true
            masterfilename = sprintf('output_EnergyMod_phi%.2f_E%.2f_Esp%.2f_uniform_quads', phioffsetE, energy0, energyspreadpercent);
        elseif uniform==false
            masterfilename = sprintf('output_EnergyMod_phi%.2f_E%.2f_Esp%.2f_quads', phioffsetE, energy0, energyspreadpercent);
        end
    
        if NoRF==true
            masterfilename= sprintf('output_noRF_EnergyMod_phi%.2f_E%.2f_Esp%.2f_uniform_quads', phioffsetE, energy0, energyspreadpercent);
            ffacE=0;
        end
        if ffac==true
            masterfilename= sprintf('output_EnergyMod_phi%.2f_E%.2f_Esp%.2f_ffac%.2f_quads', phioffsetE, energy0, energyspreadpercent,ffacE);
        end
        
        load energyMod_phase_1_24_2022_40cells30MeV
        philistE = philist;
        masterfilename;
        %% Define linac parameters
        
        freq = 2.856e9;
        dcellE = 14.7*0.0254; %distance between the cells, 14.7 inches, takes input as m
        a = 0.005; %0.5 cm
        ncellsE = 0; %length(philistE); %changed to 2 (only 2 cell cavities)
        phasebreakE = 2;
        subplotnum = 5;
        drift = .67; %m I think
        
        %% Define beam parameters
        
        npart0 = 2000;
        
        gamma0 = (energy0+938.27)/938.27; % 1.2435;
        
        dgamma0 = (energy0*energyspreadpercent/100+938.27)/938.27-1; % .03% energy spread
        
        mevion_25nA=true;
        mevion_1nA=false;
       
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
        %emit0 = .01e-6; % 3 pi mm-mrad emittance
        Qtot0 = 4.2e-15; %in C % assumes 6 uA pulsed average current
        %current for 2us period the pulse is there
        %mevion gave average current
        zposE0 = zlen0/1.8; %.03; %what is this doing
        sc = 0;
        xoffset=0; %m
        yoffset= 0; %m
        %tdiff = .001/beta0/c;
        
        %% Initialize particle distribution entering treatment room
            
        if uniform==false
            buildparticles = {
            'accuracy(6);';
            ['npart = ' num2str(npart0) ';'];
            ['sc = ' num2str(sc) ';'];
            % 'if(npart==1){';
            % ['setstartpar("beam",0,0,0,0,0,' num2str(gamma0*beta0) ',mp,-qe,' num2str(Qtot0) ');'];
            % '}';
            'if(npart > 1){';
            ['setparticles("beam",' num2str(npart0) ',mp,-qe,' num2str(Qtot0) ');'];
            ['setxdist("beam","g",0,' num2str(xrms0) ',3,3);'];
            ['setydist("beam","g",0,' num2str(yrms0) ',3,3);'];
            ['setzdist("beam","u", 0, ' num2str(zlen0) ');'];
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
        
        linactext = cell(2*ncellsE,1);
        zpos = zposE0;
        
        
        %should set up a loop for strength, length, position, do for 2 and 3
        %quads, find min divergence and size at desired distance away, 1-1.5 m
        %away
        %how it is doing length of quad (ie start and end)
        %inputfiletext=[{ ['quadrupole( "wcs","z",' num2str(zpos) ',' num2str(length_quad) ',' num2str(gq1) ');'] }];
        quadpos(2)
        inputfiletext = [buildparticles; 
            { ['quadrupole("wcs","z",' num2str(quadpos(1)) ',' num2str(length_quad) ',' num2str(gq1) ');'] }; 
            { ['quadrupole("wcs","z",' num2str(quadpos(2)) ',' num2str(length_quad) ',' num2str(gq2) ');'] };
            
            linactext; {
            ['tout(' num2str(0/c) ',' num2str((2)/beta0/c) ',' num2str((0.01)/beta0/c)  ');']; 
            };];
        %
        %'tout(' num2str((2*drift)/beta0/c) ');'
        %',' num2str((zpos+2)/beta0/c) ',' num2str((0.01)/beta0/c) 
        % Write input file
        masterfilenamein=sprintf('%s.in', masterfilename);
        fileID = fopen([masterfilenamein],'wt');
        for ii = 1:length(inputfiletext)
        fprintf(fileID,'%s \n',inputfiletext{ii});
        end
        fclose(fileID); 
        
        
        %run the GPT script
        system('bash "sim_auto.bat"');
        

        simavg = readtable(sprintf('avgfull_%s.txt',masterfilename));
        sprintf('avgfull_%s.txt',masterfilename)
        avg = table2struct(simavg,'ToScalar',true);
        times=avg.time;
        stdx=avg.stdx; %std dev in x in m
        stdy=avg.stdy;
        avgz=avg.avgz;
        
        fig=figure(qps1+quadstrength); 
        set(gcf, 'WindowStyle', 'docked');
        %figure('Visible', 'on');
        scatter(avgz,stdx*1000, 'Color', "#0072BD", 'DisplayName', 'x')
        hold on
        scatter(avgz,stdy*1000, 'Color', "red", 'DisplayName', 'y')
        %hold off
        %xline(quadpos(1),'-','DisplayName', sprintf('quad position 1 at %.2f m, %.2f T/m * %.2f m', quadpos(1),gq1,length_quad), 'LineWidth',2)
        %xline(quadpos(2),'-','DisplayName', sprintf('quad position 2 at %.2f m, %.2f T/m * %.2f m', quadpos(2),gq2,length_quad), 'LineWidth',2)
        fill([quadpos(1)-length_quad/2, quadpos(1)+length_quad/2, quadpos(1)+length_quad/2, quadpos(1)-length_quad/2], [0, 0, yrms0*1000+2, yrms0*1000+2], 'b', 'FaceAlpha',0.1,'DisplayName', sprintf('quad position 1 at %.2f m, %.2f T/m ', quadpos(1),gq1),'LineStyle',"none")
        fill([quadpos(2)-length_quad/2, quadpos(2)+length_quad/2, quadpos(2)+length_quad/2, quadpos(2)-length_quad/2], [0, 0, yrms0*1000+2, yrms0*1000+2], 'b', 'FaceAlpha',0.1,'DisplayName', sprintf('quad position 2 at %.2f m, %.2f T/m ', quadpos(2),gq2), 'LineStyle',"none")
        
        ylim([0,yrms0*1000+2])
        %xline(quadpos(3),'-','DisplayName', 'quad position 3', 'LineWidth',2)
        legend();
        xlabel('Average Z [m]');
        ylabel('Transverse Profile Size rms [mm]');
        %hold off
        title(sprintf('Transverse profile with %.0f quads',length(quadpos)), 'FontSize', 14);
        %saveas(gcf,sprintf('%sFODO.png', masterfilename))
        hold off

        % Find indices where z is within the tolerance of beamonitorpos
        indices = find(abs(avgz - beamonitorpos) <= tolerance);

        % Extract corresponding x values
        area_at_beamonitor = min(stdy(indices).*stdx(indices))*3.14;
        beammonitorarea(qps1,quadstrength)=area_at_beamonitor;
        divy_at_beamonitor=abs(stdy(max(indices)+1)-stdy(min(indices)-1))/(avgz(max(indices)+1)-avgz(min(indices)-1));
        divanglesy(qps1,quadstrength)=divy_at_beamonitor;
        divx_at_beamonitor=abs(stdx(max(indices)+1)-stdx(min(indices)-1))/(avgz(max(indices)+1)-avgz(min(indices)-1));
        divanglesx(qps1,quadstrength)=divx_at_beamonitor;
    end
    hold off
end



% Flatten the matrix and get linear indices
[x_flat, linear_indices] = sort(beammonitorarea(:));
[divy_flat, linear_indices_divy] = sort(divanglesy(:));
[divx_flat, linear_indices_divx] = sort(divanglesx(:));
% Get the top 5 smallest values and their indices
top_n = 2;
top_values = x_flat(1:top_n);
top_values_divy=divy_flat(1:top_n);
top_values_divx=divx_flat(1:top_n);
top_linear_indices = linear_indices(1:top_n);
top_linear_indices_divy = linear_indices_divy(1:top_n);
top_linear_indices_divx = linear_indices_divx(1:top_n);


% Convert linear indices to subscripts (row and column indices)
[row_indices, col_indices] = ind2sub(size(beammonitorarea), top_linear_indices);
[row_indices_divy, col_indices_divy] = ind2sub(size(divanglesy), top_linear_indices_divy);
[row_indices_divx, col_indices_divx] = ind2sub(size(divanglesx), top_linear_indices_divx);
% Get corresponding position2s and quadstrengths
corresponding_positions = position2s(row_indices);
corresponding_quadstrengths = quadstrengths(col_indices);
corresponding_positions_divy = position2s(row_indices_divy);
corresponding_quadstrengths_divy = quadstrengths(col_indices_divy);
corresponding_positions_divx = position2s(row_indices_divx);
corresponding_quadstrengths_divx = quadstrengths(col_indices_divx);

% Display results
disp('Top smallest beammonitorarea values and their corresponding position2s and quadstrengths:');
for i = 1:top_n
    fprintf('Value: %.2f mm, Position2s: %.2f m, QuadStrength: %.2f T/m\n', ...
        top_values(i)*1000, corresponding_positions(i), corresponding_quadstrengths(i));
    fprintf('Value: %.2f divergence in y, Position2s: %.2f m, QuadStrength: %.2f T/m\n', ...
        top_values_divy(i)*1000, corresponding_positions_divy(i), corresponding_quadstrengths_divy(i));
    fprintf('Value: %.2f divergence in x, Position2s: %.2f m, QuadStrength: %.2f T/m\n', ...
        top_values_divx(i)*1000, corresponding_positions_divx(i), corresponding_quadstrengths_divx(i));
end

colLabels = arrayfun(@(x) sprintf('Position= %.2f m', x), position2s, 'UniformOutput', false);
rowLabels = arrayfun(@(x) sprintf('Strength= %.2f T/m', x), quadstrengths, 'UniformOutput', false);

% Create the table
T = array2table(beammonitorarea, 'RowNames', rowLabels, 'VariableNames', colLabels)

T_divy = array2table(divanglesy, 'RowNames', rowLabels, 'VariableNames', colLabels)

T_divx = array2table(divanglesx, 'RowNames', rowLabels, 'VariableNames', colLabels)

h=figure()
h=heatmap(position2s, quadstrengths,  divanglesy)
h.Title = 'Divergence in y';
h.XLabel = 'Strength (T/m)';
h.YLabel = 'Position (m)';

h2=figure()
h2=heatmap(position2s, quadstrengths,  divanglesx)
h2.Title = 'Divergence in x';
h2.XLabel = 'Strength (T/m)';
h2.YLabel = 'Position (m)';

h3=figure()
h3=heatmap(position2s, quadstrengths,  beammonitorarea*1000*1000)
h3.Title = '~Beam Size [mm^2]';
h3.XLabel = 'Strength (T/m)';
h3.YLabel = 'Position (m)';


    %do the plotting for the quad magnetic field
%testing quad position
%data = readtable(sprintf('%s.txt',masterfilename));
% G = data.G;
% E=938.272*(G-1); %MeV
%     %G=G(~isnan(G));
%x=data.x;
%length(x)
% y=data.y;
% z=data.z; 
% Bx=data.fBx;
% By=data.fBy;
% B_tot=sqrt(By.^2+Bx.^2);
% figure;
% %quiver(x,y,Bx,By);
% %cb = colorbar;
% %cb.Label.String = 'Field Strength [T]';
% scatter(z,Bx)
% xlabel('x [m]');
% ylabel('y [m]');
% saveas(gcf,sprintf('%sQuadField.png', masterfilename))