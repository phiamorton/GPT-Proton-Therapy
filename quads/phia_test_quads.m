% % runs through end of chopper slit - started 2/28/2022
% % will include optional phase offset for section of linac
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
npos=10;
position2s= linspace(0.4, 0.8,npos);
quadstrengths1= [15,20,25]; %start with 25 to get focal length ~0.8 m %linspace(0.1,36,40);
quadstrengths2= linspace(10,25,15);
quadpos_forquadstrength_minarea=zeros(1,length(quadstrengths2));
%minareaforquadstrength=zeros(1,length(quadstrengths2));
     %quadrupole strength in the unit of T/m~~~ dimension is IMPORTANT
beamonitorpos=1; %m
tolerance = 0.005;  % Accept values within ±0.05 for the pos monitor
beammonitorarea=zeros([npos,length(quadstrengths2)]);
divanglesy=zeros([npos,length(quadstrengths2)]);
divanglesx=zeros([npos,length(quadstrengths2)]);
counter=1;
z_xmins=zeros(1,npos);
z_ymins=zeros(1,npos);
% area_at_xmins=zeros(1,npos);
area_at_ymins=zeros(1,npos);
z_focaldiffs=zeros(1,npos);
minbeamareas=zeros(1,npos);
minbeamareaspos=zeros(1,npos);
minbeamareasstrengths=zeros(1,length(quadstrengths2));
minbeamareasstrengths_table=zeros(length(quadstrengths1),length(quadstrengths2));
minbeamareasq2pos_table=zeros(length(quadstrengths1),length(quadstrengths2));
minfocaldiffs=zeros(1,length(quadstrengths2));
minfocaldiffs_table=zeros(length(quadstrengths1),length(quadstrengths2));

for quadstrength1=1:length(quadstrengths1)
    for quadstrength2=1:length(quadstrengths2)
        for qps2=1:npos
                qps2;
                length_quad = 0.2062;
                quadpos=[0.15,position2s(qps2)];
                gq1 = -quadstrengths1(quadstrength1); %~36kG/m + focuses in x and - focuses in y 
                %set first quad to focus~0.8m but sweeping strength
                gq2 = quadstrengths2(quadstrength2); %add in second quad and adjust strength to look at min in x and y and minimize the diff in z between the two 
                %pick 2 other initial magnet strengths and find 
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
                dcellE = 14.7*0.0254; %distance between the cells, 14.7 inches, given as m
                a = 0.005; %0.5 cm
                ncellsE = 0; %length(philistE); %changed to 0 for quads only
                phasebreakE = 2;
                subplotnum = 5;
                drift = .67; %m I think
               
        
                %% Define beam parameters
        
                npart0 = 2000;
        
                gamma0 = (energy0+938.27)/938.27; % 1.2435;
        
                dgamma0 = (energy0*energyspreadpercent/100+938.27)/938.27-1; % .03% energy spread
                t_bunch= 2*10^(-6); %2 us
                c = 2.998e8; %m/s
                beta0= sqrt(1-1/(gamma0^2));
        
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
                isocenter = zposE0+dcellE/2; %m
                %tdiff = .001/beta0/c;
                
        
                mevion_25nA=true;
                mevion_1nA=false;
        
                if mevion_1nA==true
                    divangx0 = (3.794-3.496)/120; %.58; %change in x [mm] over 12 cm
                    divangy0 = (4.299-4.007)/120; % .67;
                    xiso=3.604/1000 ;%m 4.9mm
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
                    xiso=5.05/1000; %m 4.9mm
                    yiso = 6.23/1000 ;%m 6mm
                    divangx0*isocenter;
                    divangy0*isocenter;
                    xrms0 = xiso-divangx0*isocenter;
                    yrms0= yiso-divangy0*isocenter;
                    %remember, div in GPT MUST be rad/m ie div/sizerms
                    %this actually gives isocenter yrms0=6.28 and xrms0=5.09
                   
                    %based on mevion numbers this will be ~3-4mm at 1nA or 5-6mm at 25 nA
                    
                    %divergence of beam
                    
                end 
                %beta0 = .5944; %v/c? 
                
        
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
                    ['addxdiv("beam",0,' num2str(divangx0/xiso) ');'];
                    ['addydiv("beam",0,' num2str(divangy0/yiso) ');'];
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
                quadpos(2);
                inputfiletext = [buildparticles; 
                    { ['quadrupole("wcs","z",' num2str(quadpos(1)) ',' num2str(length_quad) ',' num2str(gq1) ');'] }; 
                    { ['quadrupole("wcs","z",' num2str(quadpos(2)) ',' num2str(length_quad) ',' num2str(gq2) ');'] };
        
                    linactext; {
                    ['tout(' num2str(0/c) ',' num2str((2)/beta0/c) ',' num2str((0.01)/beta0/c)  ');']; 
                    };];
                %['tout(' num2str(0/c) ',' num2str((2)/beta0/c) ',' num2str((0.01)/beta0/c)  ');']; 
                %['tout(' num2str((isocenter-.06)/beta0/c) ',' num2str((isocenter+.06)/beta0/c) ',' num2str((.03)/beta0/c)  ');'];  
                %use line above if checking size at isocenter
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
                counter=counter+1;
                %fig=figure(counter); 
                %qps2+quadstrength2
                %set(gcf, 'WindowStyle', 'docked');
                % figure('Visible', 'off');
                % scatter(avgz,stdx*1000, 'Color', "#0072BD", 'DisplayName', 'x')
                % hold on
                % scatter(avgz,stdy*1000, 'Color', "red", 'DisplayName', 'y')
                % %hold off
                % %xline(quadpos(1),'-','DisplayName', sprintf('quad position 1 at %.2f m, %.2f T/m * %.2f m', quadpos(1),gq1,length_quad), 'LineWidth',2)
                % %xline(quadpos(2),'-','DisplayName', sprintf('quad position 2 at %.2f m, %.2f T/m * %.2f m', quadpos(2),gq2,length_quad), 'LineWidth',2)
                % fill([quadpos(1)-length_quad/2, quadpos(1)+length_quad/2, quadpos(1)+length_quad/2, quadpos(1)-length_quad/2], [0, 0, yrms0*2*1000+2, yrms0*2*1000+2], 'b', 'FaceAlpha',0.1,'DisplayName', sprintf('quad position 1 at %.2f m, %.2f T/m ', quadpos(1),gq1),'LineStyle',"none")
                % fill([quadpos(2)-length_quad/2, quadpos(2)+length_quad/2, quadpos(2)+length_quad/2, quadpos(2)-length_quad/2], [0, 0, yrms0*2*1000+2, yrms0*2*1000+2], 'b', 'FaceAlpha',0.1,'DisplayName', sprintf('quad position 2 at %.2f m, %.2f T/m ', quadpos(2),gq2), 'LineStyle',"none")
                % 
                % ylim([0,yrms0*2*1000+2])
                % %xline(quadpos(3),'-','DisplayName', 'quad position 3', 'LineWidth',2)
                % legend();
                % xlabel('Average Z [m]');
                % ylabel('Transverse Profile Size rms [mm]');
                % %hold off
                % title(sprintf('Transverse profile with %.0f quads',length(quadpos)), 'FontSize', 14);
                %saveas(gcf,sprintf('%sFODO.png', masterfilename))
                %hold off
        
                % Find indices where size in x,y is min
                %stdx(round(length(stdx)/2):length(stdx))
                %index = find(array >= quadpos2, 1);
                index_zpastquad=min(find(avgz>=(position2s(qps2)+length_quad/2)));
                avgz(index_zpastquad);
                min_stdx_past1m=min(stdx(index_zpastquad:length(stdx)));
                min_stdy_past1m=min(stdy(index_zpastquad:length(stdy)));

                indicesx = find(stdx==min_stdx_past1m); %need to adapt for pos>1m
                indicesy= find(stdy==min_stdy_past1m);
                area_at_xmin = min(stdy(indicesx).*stdx(indicesx))*3.14;
                area_at_ymin = min(stdy(indicesy).*stdx(indicesy))*3.14;
                
                if area_at_ymin < area_at_xmin
                    minbeamarea=area_at_ymin;
                    minbeamareaspos(qps2)=avgz(indicesx);
                end
                if area_at_xmin <= area_at_ymin
                    minbeamarea=area_at_xmin;
                    minbeamareaspos(qps2)=avgz(indicesy);
                end
                minbeamareas(qps2)=minbeamarea;
                %create minbeamareas array and store min area
                %then plot (for each strength) position vs min area
                %for each strength find overall min and plot quad strength
                %vs global min 

                z_xmin=avgz(indicesx);
                z_ymin=avgz(indicesy);
                z_focaldiff=abs(z_xmin-z_ymin);
                z_focaldiffs(qps2)=z_focaldiff;

                % Extract corresponding x values
                
                z_xmins(qps2)=z_xmin;
                z_ymins(qps2)=z_ymin;

                area_at_xmins(qps2)=area_at_xmin;

                area_at_ymins(qps2)=area_at_ymin;

                %beammonitorarea(qps2,quadstrength2)=area_at_xmin;
                % divy_at_beamonitor=abs(stdy(max(indices)+1)-stdy(min(indices)-1))/(avgz(max(indices)+1)-avgz(min(indices)-1));
                % divanglesy(qps2,quadstrength2)=divy_at_beamonitor;
                % divx_at_beamonitor=abs(stdx(max(indices)+1)-stdx(min(indices)-1))/(avgz(max(indices)+1)-avgz(min(indices)-1));
                % divanglesx(qps2,quadstrength2)=divx_at_beamonitor;
            
        end
        % figure()
        % plot(position2s,minbeamareas)

        indices_zdiffmin=find(z_focaldiffs==(min(z_focaldiffs)));
        if length(indices_zdiffmin)>1
            indices_zdiffmin=indices_zdiffmin(1)
        end
        %quadpos_forquadstrength_minarea(quadstrength2);
        position2s(indices_zdiffmin);
        quadpos_forquadstrength_minarea(quadstrength2)=position2s(indices_zdiffmin);
        minfocaldiffs(quadstrength2)=min(z_focaldiffs);
        % minareaforquadstrength(quadstrength2)=(area_at_ymins(indices_zdiffmin)*1000*1000+area_at_xmins(indices_zdiffmin)*1000*1000)/2
        sprintf('min z diff =%.2f m with quad 2 position %.2f m and quads strengths %.2f T/m (Q1) %.2f T/m (Q2) and minimum area ~%.4f mm^2 at position %.2f m ',min(z_focaldiffs), position2s(indices_zdiffmin), quadstrengths1(quadstrength1), quadstrengths2(indices_zdiffmin), minbeamareas(indices_zdiffmin)*1000*1000,minbeamareaspos(indices_zdiffmin))
        figure('Visible','on')
        yyaxis left 
        plot(position2s,minbeamareas*1000*1000)
        % area_at_xmins;
        hold on
        % plot(position2s,area_at_ymins*1000*1000);
        xlabel('quad 2 position [m]')
        ylabel('minimum area [mm^2]')
        title(sprintf('area for quad 1 pos %.2f quad 2 pos %.2f for strength %.2f T/m and %.2f T/m', quadpos(1),quadpos(2), quadstrengths2(quadstrength2),quadstrengths1(quadstrength1)))
        yyaxis right
        %plot(position2s, z_xmins)
        %hold on
        plot(position2s, z_focaldiffs)
        %xlabel('quad 2 position [m]')
        ylabel('difference in x and y focal point [m]')
        title(sprintf('focal point difference for strength %.2f T/m (Q1) and %.2f T/m (Q2)', quadstrengths1(quadstrength1),quadstrengths2(quadstrength2)))
        hold off
        minbeamareasstrengths(quadstrength2)=(minbeamareas(indices_zdiffmin));
        minbeamareasstrengths_table(quadstrength1,quadstrength2)=(minbeamareas(indices_zdiffmin));
        minbeamareasq2pos_table(quadstrength1,quadstrength2)=(position2s(indices_zdiffmin));
        minfocaldiffs_table(quadstrength1,quadstrength2)=min(z_focaldiffs);
    end
    hold off
    figure(quadstrength1)
    yyaxis left
    plot(quadstrengths2, quadpos_forquadstrength_minarea, 'DisplayName','Optimal quad position')
    hold on
    plot(quadstrengths2, minfocaldiffs, 'DisplayName', 'minimum focal difference')
    title('Quad 2 position for minimum difference in focal points')
    title(sprintf('Q2 position for strength %.2f T/m (Q1) and %.2f T/m (Q2)', quadstrengths1(quadstrength1),quadstrengths2(quadstrength2)))
    xlabel('Quad Strength [T/m]')
    ylabel('Optimal Position/Focal Difference [m]')
    
    yyaxis right
    plot(quadstrengths2, minbeamareasstrengths*1000*1000,"DisplayName",'Min beam area')
    title('Quad 2 strength vs approx area at focal point')
    %xlabel('Quad Strength [T/m]')
    ylabel('Minimum area [mm^2]')
    legend()
    %legend()
    hold off

end
% 
% rowLabels = arrayfun(@(x) sprintf('Strength= %.2f T/m', x), quadstrengths1, 'UniformOutput', false);
% colLabels = arrayfun(@(x) sprintf('Strength= %.2f T/m', x), quadstrengths2, 'UniformOutput', false);
% % 
% % % Create the table
% T = array2table(minbeamareasstrengths_table*1000*1000, 'RowNames', rowLabels, 'VariableNames', colLabels)
h=figure()
h=heatmap(quadstrengths1,quadstrengths2, minbeamareasstrengths_table*1000*1000)
h.Title = 'Min area for quad strengths';
h.XLabel = 'Strength Q1 (T/m)';
h.YLabel = 'Strength Q2 (T/m)';

h=figure()
h=heatmap(quadstrengths1,quadstrengths2, minbeamareasq2pos_table)
h.Title = 'Q2 position corresponding to min area';
h.XLabel = 'Strength Q1 (T/m)';
h.YLabel = 'Strength Q2 (T/m)';

h=figure()
h=heatmap(quadstrengths1,quadstrengths2, minfocaldiffs_table)
h.Title = 'Minimum focal point difference [m]';
h.XLabel = 'Strength Q1 (T/m)';
h.YLabel = 'Strength Q2 (T/m)';

% % Flatten the matrix and get linear indices
% [x_flat, linear_indices] = sort(beammonitorarea(:));
% [divy_flat, linear_indices_divy] = sort(divanglesy(:));
% [divx_flat, linear_indices_divx] = sort(divanglesx(:));
% % Get the top 5 smallest values and their indices
% top_n = 2;
% top_values = x_flat(1:top_n);
% % top_values_divy=divy_flat(1:top_n);
% % top_values_divx=divx_flat(1:top_n);
% top_linear_indices = linear_indices(1:top_n);
% % top_linear_indices_divy = linear_indices_divy(1:top_n);
% % top_linear_indices_divx = linear_indices_divx(1:top_n);
% 
% 
% % Convert linear indices to subscripts (row and column indices)
% [row_indices, col_indices] = ind2sub(size(beammonitorarea), top_linear_indices);
% % [row_indices_divy, col_indices_divy] = ind2sub(size(divanglesy), top_linear_indices_divy);
% % [row_indices_divx, col_indices_divx] = ind2sub(size(divanglesx), top_linear_indices_divx);
% % Get corresponding position2s and quadstrengths
% corresponding_positions = position2s(row_indices);
% corresponding_quadstrengths = quadstrengths2(col_indices);
% % corresponding_positions_divy = position2s(row_indices_divy);
% % corresponding_quadstrengths_divy = quadstrengths2(col_indices_divy);
% % corresponding_positions_divx = position2s(row_indices_divx);
% % corresponding_quadstrengths_divx = quadstrengths2(col_indices_divx);
% 
% % Display results
% disp('Top smallest beammonitorarea values and their corresponding position2s and quadstrengths:');
% for i = 1:top_n
%     fprintf('Value: %.2f mm, Position2s: %.2f m, QuadStrength: %.2f T/m\n', ...
%         top_values(i)*1000, corresponding_positions(i), corresponding_quadstrengths(i));
%     % fprintf('Value: %.2f divergence in y, Position2s: %.2f m, QuadStrength: %.2f T/m\n', ...
%     %     top_values_divy(i)*1000, corresponding_positions_divy(i), corresponding_quadstrengths_divy(i));
%     % fprintf('Value: %.2f divergence in x, Position2s: %.2f m, QuadStrength: %.2f T/m\n', ...
%     %     top_values_divx(i)*1000, corresponding_positions_divx(i), corresponding_quadstrengths_divx(i));
% end
% 
% rowLabels = arrayfun(@(x) sprintf('Position= %.2f m', x), position2s, 'UniformOutput', false);
% colLabels = arrayfun(@(x) sprintf('Strength= %.2f T/m', x), quadstrengths2, 'UniformOutput', false);
% 
% % Create the table
% T = array2table(beammonitorarea, 'RowNames', rowLabels, 'VariableNames', colLabels)
% % 
% % T_divy = array2table(divanglesy, 'RowNames', rowLabels, 'VariableNames', colLabels)
% % 
% % T_divx = array2table(divanglesx, 'RowNames', rowLabels, 'VariableNames', colLabels)
% % 
% % h=figure()
% % h=heatmap(quadstrengths2, position2s,   divanglesy)
% % h.Title = 'Divergence in y';
% % h.XLabel = 'Strength (T/m)';
% % h.YLabel = 'Position (m)';
% % 
% % h2=figure()
% % h2=heatmap(quadstrengths2, position2s,   divanglesx)
% % h2.Title = 'Divergence in x';
% % h2.XLabel = 'Strength (T/m)';
% % h2.YLabel = 'Position (m)';
% 
% h3=figure()
% h3=heatmap(quadstrengths2, position2s,  beammonitorarea*1000*1000)
% h3.Title = '~Beam Size [mm^2]';
% h3.XLabel = 'Strength (T/m)';
% h3.YLabel = 'Position (m)';


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
