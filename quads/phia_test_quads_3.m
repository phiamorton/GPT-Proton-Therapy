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
npos=5;
position2s= linspace(0.45, 0.8,npos);
quadstrengths1= linspace(15,25,30); %start with 25 to get focal length ~0.8 m %linspace(0.1,36,40);
quadstrengths2= linspace(15,35,10);
%quadstrengths3= linspace(15,35,10);
position3s= linspace(0.75, 1.2,npos);
%quadrupole strength in the unit of T/m~~~ dimension is IMPORTANT
counter=1;
set(0,'DefaultFigureWindowStyle','docked') 

for quadstrength1=1:length(quadstrengths1)
    for quadstrength2=1:length(quadstrengths2)
        for qps2=1:length(position2s)
            position3s= linspace(position2s(qps2)+length_quad/2, 1,npos);
            for quadstrength3=1:length(quadstrengths3)
                for qps3=1:length(position3s)
                        length_quad = 0.2062;
                        quadpos=[0.15,position2s(qps2),position3s(qps3)];
                        gq1 = -quadstrengths1(quadstrength1); %~36kG/m + focuses in x and - focuses in y 
                        %set first quad to focus~0.8m but sweeping strength
                        gq2 = quadstrengths2(quadstrength2); %add in second quad and adjust strength to look at min in x and y and minimize the diff in z between the two 
                        %set 3rd quad strength
                        gq3 = -quadstrengths3(quadstrength3);
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
                            
                        end 
                        
                
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
                
               
                        inputfiletext = [buildparticles; 
                            { ['quadrupole("wcs","z",' num2str(quadpos(1)) ',' num2str(length_quad) ',' num2str(gq1) ');'] }; 
                            { ['quadrupole("wcs","z",' num2str(quadpos(2)) ',' num2str(length_quad) ',' num2str(gq2) ');'] };
                            { ['quadrupole("wcs","z",' num2str(quadpos(3)) ',' num2str(length_quad) ',' num2str(gq3) ');'] };
                
                            linactext; {
                            ['tout(' num2str(0/c) ',' num2str((2)/beta0/c) ',' num2str((0.01)/beta0/c)  ');']; 
                            };];
                        
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
                    
                        area=pi*stdx.*stdy*1000*1000;
                        minarea=min(area);
                        indexmin=find(area==minarea);
                        z_for_min=avgz(indexmin);
                        if minarea<10 %only show configs with small area
                            figure('Visible','on');
                            scatter(avgz,stdx*1000, 'Color', "#0072BD", 'DisplayName', 'x')
                            hold on
                            scatter(avgz,stdy*1000, 'Color', "red", 'DisplayName', 'y')
                            xline(z_for_min,'DisplayName', 'Minimum area')
                            fill([quadpos(1)-length_quad/2, quadpos(1)+length_quad/2, quadpos(1)+length_quad/2, quadpos(1)-length_quad/2], [0, 0, yrms0*2*1000+2, yrms0*2*1000+2], 'b', 'FaceAlpha',0.1,'DisplayName', sprintf('quad position 1 at %.2f m, %.2f T/m ', quadpos(1),gq1),'LineStyle',"none")
                            fill([quadpos(2)-length_quad/2, quadpos(2)+length_quad/2, quadpos(2)+length_quad/2, quadpos(2)-length_quad/2], [0, 0, yrms0*2*1000+2, yrms0*2*1000+2], 'b', 'FaceAlpha',0.1,'DisplayName', sprintf('quad position 2 at %.2f m, %.2f T/m ', quadpos(2),gq2), 'LineStyle',"none")
                            fill([quadpos(3)-length_quad/2, quadpos(3)+length_quad/2, quadpos(3)+length_quad/2, quadpos(3)-length_quad/2], [0, 0, yrms0*2*1000+2, yrms0*2*1000+2], 'b', 'FaceAlpha',0.1,'DisplayName', sprintf('quad position 3 at %.2f m, %.2f T/m ', quadpos(3),gq3), 'LineStyle',"none")
                            ylim([0,yrms0*2*1000+2])
                            legend();
                            xlabel('Average Z [m]');
                            ylabel('Transverse Profile Size rms [mm]');
                            title(sprintf('Min area %.2f mm^2', minarea))
                            hold off
                            filename=sprintf('transverseprof_strength_%.2f_Tpm_Q1_and_%.2f_Tpm_Q2_at_%.2f_m_and_%.2f_Tpm_Q3_at_%.2f_m', quadstrengths1(quadstrength1),quadstrengths2(quadstrength2),qps2,quadstrengths3(quadstrength3),qps3);
                            saveas(gcf,filename,'png')
                        end
                end
            end
        end
    end
end