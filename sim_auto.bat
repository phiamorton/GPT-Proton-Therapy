#/Applications/MATLAB_R2024b.app/bin/matlab -batch -nodisplay phia_test_EMod_SpreadBragg.m #creates input files
energyspreadpercent='0.03'
for phase in 0.00 # 0.33        0.66        0.99        1.32        1.65        1.98        2.31        2.65        2.98      3.14   3.31        3.64        3.97         4.30        4.63        4.96        5.29        5.62        5.95        6.28
do
	echo $energyspreadpercent 
	gpt -o phia_simulationsEnergyMod_phi${phase}_E${energyspreadpercent}.gdf phia_simulationsEnergyMod_phi${phase}_E${energyspreadpercent}.in               
	gdf2his -o phia_simulationsEnergyMod_phi${phase}_E${energyspreadpercent}hist.gdf phia_simulationsEnergyMod_phi${phase}_E${energyspreadpercent}.gdf G 0.00001
	gdf2a -w 16 -o  phia_simulationsEnergyMod_phi${phase}_E${energyspreadpercent}hist.txt phia_simulationsEnergyMod_phi${phase}_E${energyspreadpercent}hist.gdf
	ex -s -c '1d3|x' phia_simulationsEnergyMod_phi${phase}_E${energyspreadpercent}hist.txt
done
#/Applications/MATLAB_R2024b.app/bin/matlab -batch -nodisplay EnergyDensity_BraggCurve.m 
#/Applications/MATLAB_R2024b.app/bin/matlab -batch -nodisplay dose_calcs.m
#/Applications/MATLAB_R2024b.app/bin/matlab -batch -nodisplay ImageEnergyReconstruction.m

