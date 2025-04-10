for phase in 0.00 # 0.33        0.66        0.99        1.32        1.65        1.98        2.31        2.65        2.98      3.14   3.31        3.64        3.97         4.30        4.63        4.96        5.29        5.62        5.95        6.28
do
	gpt -o phia_simulationsEnergyMod_phi${phase}.gdf phia_simulationsEnergyMod_phi${phase}.in               
	gdf2his -o phia_simulationsEnergyMod_phi${phase}hist.gdf phia_simulationsEnergyMod_phi${phase}.gdf G 0.00001
	gdf2a -w 16 -o  phia_simulationsEnergyMod_phi${phase}hist.txt phia_simulationsEnergyMod_phi${phase}hist.gdf
	ex -s -c '1d3|x' phia_simulationsEnergyMod_phi${phase}hist.txt
done
