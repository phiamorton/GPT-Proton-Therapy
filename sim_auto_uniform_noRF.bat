energyspreadpercent='0.03'
energy='228.50'
for phase in 0.00 # 0.33        0.66        0.99        1.32        1.65        1.98        2.31        2.65        2.98      3.14   3.31        3.64        3.97         4.30        4.63        4.96        5.29        5.62        5.95        6.28
do
	echo $energyspreadpercent
	filepath=output_noRF_EnergyMod_phi${phase}_E${energy}_Esp${energyspreadpercent}_uniform 
	gpt -o ${filepath}.gdf ${filepath}.in               
	gdf2his -o ${filepath}hist.gdf ${filepath}.gdf G 0.00001
	gdf2a -w 16 -o ${filepath}.txt ${filepath}hist.gdf
	ex -s -c '1d3|x' ${filepath}.txt
done

