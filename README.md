Simulating the phase offset between 2 cavities for tunable proton therapy.

Creating the GPT input file in phia_test_EMod_spreadBragg.m,
  relies on the E and H field and linaciris files
  
Using sim_auto.bat to run GPT and write output files for all the different phase offset values.

Energy distribution vs simulated particles in plots_automated.m or as Bragg curve (not scaled) in Bragg_curve.m

Automated image creation (phantom water sample) pulling data from the GPT output in EnergyDensity_BraggCurve.m,
  data = readtable(sprintf('phia_simulationsEnergyMod_phi%.2fhist.txt',phase));
  
Image energy reconstruction from the output photos in ImageEnergyReconstruction.m

Gif from braggpeakgif.m
