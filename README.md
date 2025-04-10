Simulating the phase offset between 2 cavities for tunable proton therapy.

Creating the GPT input file in phia_test_EMod_spreadBragg.m,
  relies on the E and H field and linaciris gdf files
  
Using sim_auto.bat to run GPT and write output files for all the different phase offset values. "bash sim_auto.bat"

Energy distribution vs simulated particles in plots_automated.m or as Bragg curve (not scaled) in Bragg_curve.m

Automated image creation (phantom water sample) pulling data from the GPT output files created from EnergyDensity_BraggCurve.m,
  'data = readtable(sprintf('phia_simulationsEnergyMod_phi%.2fhist.txt',phase));'
  
Image energy reconstruction from the output photos in ImageEnergyReconstruction.m, right now the lateral size of the beam in the images is set arbitrarily and off axis beams have not been tested.

Gif from braggpeakgif.m to show water phantom evolution as teh phase difference in the 2 cavity design changes (not all the necessary files may be on this git and you may need to go through the workflow for all phases to create the necessary inputs.
