Simulating the phase offset between 2 cavities for tunable proton therapy.

Introductory materials about proton therapy and some of the equations used: https://gray.mgh.harvard.edu/media/com_dpattachments/attachments/com_content.article/Techniques-of-Proton-Radiotherapy-04-Basics.pdf 

Creating the GPT input file in phia_test_EMod_spreadBragg.m,
  relies on the E and H field and linaciris gdf files
  
Using sim_auto.bat to run GPT and write output files for all the different phase offset values. "bash sim_auto.bat"

More information can be found in the General Particle Tracer (GPT) manual: https://wiki.jlab.org/ciswiki/images/4/42/UserManual.pdf  

Energy distribution vs simulated particles in plots_automated.m or as Bragg curve (not scaled) in EnergyDensity_BraggCruve.m, takes 
'comparison_noE = readtable(sprintf('phia_simulationsEnergyMod_phi0.00_0.03Espread_nominalhist.txt'));' as a reference beam with RF into the cavities

Automated image creation (phantom water sample) pulling data from the histogram GPT output files created from phia_test_EMod_spreadBragg.m,
  'data = readtable(sprintf('phia_simulationsEnergyMod_phi%.2fhist.txt',phase));'
  
Image energy reconstruction from the output photos in ImageEnergyReconstruction.m, right now the lateral size of the beam in the images is set arbitrarily and off axis beams have not been tested.

Gif from braggpeakgif.m to show water phantom evolution as teh phase difference in the 2 cavity design changes (not all the necessary files may be on this git and you may need to go through the workflow for all phases to create the necessary inputs.
