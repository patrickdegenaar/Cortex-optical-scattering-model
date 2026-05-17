# Cortex-optical-scattering-model
Monte Carlo model of how light scatters through the brain

For many reasons, but particularly optogenetics, we need to have a model of how light traverses through brain tissue in order to determine the irradiance at a given point. This model uses a Monte Carlo approach - simulating photon by photon to build up an irradiance profile.

It has 3 stages:
1. Main_1_settingsGenerator: This creates the various variables required for the Monte Carlo simulation. Particularly, it creates scattering mean free path and angular distribution look-up tables to speed up the simulation.
2. Main_2_scatterSimulation: This performs the actual Monte Carlo simulation, randomising an initial emission direction (from the emission look-up table distribution), and then at the mean free path, it works out the scattering direction based on Mie and Rayleigh scattering behaviours, and absorption. This process repeats for a defined number of photons within a defined 3D matrix of a given resolution
3. Main_3_analysis: This performs actual analysis - using a 3D Gaussian to smooth the Monte Carlo result and then normalising. This normalised distribution can then be multiplied by whatever the emission of the light source is.

light sources include: 
LED (Lambertian)
collimated LED (Tyndall, pseudo-collimated)
Light fibre (collimated)

Please email me: patrick.degenaar@newcastle.ac.uk if you have any comments

Publications to reference:
Dong N, Berlinguer-Palmini R, Soltan A, Ponon N, O'Neil A, Travelyan A, Maaskant P, Degenaar P, Sun X. Opto-electro-thermal optimization of photonic probes for optogenetic neural stimulation. J Biophotonics. 2018 Oct;11(10):e201700358. doi: 10.1002/jbio.201700358. Epub 2018 Jul 16. PMID: 29603666.

N. Dong et al., "Optogenetic Multiphysical Fields Coupling Model for Implantable Neuroprosthetic Probes," in IEEE Access, vol. 12, pp. 129160-129172, 2024, doi: 10.1109/ACCESS.2024.3441571.
keywords: {Probes;Optogenetics;Biomedical optical imaging;Absorption;Neuroprostheses;Optical fibers;Stimulated emission;Neural engineering;Optogenetics;prosthetic brain implants;multi-physical fields coupling model;tissue optics;bioheat transfer;opto-neuro interaction},


https://pubmed.ncbi.nlm.nih.gov/29603666/
