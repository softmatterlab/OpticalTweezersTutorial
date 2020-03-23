# Section 5 Optical tweezers in vacuum

This folder contains the codes to analyze the data obtained from the experiments with optical tweezers in vacuum.



**From overdamped to underdamped**


Run [overdamped_to_uderdamped_Fig35.m](overdamped_to_underdamped/programs/overdamped_to_uderdamped_Fig35.m)
We plot the measured PSD for the z-axis
motion, normalized to the damping rate 
γ_0/(2π), at three different pressures for a levitated
particle (radius a = 68 nm, laser power P≋150mW) that is overdamped (blue line),
critically damped (orange line), and underdamped (green line). The dashed black lines are
least-square fits to equation (164) (1000 mbar) or equation (163) (60 and 2.5 mbar). The
colored vertical solid lines indicate the roll-off frequencies 
Ω_c/2π. For the time traces at
60 and 2.5 mbar, the spectra contain leakage signals from the other oscillation axes above
100 kHz. Data from E. Hebestreit, “Thermal properties of levitated nanoparticles,” Ph.D. thesis, ETH Zurich (2017).

![alt text](/sec_5_optical_tweezers_vacuum/overdamped_to_underdamped/figures/Fig35.jpg 
"overdamped to underdamped ")



**Gas composition**

Run [composition_gas_Fig39.m](gas_composition/programs/composition_gas_Fig39.m)Gas composition in a vaccum chamber. Gas composition in the vacuum chamber
at pressures below 10􀀀3 mbar plotted over time. The solid lines show the partial pressures of
different gas species measured with a residual gas analyzer (RGA). The sum of all the partial
pressures yields the total pressure at the RGA (dashed line). Due to the configuration of the
vacuum system and the reduced pumping speed at the RGA, the pressure at the RGA deviates
from the pressure in the main vacuum chamber (dash-dot line). The initial rise in partial
pressures after turning on the RGA is attributed to the warm-up process and desorption of
gases from the filament of the RGA. Reproduced from Ref. [467].

![alt text](/sec_5_optical_tweezers_vacuum/gas_composition/figures/pressure.jpg "Gas composition")
***
 
 
 **Particle loading**
 
Run [particle_loading_Fig40.m](particle_loading/programs/particle_loading_Fig40.m) Particle loading with a nebulizer. The particle is loaded by spraying a solution
of nanoparticles through a nozzle which is placed above the focus. (a) Positioning of the
nozzle in the vacuum chamber. (b) Nozzle to funnel the falling particles towards the focus
of the trapping laser. (c) Histogram of brightness observed with a camera from the side.
The inset shows the brightness over a wider range. Reproduced from J. Gieseler, “Dynamics of optically levitated nanoparticles in high vacuum,” Ph.D. thesis, Universitat Politècnica de
Catalunya (2014). (a,b) and
adapted from F. Ricci, “Levitodynamics toward force nano-sensors in vacuum,” Ph.D. thesis, Universitat Politècnica de Catalunya
(2019).

![alt text](particle_loading/figures/particle_loading.jpg "Particle loading")

 
This repository includes the analysis and plot codes for the the article **Optical Tweezers: A comprehennsive Tutorial  from Calibration to Applications** by *Jan Gieseler, Juan Ruben Gomez-Solano, [Alessandro Magazù](http://softmatterlab.org/people/alessandro-magazzu/),  Isaac Pérez Castillo [Laura Pérez García](http://softmatterlab.org/people/laura-perez-garcia/), Martha Gironella-Torrent, Xavier Viader-Godoy, Felix Ritort, Giusceppe Pesce, Alejandro V. Arzola, Karen Volke-Sepulveda and [Giovanni Volpe](http://softmatterlab.org/people/giovanni-volpe/)*. 
 
