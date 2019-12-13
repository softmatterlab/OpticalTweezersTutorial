# Section 4.5 Statistical Physics

This folder contains the codes to analyze the data obtained from the double well potential experiment.
 



**Kramers transitions and stochastic resonance**

Run [figures_double_well.m](figures_double_well.m)  to obtain (a) Trajectory of a brownian bead in a bistable optical potential. A particle
with radius rp = 0:48 m in an aqueous solution at temperature T = 22C is subjected to the bistable
potential generated with two optical tweezers with orthogonal polarizations separated by a distance
d = 0:8 m along x􀀀 direction, whereas the separation between the equilibrium points, A and B, in
the resulting bistable potential is dAB = 0:71 m (see Fig. 22). The particle spend most of the time
near the equilibrium points A or B, but from time to time it is thermally activated to the saddle point
in S with a mean rate defined by Eqs. (36) and (38). (c-d) Recosntruction of the 2-D bistable optical
potential. (b)Reconstruction of the force field measured with FORMA (arrows) and of the potential
energy measured with the potential analysis (background colour). (c) x and y profiles of the optical
potential along the critical points. FORMA identifies two wells (full circles in A and B) and one saddle
(empty circle in B ) points, and measures their stiffness (dashed lines) along x-y-directions in (a). The
separation between equilibrium points along x􀀀 direction is dAB = 0:71 m. Probability per unit time
of the residence time in a bistable potential. (d) Probability for well A and (b) probability for well B, for
a brownian particle in the bistable optical potential shown in Fig. 22. Red curves show the expoential
probability function defined by Eq. (39), with the mean residence time estimated from the experimental
data:  A = 3:03 s and B = 0:67 s. The probability in wells A and B were estimated with 900 and
1000 crossing events correspondingly.


The data file is contained in tyhe following link: https://drive.google.com/open?id=13Q3KOchO9b2qjyguuXqS_kV-4P0tPaxd
[Determination_of_Zero_Shear_Viscosities_1stimage.m](Determination_of_Zero_Shear_Viscosities_1stimage.m) calculates the velocity autocorrelation of a particle trapped by optical tweezers in a 
simple fluid and determine viscosity from decay time. The files used for this script are [Water.txt](PassiveMicrorheologyData/Water.txt),  [PNP.txt](PassiveMicrorheologyData/PNP.txt) and  [CPyCl4mM.txt](PassiveMicrorheologyData/CPyCl4mM.txt)contained in **[PassiveMicrorheologyData](PassiveMicrorheologyData/)**. The first column of these files contains the particle position (a), while the second has the corresponding times. 



![alt text](https://github.com/LauraPerezG/tweezers_AOP_tutorial/blob/merge_26nov_ales_lau/sec_4_5_statistical_physics/double_well/figures/double_well.jpg
"Double- well potential")





***


 
This repository includes the analysis and plot codes for the the article **Optical Tweezers: A comprehennsive Tutorial  from Calibration to Applications** by *Jan Gieseler, Juan Ruben Gomez-Solano, [Alessandro Magazù](http://softmatterlab.org/people/alessandro-magazzu/), Isaac Castillo, [Laura Pérez García](http://softmatterlab.org/people/laura-perez-garcia/), Martha Gironella-Torrent, Xavier Viader-Godoy, Felix Ritort, Giusceppe Pesce, Alejandro V. Arzola, Karen Volke-Sepulveda and [Giovanni Volpe](http://softmatterlab.org/people/giovanni-volpe/)*. 
 
