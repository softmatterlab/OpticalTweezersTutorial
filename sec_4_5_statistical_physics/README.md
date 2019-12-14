# Section 4.5 Statistical Physics

This folder contains the codes to analyze the data obtained from the double well potential experiment.
 



**Kramers transitions and stochastic resonance**

Run [figures_double_well.m](double_well/figures_double_well.m)  to obtain (a) Trajectory of a brownian bead in a bistable optical potential. A particle
with radius rp = 0:48 um in an aqueous solution at temperature T = 22C is subjected to the bistable
potential generated with two optical tweezers with orthogonal polarizations separated by a distance
d = 0.8 um along x- direction, whereas the separation between the equilibrium points, A and B, in
the resulting bistable potential is dAB = 0.71 um (see Fig. 22). The particle spend most of the time
near the equilibrium points A or B, but from time to time it is thermally activated to the saddle point
in S with a mean rate defined by Eqs. (36) and (38). (c-d) Recosntruction of the 2-D bistable optical
potential. (b)Reconstruction of the force field measured with FORMA (arrows) and of the potential
energy measured with the potential analysis (background colour). (c) x and y profiles of the optical
potential along the critical points. FORMA identifies two wells (full circles in A and B) and one saddle
(empty circle in B ) points, and measures their stiffness (dashed lines) along x-y-directions in (a). The
separation between equilibrium points along x-direction is dAB = 0:71 um. Probability per unit time
of the residence time in a bistable potential. (d) Probability for well A and (b) probability for well B, for
a brownian particle in the bistable optical potential shown in Fig. 22. Red curves show the expoential
probability function defined by Eq. (39), with the mean residence time estimated from the experimental
data: τA = 3.03 s and τB = 0.67 s. The probability in wells A and B were estimated with 900 and
1000 crossing events correspondingly.


The data file is contained in tyhe following link: https://drive.google.com/open?id=13Q3KOchO9b2qjyguuXqS_kV-4P0tPaxd




![alt text](https://github.com/LauraPerezG/tweezers_AOP_tutorial/blob/merge_26nov_ales_lau/sec_4_5_statistical_physics/double_well/figures/double_well.jpg
"Double- well potential")

**Fluctuation-dissipation relation for non-equilibrium steady states**

Run [NESSFDT.m](Fluctuation_dissipation_relation_for_NESS/NESSFDT.m) to obtain (a) Typical trajectories Θ(t) of a colloidal particle in a NESS, defined over the inverval [0,∞). The dashed line represents the mean drift, <θ(t)><sub>0</sub> = 2πjt+ const :, due to the non-zero probability current j induced by thenon-conservative term F. Inset: example of trajectory θ(t) defined over [0, 2π) (b) NESS probability
density function of θ, defined over [0, 2π) (bars) and reconstructed potential energy U(θ) (solid
line). The arrow indicates the direction of the non-conservative force f<sub>0</sub>, which shifts the maximum
of ρ<sub>NESS</sub>(θ) to the right relative to the minimum of U(θ). Inset: comparison between the NESS
distribution ρ<sub>NESS</sub>(θ) and the equilibrium one ρ<sub>eq</sub> (θ) (F = 0, sharp peak). (c) Correlation function between the observable Q(θ) = sin θ and the variable V (θ) = sin θ. (d) NESS correlation functions
involved in the integral form of the generalized FDT (47): C(0) - C(t) (dotted-dashed line), B(t )
(dashed line) and C(0) - C(t ) - B(t ) (solid line). Inset: estimate of the response function R(t ) by
means of the time derivative of [C(0) -C(t )]/(kBT ) (dotted-dashed line), and taking into account the
corrective term, [C(0) -C(t ) - B(t )]/(k<sub>B</sub>T ) (solid line).


The data file is contained in the data file [trajNESS_data.txt](data/trajNESS_data.txt).


![alt text](https://github.com/LauraPerezG/tweezers_AOP_tutorial/blob/merge_26nov_ales_lau/sec_4_5_statistical_physics/Fluctuation_dissipation_relation_for_NESS/figures/fluctuation_dissipation_NESS.jpg
"Fluctuation dissipation NESS")




***


 
This repository includes the analysis and plot codes for the the article **Optical Tweezers: A comprehennsive Tutorial  from Calibration to Applications** by *Jan Gieseler, Juan Ruben Gomez-Solano, [Alessandro Magazù](http://softmatterlab.org/people/alessandro-magazzu/), Isaac Castillo, [Laura Pérez García](http://softmatterlab.org/people/laura-perez-garcia/), Martha Gironella-Torrent, Xavier Viader-Godoy, Felix Ritort, Giusceppe Pesce, Alejandro V. Arzola, Karen Volke-Sepulveda and [Giovanni Volpe](http://softmatterlab.org/people/giovanni-volpe/)*. 
 
