# Section 4.5 Statistical Physics

This folder contains the codes to analyze the data obtained from the double well potential experiment.
 



## Kramers transitions and stochastic resonance

Run [double_well_Fig27.m](double_well/programs/double_well_Fig27.m) in subfolder [double_well/programs/](double_well/programs/) to obtain (a) Trajectory of a brownian bead in a bistable optical potential. A particle
with radius rp = 0.48 um in an aqueous solution at temperature T = 22C is subjected to the bistable
potential generated with two optical tweezers with orthogonal polarizations separated by a distance
d = 0.8 um along x- direction, whereas the separation between the equilibrium points, A and B, in
the resulting bistable potential is dAB = 0.71 um (see Fig. 22). The particle spend most of the time
near the equilibrium points A or B, but from time to time it is thermally activated to the saddle point
in S with a mean rate defined by Eqs. (36) and (38). (c-d) Recosntruction of the 2-D bistable optical
potential. (b)Reconstruction of the force field measured with FORMA (arrows) and of the potential
energy measured with the potential analysis (background colour). (c) x and y profiles of the optical
potential along the critical points. FORMA identifies two wells (full circles in A and B) and one saddle
(empty circle in B ) points, and measures their stiffness (dashed lines) along x-y-directions in (a). The
separation between equilibrium points along x-direction is dAB = 0.71 um. Probability per unit time
of the residence time in a bistable potential. (d) Probability for well A and (b) probability for well B, for
a brownian particle in the bistable optical potential shown in Fig. 22. Red curves show the exponential
probability function defined by Eq. (39), with the mean residence time estimated from the experimental
data: τA = 3.03 s and τB = 0.67 s. The probability in wells A and B were estimated with 900 and
1000 crossing events correspondingly.


The data file is contained in tyhe following link: https://drive.google.com/open?id=13Q3KOchO9b2qjyguuXqS_kV-4P0tPaxd




![alt text](/sec_4_5_statistical_physics/double_well/figures/double_well.jpg
"Double- well potential")

## Fluctuation-dissipation relation for non-equilibrium steady states

Run [NESSFDT_Fig28.m](Fluctuation_dissipation_relation_for_NESS/programs/NESSFDT_Fig28.m) in subfolder  [Fluctuation_dissipation_relation_for_NESS/programs/](Fluctuation_dissipation_relation_for_NESS/programs) to obtain (a) Typical trajectories Θ(t) of a colloidal particle in a NESS, defined over the inverval [0,∞). The dashed line represents the mean drift, <θ(t)><sub>0</sub> = 2πjt+ const :, due to the non-zero probability current j induced by thenon-conservative term F. Inset: example of trajectory θ(t) defined over [0, 2π) (b) NESS probability
density function of θ, defined over [0, 2π) (bars) and reconstructed potential energy U(θ) (solid
line). The arrow indicates the direction of the non-conservative force f<sub>0</sub>, which shifts the maximum
of ρ<sub>NESS</sub>(θ) to the right relative to the minimum of U(θ). Inset: comparison between the NESS
distribution ρ<sub>NESS</sub>(θ) and the equilibrium one ρ<sub>eq</sub> (θ) (F = 0, sharp peak). (c) Correlation function between the observable Q(θ) = sin θ and the variable V (θ) = sin θ. (d) NESS correlation functions
involved in the integral form of the generalized FDT (47): C(0) - C(t) (dotted-dashed line), B(t )
(dashed line) and C(0) - C(t ) - B(t ) (solid line). Inset: estimate of the response function R(t ) by
means of the time derivative of [C(0) -C(t )]/(kBT ) (dotted-dashed line), and taking into account the
corrective term, [C(0) -C(t ) - B(t )]/(k<sub>B</sub>T ) (solid line).


The data file is contained in the data file [trajNESS_data.txt](Fluctuation_dissipation_relation_for_NESS/data/trajNESS_data.txt).

![alt text](/sec_4_5_statistical_physics/Fluctuation_dissipation_relation_for_NESS/figures/fluctuation_dissipation_NESS.jpg "Fluctuation-dissipation")

**Stochastic thermodynamics**


Run [stochastic_thermodynamics_Fig29.m](first_law_stochastic_thermodynamics/programs/stochastic_thermodynamics_Fig29.m)  in subfolder [first_law_stochastic_thermodynamics/programs/](first_law_stochastic_thermodynamics/programs/)  (a) Upper panel: Portion of the
trajectory x(t) (blue solid line) of a colloidal bead (2a = 2:73 um), dragged by moving
an optical tweezers (k = 1:15 pN μm^-1, v = 0:520 μms^-1) through water (T = 295K).
The dashed line represents the time evolution of the minimum of the trapping potential:
λ (t) = vt + λ(0). Lower panel: time evolution of λ(t) - x(t). (b) Dependence of the mean
work <W_τ (ο)>, mean heat dissipated into the aqueous medium  <Q_τ (□)>, and mean variation
of potential energy <ΔU_τ (◇)>,as a function of the measurement time τ. The solid line
represents the values given by equation (135). (c) Dependence of the standard deviation
of the work (ο), of the heat dissipated into the aqueous medium (□), and of the variation
of potential energy (◇) as a function of τ. The solid line represents the values given by the
square roots of equation (136). (d) Distribution of the work (in units of kBT) for different
values of the measurement time τ. From inner to outer curves: τ= 0:034 s, 0:069 s, 0:103 s,
0:138 s, 0:172 s, and 0:207 s. (e) Distribution of the dissipated heat (in units of kBT) for
different values τ. Same color code as in (d). Inset: verification of the integral fluctuation
theorem (141) for the total entropy production. (f) Distribution of the variation of potential
energy (in units of kBT) for different values τ. Same color code as in (d).

This program computes the work (Work), the heat (Heat), the
variation of the potential energy (Energy) and the total entropy production (totEntropy), 
measured over a time interval of duration tau, for a spherical particle (radius r) driven 
in water (viscosity gamma, temperature T0) by an optical trap (stiffness k) moving in  
a straightline at constant velocity v. The stochastic thermodynamics quantities are determined  
from the time (t) evolution of the particle position (x) and
the time evolution of the position of the optical trap (xtrap). 
To this end, first the numerical values of k and v are determined from x and xtrap,
respectively, and then they are used as input parameters within the formalism of Stochastic
Thermodynamics. The probability density functions of the stochastic
thermodynamic quantities are also computed, as well as their mean values
and standars deviations. Finally, the integral fluctuation theorem is
checked for the total entropy production at different measurement times
tau.

![alt text](/sec_4_5_statistical_physics/first_law_stochastic_thermodynamics/figures/stochastic_thermodynamics.jpg
"Stochastic thermodynamics")




***


 
This repository includes the analysis and plot codes for the the article **Optical Tweezers: A comprehensive Tutorial  from Calibration to Applications** by *[Jan Gieseler](https://scholar.google.com.ar/citations?user=6OKJlNgAAAAJ&hl=en), [Juan Ruben Gomez-Solano](https://www.fisica.unam.mx/es/personal.php?id=639), [Alessandro Magazù](http://softmatterlab.org/people/alessandro-magazzu/),  [Isaac Pérez Castillo](https://scholar.google.com.mx/citations?user=58GAc80AAAAJ&hl=en), [Laura Pérez García](http://softmatterlab.org/people/laura-perez-garcia/), [Martha Gironella-Torrent](https://scholar.google.com/citations?user=tITfJqkAAAAJ&hl=en), [Xavier Viader-Godoy](https://scholar.google.com/citations?user=dTLMJy0AAAAJ&hl=en), [Felix Ritort](http://ffn.ub.es/ritort/), [Giuseppe Pesce](https://scholar.google.com/citations?user=Sf4mmT8AAAAJ&hl=en), [Alejandro V. Arzola](https://orcid.org/0000-0002-4860-6330), [Karen Volke-Sepulveda](https://www.fisica.unam.mx/es/personal.php?id=27) and [Giovanni Volpe](http://softmatterlab.org/people/giovanni-volpe/)*.
