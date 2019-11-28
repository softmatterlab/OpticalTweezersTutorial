# Section 3: Calibration
 This folder contains the codes to implement the different calibration methods discussed in the text.
 
 
**Potential and Equipartition Method**

Run potential_analysis.m to obtain the values of the stiffness given by the potential method. This will call the functions, [prob_dist_energy.m](prob_dist_energy.m),  pot_lfit.m and pot_nlfit.m contained in the folder **pot**.
Run eq_analysis.m to obtain the values of the stiffness given by the equipartition method. This will call the function eq1d.m contained in the folder **eq**.
To obtain the corresponding figure of the fitting run pot_figs_table. 


**Mean Square Displacement**


**Autocorrelation Function**

**Power Spectrum Density**


**FORMA**

**Bayesian**
  


| Method        | Input variables           | Output variables |
| ------------- |-------------| -------------|
|Drift method  | ... |  ...|
| Potential     |  x: Time series of the position <br> T:temperature  <br>  n_b: number of bins| x_bins: position of the histogram bins <br> mU: mean potential <br> EU: standar daviation of the potential <br> k_pot: stiffnes <br> Ek_pot:standar deviation of the stiffness <br> m_hist: mean histogram <br> E_hist: standard deviation of the histogram <br> h_0: ? <br> x_eq : equilibrium position |
| Equipotential   | x: Time series of the position <br> T:temperature  <br>  deltax: error in the position detection |  k_eq:stiffness <br> Ek_eq: standard deviation of the stiffness|
| Mean square displacement | x: Time series of the position <br> T:temperature  <br>  deltat: time between acquired frames <br> maxlag: lag integer |k_msd:stiffness <br> Ek_msd:standard deviation of the stiffness <br> D_msd: diffussion <br> ED_msd <br> tau:lag times <br> mmsd: mean mean squared displacement <br> Emsd: standard deviation of the mean square displacement <br> indc:? |
| Autocorrelation function | ... |  ...|
| PSD | ... |  ...|
|FORMA  | ... |  ...|


 
This repository includes the analysis and plot codes for the the article **Optical Tweezers: A comprehennsive Tutorial  from Calibration to Applications** by *Jan Gieseler, Juan Ruben Gomez-Solano, Alessandro Magazù, Isaac Castillo, Laura Pérez García, Martha Gironella-Torrent, Xavier Viader-Godoy, Felix Ritort, Giusceppe Pesce, Alejandro V. Arzola, Karen Volke-Sepulveda and Giovanni Volpe*. 
 
