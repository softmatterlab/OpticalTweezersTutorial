# Section 3: Calibration
 This folder contains the codes to implement the different calibration methods discussed in the text.
 
How to use the calibration section:

**Potential and Equipartition Method**

Run [pot_analysis.m](pot_analysis.m) to obtain the values of the stiffness given by the potential method. This will call the functions, [prob_dist_energy.m](pot/prob_dist_energy.m), [pot_lfit.m](pot/pot_lfit.m)  and   [pot_nlfit.m](pot/pot_nlfit.m) contained in the folder **[pot](pot/)**.
Run eq_analysis.m to obtain the values of the stiffness given by the equipartition method. This will call the function [eq1d.m](eq/eq1d.m) contained in the folder **[eq](eq/)**.
To obtain the corresponding figure of the fitting run [pot_figs_table.m](pot_figs_table.m), this will call function  [plotsub_pot.m](pot/plotsub_pot.m). 


**Mean Square Displacement**

Run [msd_analysis.m](msd_analysis.m) to obtain the values of the stiffness and difussion coeficient given by the mean square displacement analysis method. This will call the functions, [msd_nlfit.m](msd/msd_nlfit.m) contained in the folder **[msd](msd/)**.


To obtain the corresponding figure of the fitting run [msd_figs_table.m](msd_figs_table.m), this will call function  [plotsub_msd.m](msd/plotsub_msd.m).

**Autocorrelation Function**

Run [acf_analysis.m](acf_analysis.m) to obtain the values of the stiffness and difussion coeficient given by the mean square displacement analysis method. This will call the functions, [acf_nlfit.m](acf/acf_nlfit) and [acf_lfit.m](acf/acf_lfit) contained in the folder  **[acf](acf/)**.


To obtain the corresponding figure of the fitting run [acf_figs_table.m](msd_figs_table.m), this will call function  [plotsub_acf.m](acf/plotsub_acf.m).


**Power Spectrum Density**

Run [psd_analysis.m](psd_analysis.m) to obtain the values of the stiffness and difussion coeficient given by the mean square displacement analysis method. This will call the functions, [psd_nlfit.m](acf/acf_nlfit) and [acf_lfit.m](acf/acf_lfit) contained in the folder  **[psd](psd/)**.


To obtain the corresponding figure of the fitting run [acf_figs_table.m](msd_figs_table.m), this will call function  [plotsub_acf.m](acf/plotsub_acf.m).

**FORMA**

**Bayesian**
  
**WLSICE, non linear fit of correlated data**

The functions [msd_nlfit.m](msd/msd_nlfit) will call function [wlsice.m](wlsice/wlsice.m) to perform a non-linear fitting considering the whole covariance matrix.  The target function that [wlsice.m](wlsice/wlsice.m) uses depends on the fitting method. This is especified by the 





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
 
