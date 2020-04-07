# Section 3: Calibration
 This folder contains the codes to implement the different calibration methods discussed in the text.
 
How to use the calibration section:

## Potential and Equipartition Method

Run [pot_analysis.m](pot/pot_analysis.m) to obtain the values of the stiffness given by the potential method. This will call the functions, [prob_dist_energy.m](pot/prob_dist_energy.m), [pot_lfit.m](pot/pot_lfit.m)  and   [pot_nlfit.m](pot/pot_nlfit.m) contained in the folder **[pot](pot/)**.


Run eq_analysis.m to obtain the values of the stiffness given by the equipartition method. This will call the function [eq1d.m](eq/eq1d.m) contained in the folder **[eq](eq/)**.

To obtain the corresponding figure of the fitting run [pot_figs_table.m](pot/pot_figs_table.m), this will call function  [plotsub_pot.m](pot/plotsub_pot.m) which will show the fit of each experiment in a subplot. 

![alt text](https://github.com/LauraPerezG/tweezers_AOP_tutorial/blob/merge_26nov_ales_lau/sec3_calibration/figures/POT.jpg "Fit for Potential and Equipartition method")

## Mean Square Displacement

Run [msd_analysis.m](msd/msd_analysis.m) to obtain the values of the stiffness and difussion coeficient given by the mean square displacement analysis method. This will call the functions, [msd_nlfit.m](msd/msd_nlfit.m) contained in the folder **[msd](msd/)**.


To obtain the corresponding figure of the fitting run [msd_figs_table.m](msd/msd_figs_table.m), this will call function  [plotsub_msd.m](msd/plotsub_msd.m) which will show the fit of each experiment in a subplot. 

![alt text](https://github.com/LauraPerezG/tweezers_AOP_tutorial/blob/merge_26nov_ales_lau/sec3_calibration/figures/MSD.jpg "Fit for Mean Square Displacement")


## Autocorrelation Function

Run [acf_analysis.m](acf/acf_analysis.m) to obtain the values of the stiffness and difussion coeficient given by the mean square displacement analysis method. This will call the functions, [acf_nlfit.m](acf/acf_nlfit.m) and [acf_lfit.m](acf/acf_lfit.m) contained in the folder  **[acf](acf/)**.


To obtain the corresponding figure of the fitting run [acf_figs_table.m](acf/acf_figs_table.m), this will call function  [plotsub_acf.m](acf/plotsub_acf.m) which will show the fit of each experiment in a subplot. 
Inline-style: 
![alt text](https://github.com/LauraPerezG/tweezers_AOP_tutorial/blob/merge_26nov_ales_lau/sec3_calibration/figures/ACF.jpg "Fit for Autocorrelation Function")

## Power Spectrum Density

Run [psd_analysis.m](psd/psd_analysis.m) to obtain the values of the stiffness and difussion coeficient given by the mean square displacement analysis method. This will call the functions, [psd_nlfit.m](psd/psd_nlfit.m) and [psd_lfit.m](psd/psd_lfit.m) contained in the folder  **[psd](psd/)**.


To obtain the corresponding figure of the fitting run [psd_figs_table.m](psd/psd_figs_table.m), this will call function  [plotsub_psd.m](psd/plotsub_psd.m) which will show the fit of each experiment in a subplot. 
![alt text](https://github.com/LauraPerezG/tweezers_AOP_tutorial/blob/merge_26nov_ales_lau/sec3_calibration/figures/PSD.jpg "Fit for Power Spectrum Density method")

## FORMA

Run [forma_analysis.m](forma/forma_analysis.m) to obtain the values of the stiffness and difussion coeficient given by FORMA.  This will call the functions, [forma1d.m](forma/forma1d.m) containded in the folder  **[forma](forma/)**.


## Bayesian
  
Run [bayesian_analysis.m](bayesian/bayesian_analysis.m) to obtain the values of the stiffness and difussion coeficient given by Bayesian analysis.  This will call the functions, [bayesian1d.m](bayesian/bayesian1d.m) containded in the folder  **[bayesian](bayesian/)**.


## WLSICE, non linear fit of correlated data

The functions [msd_nlfit.m](msd/msd_nlfit.m), [acf_nlfit.m](acf/acf_nlfit.m) and [acf_lfit.m](acf/acf_lfit.m) will call function [wlsice.m](statistics_func/wlsice.m) to perform a non-linear fitting considering the whole covariance matrix.  The target function that [wlsice.m](statistics_func/wlsice.m) uses depends on the fitting method. This is especified by the last input of  [wlsice.m](statistics_func/wlsice.m) function (*opt*).

opt=1 : ACF, linear fit

opt=2 : ACF, non-linear fit

opt=3 :  MSD, non-linear fit

 Follow [this link](https://www.nature.com/articles/s41598-018-24983-y). 
 to the the article by  Fogelmark, K., Lomholt, M.A., Irbäck, A. et al. Fitting a function to time-dependent ensemble averaged data. Sci Rep 8, 6984 (2018). About the wlsice method. 
[//]: # (| Method        | Input variables           | Output variables |
| ------------- |-------------| -------------|
|Drift method  | ... |  ...|
| Potential     |  x: Time series of the position <br> T:temperature  <br>  n_b: number of bins| x_bins: position of the histogram bins <br> mU: mean potential <br> EU: standar daviation of the potential <br> k_pot: stiffnes <br> Ek_pot:standar deviation of the stiffness <br> m_hist: mean histogram <br> E_hist: standard deviation of the histogram <br> h_0: ? <br> x_eq : equilibrium position |
| Equipotential   | x: Time series of the position <br> T:temperature  <br>  deltax: error in the position detection |  k_eq:stiffness <br> Ek_eq: standard deviation of the stiffness|
| Mean square displacement | x: Time series of the position <br> T:temperature  <br>  deltat: time between acquired frames <br> maxlag: lag integer |k_msd:stiffness <br> Ek_msd:standard deviation of the stiffness <br> D_msd: diffussion <br> ED_msd <br> tau:lag times <br> mmsd: mean mean squared displacement <br> Emsd: standard deviation of the mean square displacement <br> indc:? |
| Autocorrelation function | ... |  ...|
| PSD | ... |  ...|
|FORMA  | ... |  ...|)

***

This repository includes the analysis and plot codes for the the article **Optical Tweezers: A comprehensive Tutorial  from Calibration to Applications** by *[Jan Gieseler](https://scholar.google.com.ar/citations?user=6OKJlNgAAAAJ&hl=en), [Juan Ruben Gomez-Solano](https://www.fisica.unam.mx/es/personal.php?id=639), [Alessandro Magazù](http://softmatterlab.org/people/alessandro-magazzu/),  [Isaac Pérez Castillo](https://scholar.google.com.mx/citations?user=58GAc80AAAAJ&hl=en), [Laura Pérez García](http://softmatterlab.org/people/laura-perez-garcia/), [Martha Gironella-Torrent](https://scholar.google.com/citations?user=tITfJqkAAAAJ&hl=en), [Xavier Viader-Godoy](https://scholar.google.com/citations?user=dTLMJy0AAAAJ&hl=en), [Felix Ritort](http://ffn.ub.es/ritort/), [Giuseppe Pesce](https://scholar.google.com/citations?user=Sf4mmT8AAAAJ&hl=en), [Alejandro V. Arzola](https://orcid.org/0000-0002-4860-6330), [Karen Volke-Sepulveda](https://www.fisica.unam.mx/es/personal.php?id=27) and [Giovanni Volpe](http://softmatterlab.org/people/giovanni-volpe/)*.
