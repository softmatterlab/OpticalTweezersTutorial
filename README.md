# tweezers_AOP_tutorial
We provide the set of scripts and functions described in the AOP tutorial.



| Method        | Input variables           | Output variables |
| ------------- |-------------| -------------|
|Drift method  | ... |  ...|
| Potential     |  x: Time series of the position <br> T:temperature  <br>  n_b: number of bins| x_bins: position of the histogram bins <br> mU: mean potential <br> EU: standar daviation of the potential <br> k_pot: stiffnes <br> Ek_pot:standar deviation of the stiffness <br> m_hist: mean histogram <br> E_hist: standard deviation of the histogram <br> h_0: ? <br> x_eq : equilibrium position |
| Equipotential   | x: Time series of the position <br> T:temperature  <br>  deltax: error in the position detection |  k_eq:stiffness <br> Ek_eq: standard deviation of the stiffness|
| Mean square displacement | x: Time series of the position <br> T:temperature  <br>  deltat: time between acquired frames <br> maxlag: lag integer |k_msd:stiffness <br> Ek_msd:standard deviation of the stiffness <br> D_msd: diffussion <br> ED_msd <br> tau:lag times <br> mmsd: mean mean squared displacement <br> Emsd: standard deviation of the mean square displacement <br> indc:? |
| Autocorrelation function | ... |  ...|
| PSD | ... |  ...|
|FORMA  | ... |  ...|

