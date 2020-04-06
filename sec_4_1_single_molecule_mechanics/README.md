# Section 4.1: Single molecule mechanics
 This folder contains the codes to analize the expeirments of  a single molecule pulling experiment discussed in the text.
 
How to use the Single molecule mechanics section:

**Free energy of formation of a DNA hairpin**

Run [work_distributions](program/work_distributions_Fig16.m) in folder **[programs](programs/)** to obtain  (A) Several pulling trajectories of the hairpin. In red are plotted the unfolding cycles and in blue the folding ones. The black dots indicate the initial ( L<sub>0</sub> ; f<sub>min</sub>) and final L<sub>1</sub> 
 ; f<sub>max</sub>) points considering for computing the work, WFU. (B) Probability density distributions of the work, computed as in equation (6). The point at which the probabilities cross is the equilibrium difference of free energy of the hairpin. The red and blue lines shown are the interpolations of the bars.

The files contained in **[data](data/)** '' folder  u,f for the unfolding(stretching) and
folding(releasing) trajectories. They contain 4 columns:
$1=counter $2=trap distance(nm) $3=time (s) $4=force(pN)
The files have been sorted in order to have the trap distance order toward
increasing values.

The program will run the correspondent analysis and  save the output in folder **[results](results/)**, this is useful if you want to use another program for doing th plot. 


![alt text](figures/work_distributions.jpg)


***


 
This repository includes the analysis and plot codes for the the article **Optical Tweezers: A comprehensive Tutorial  from Calibration to Applications** by *[Jan Gieseler](https://scholar.google.com.ar/citations?user=6OKJlNgAAAAJ&hl=en), [Juan Ruben Gomez-Solano](https://www.fisica.unam.mx/es/personal.php?id=639), [Alessandro Magazù](http://softmatterlab.org/people/alessandro-magazzu/),  [Isaac Pérez Castillo](https://scholar.google.com.mx/citations?user=58GAc80AAAAJ&hl=en), [Laura Pérez García](http://softmatterlab.org/people/laura-perez-garcia/), [Martha Gironella-Torrent](https://scholar.google.com/citations?user=tITfJqkAAAAJ&hl=en), [Xavier Viader-Godoy](https://scholar.google.com/citations?user=dTLMJy0AAAAJ&hl=en), [Felix Ritort](http://ffn.ub.es/ritort/), [Giuseppe Pesce](https://scholar.google.com/citations?user=Sf4mmT8AAAAJ&hl=en), [Alejandro V. Arzola](https://orcid.org/0000-0002-4860-6330), [Karen Volke-Sepulveda](https://www.fisica.unam.mx/es/personal.php?id=27) and [Giovanni Volpe](http://softmatterlab.org/people/giovanni-volpe/)*.
