# pdeoptim

These are the MATLAB files accompanying the manuscript  

```
   A spectral-in-time Newton--Krylov method 
   for nonlinear PDE-constrained optimization
   by Stefan Guettel and John W. Pearson, 2020.
```

There are five driver files that generate all the figures and tables in that paper:

* **example_Schloegl_FD_eig.m** - Figure 3.1
* **example_Schloegl_spectral_beta.m** - Figure 3.2 and Table 3.1
* **example_Schloegl_FD_n.m** - Figure 3.3 and Table 3.2
* **example_Schloegl_FD_32_err.m** - Figure 3.4
* **example_Reactdiff_FD_n.m** - Figure 3.5 and Table 3.3

The produced tables will not contain any timings as, to obtain some 
statistical reliability, these are produced from ten consecutive runs, 
which makes the code run rather slow. If you want to produce timings on 
your computer, change the `runs` variable at the beginning of the scripts
from `0` to `10`. This then corresponds exactly to the setting used to 
produce the timings in the paper. The timings in the paper were obtained 
in MATLAB 2019A on a Windows 10 laptop with 8 GB RAM and an Intel(R) 
Core(TM) i7-8650U CPU running at 2 GHz and will vary on different architectures.

The solver code requires Chebfun in the MATLAB path. When first running
any of the example files, the code will try to download and install
Chebfun automatically. If this fails for some reason, please download 
Chebfun v4.2.2889 from https://www.chebfun.org/download/chebfun_v4.2.2889.zip 
and add the folder 'chebfun' in this directory with the other files. This 
particular Chebfun version was also used in our codes for a previous paper 
on linear optimal control problems available at http://guettel.com/dccontrol. 
We use the same for reproducibility and to facilitate comparisons between 
both solvers.

