## Emergence of multifrequency activity in a laminar neural mass model

Welcome! 
The material you will find here provides the main tools to reproduce the results shown in in the paper: [Aristides et al. , *Emergence of multifrequency activity in a laminar neural mass model* (2025)].
If you use this material or any part of it in your work, please cite this paper to give proper credit.

### Codes: 
- A minimal Python code lanmm.py numerically solves the LaNMM model and calculates the power spectral density with np.fft.fft.
- A minimal Julia code lanmm_ly calculates the Lyapunov exponents of the LaNMM as a function of the external inputs to the pyramidal populations using Julia's ChaosTools.
- A minimal AUTO code bifur_lanmm.py returns the bifurcation diagram of the LaNMM as a function of the external input to the pyramidal population 1 and plots the results with Python's matplotlib. This script relies on the AUTO files c.lanmm and lanmm.f90.




