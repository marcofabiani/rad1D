# rad1D

rad1D is a 1 dimensional radiation model for liquid rocket engines. 

The model computes the radiative heat flux along the thrust chamber wall as
$$q = \sigma \frac{1}{\left(\frac{1}{\varepsilon_\mathrm{gas}+\frac{1}{\varepsilon_\mathrm{wal}}-1}\right)}\left(T_\mathrm{gas}^4-T_\mathrm{wall}^4\right) $$

For each nozzle abscissa the pressure, temperature and composition are computed using rocketcea. The gas emissivity is computed using a custom weighted-sum-of-gray-gases model developed specifically for liquid rocket engines [1,2]. 


Requires the rocketcea package and numpy


## References
[1] Marco Fabiani, Mario Tindaro Migliorino, Daniele Bianchi, and Francesco Nasuti, 
"Spectral and Global Radiative Heat Transfer Models for Liquid Propellant Rocket Engines, " 
JPP, Vol. 41, No. 5 (2025), pp. 650-664 doi: doi/abs/10.2514/1.B39892

[2] Marco Fabiani, Mario Tindaro Migliorino, Daniele Bianchi, and Francesco Nasuti, 
"Reduced-order Models for Radiative Heat Loads Estimation in Liquid Rocket Engines,"
IAC 2025, Sydney

