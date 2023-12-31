# high-z-hubble-diagram

This is a Python implementation of using quasars as standard candles described in ["Quasars as standard candles III. Validation of a new sample for cosmological studies"](https://doi.org/10.1051/0004-6361/202038899), A&A 2020.
The catalog used in this script is available at [http://cdsarc.u-strasbg.fr/viz-bin/cat/J/A+A/642/A150](http://cdsarc.u-strasbg.fr/viz-bin/cat/J/A+A/642/A150).  

The steps of this code are as follows:  
1. Perform a linear fit on the UV and x-ray flux data to find the log-linear relationship.
2. Use fitted parameters to calculate luminosity distance $d_L$.
3. From $d_L$ calculate distance modulus DM and e_DM.
4. Plot the Hubble diagram of DM vs. z and fit the diagram with MCMC.
