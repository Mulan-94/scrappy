### Introduction
This is not fully-fledged Python script. It is a simple shim availed to create a depolarization ratio map.

The depolarization ratio is the ratio of fractional polarisation at a shorter wavelength, to that at a longer wavelength. The fractional polarisation, $P$,  is defined using Stokes $\mathit{I, Q, U}$  terms as:

$$
\begin{equation*}
P = \frac{\sqrt{Q^{2} + U^{2}}}{\mathit{I}}
\end{equation*}
$$

Given fractional polarization at different frequencies, $\nu_{max}$ and $\nu_{min}$, the depolarization ratio is given as:

$$
\begin{equation*}
\text{Depolarisation ratio}=\frac{P\text{ at }\nu_{max}}{P\text{ at }\nu_{min}}
\end{equation*}
$$


<!-- ------
### Example -->


<!-- ------
### Outputs -->

------

### Docstring Documentation

::: src.rmsynthesis.rmmap_x2.depolarisation_ratio
    handler: python