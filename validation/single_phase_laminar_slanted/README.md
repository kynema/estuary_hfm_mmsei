# Slanted Poiseuille Flow

This is a laminar channel flow in a single phase with a channel that has been slanted at an angle $\theta$ as measured in the $x$-direction.  The flow is represented as a half-channel using a symmetry plane.  The analytical solution for velocity is applied as Dirichlet inflow.  The outflow is Dirichlet on pressure, and a pressure gradient naturally forms.

Each subdirectory in the `cases` directory contains an input file with two lines:

~~~
FILE = ../base-slanted-poiseuille.inp
amr.n_cell = Nx Ny Nz 
~~~

These directories are for storing `plt` and `chk` files at different grid resolutions for error convergence analysis. The input file `base-slanted-poiseuille.inp` is shared by all cases and can be updated using values outputed from the pre-processing script `python/SlantedChannelConfig.py`.

## Pre and Post-Processing
There are two associated python pre- and post-processing scripts. These scripts require the following python packages:

~~~
numpy matplotlib yt
~~~ 

Pre-processing is handled by the `SlantedChannelConfig.py` script, which prints out input parameters for the simulation provided an angle $\theta$ in degrees or radians. Use the "help" option to see a full list of optional parameters and example usage:
~~~
python python/SlantedChannelConfig.py -h
~~~

Post processing is handled by running each cell of the `SlantedChannelPostProcess.ipynb` notebook. Note that this will only look for data in the `cases/slanted-ibfm-###` directories, and by default grabs the last `plt` file. 

## Analytical Solution

Provided a 2D domain in the x-z plane with periodic boundaries in $y$,the analytical solution for fully-developed parabolic (Poiseuille) channel flow in a tilted, vertically-shifted configuration is:

$$u(r) = U_{\max} \left(1 - \left(\frac{2r}{H}\right)^2\right)$$

where:
- $r$ = perpendicular distance from centerline
- $U_{\max}$ = maximum velocity (at centerline, $r=0$)
- $H$ = pipe diameter
- $u(r) = 0$ for $r > H/2$ (at pipe wall)

Here, $r$ is measured in the x-z plane perpendicular to the tilted centerline as:

$$r = |x\sin\theta - (z-s)\cos\theta|.$$

The centerline is shifted vertically by $s$ in the $z$-direction, and rotated at an angle $\theta$ about the point $(0, 0, s)$ in x-z plane:

$$\mathbf{C}(t) = [t\cos\theta, \quad 0, \quad s + t\sin\theta], \quad t > 0$$

Flow is along the tilted centerline in the x-z plane:

$$\mathbf{u}(x,z) = u(r) \cdot (\cos\theta, \sin\theta)$$

$$u_x = u(r)\cos\theta, \quad u_z = u(r)\sin\theta$$

