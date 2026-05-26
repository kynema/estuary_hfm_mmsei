# Slanted Poiseuille Flow

This is a laminar channel flow in a single phase with a channel that has been slanted at an angle $\theta$ as measured in the $x$-direction.  The flow is represented as a half-channel using a symmetry plane.  The analytical solution for velocity is applied as Dirichlet inflow.  The outflow is Dirichlet on pressure, and a pressure gradient naturally forms.

## Pre and Post-Processing
There are associated python configuration and post-processing scripts. The configuration script prints out input parameters for the simulation provided the angle $\theta$ in degrees of radians. To see some examples:
~~~
cd python
python SlantedChannelConfig.py -h
~~~

## Analytical Solution

Provided a 2D domain in the $x$-$z$ plane with periodic boundaries in $y$,the analytical solution for fully-developed parabolic (Poiseuille) channel flow in a tilted, vertically-shifted configuration is:

$$u(r) = U_{\max} \left(1 - \left(\frac{2r}{H}\right)^2\right)$$

where:
- $r$ = perpendicular distance from centerline
- $U_{\max}$ = maximum velocity (at centerline, $r=0$)
- $H$ = pipe diameter
- $u(r) = 0$ for $r > H/2$ (at pipe wall)

Here, $r$ is measured in the $x$-$z$ plane perpendicular to the tilted centerline as:

$$r = |(z-s)\cos\theta - x\sin\theta|.$$

The centerline is shifted vertically by $s$ in the $z$-direction, and rotated at an angle $\theta$ about the point $(0, 0, s)$ in $x$-$z$ plane:
$$\mathbf{C}(t) = (t\cos\theta, 0, s + t\sin\theta)$$

Flow is along the tilted centerline in the x-z plane:

$$\mathbf{u}(x,z) = u(r) \cdot (\cos\theta, \sin\theta)$$

$$u_x = u(r)\cos\theta, \quad u_z = u(r)\sin\theta$$

