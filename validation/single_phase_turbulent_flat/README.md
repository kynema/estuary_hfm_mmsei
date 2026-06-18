# Single Phase Turbulence in Flat Channel

This validation case simulates turbulent channel flow at three stress Reynold's numbers $\text{Re}_\tau = $ 180, 395, and 934 corresponding to the cases described in [Kim, Moin & Moser (1986)](https://doi.org/10.1017/S0022112087000892), [Moser, Kim & Mansour (1999)](https://doi.org/10.1063/1.869966), and [Hoyas & Jimenez (2006)](https://doi.org/10.1063/1.2162185), respectively. A first DNS is performed at $\text{Re}_\tau = $ 180, then the LES models implementation is tested at higher Reynolds numbers. 

For all cases, the configurartion is periodic in the $x$ and $y$ directions and wall boundaries are imposed in the $z$ direction using either `ChannelBuilder` for IBFM, or AMReX's boundary conditions. _Note: this differs from the original sources where periodicty was in the $x$ and $z$ directions and wall boundaries were imposed in the $y$ direction. Solutions are compared using the original source coordinate system, which requires mapping `kynema-sgf` solutions from $v_\text{sgf} \rightarrow w_\text{data}$ and $w_\text{sgf} \rightarrow v_\text{sol}$, as shown below._ 

![Coordinate System](figures/CoordinateSystem.png)

## DNS results at Re = 180
The channel half width is set to $\delta = 0.005$ m. For simulations without IB, the upper and lower walls in the $z$ direction are set at $\pm \delta$, respectively. For simulations with IB, the computational domain is extended $\pm 0.0052$ m in $z$ with IB intersecting the domain at $\pm \delta$. A background pressure gradient is imposed in the $x$ direction to compensate for wall friction. The base grid and two levels of refinement is described in the table below:


|           | $x$   | $y$   | $z$                |
|---------- | ----- | ----- | ------------------ |
|Domain Size| 6.24 $\delta$ | 3.12 $\delta$ | 2.0 $\delta$ or 2.08 $\delta$ (IB) |
| Level 0   | 384   | 192   | ___ or 128 (IB) |
| Level 1   | 768   | 384   | ___ or 256 (IB) |
| Level 2   | 1536  | 768   | ___ or 512 (IB) |

The fluid in the simulation is air at ambient pressure and a temperature of 750 K. The physical properties used in the simulation are provided in the table below:

| Pressure (Pa) | Temperature (K) | Density  | $\mu$      | $\nu$      |
| ------------- | --------------- | -------- | ---------- | ---------- |
| 101325.0      | 750.0           | 0.468793 | 3.57816e-5 | 7.63271e-5 |

The characteristics of the flows are reported below:

| $\text{Re}_\tau$ | $u_\tau$ | $\tau_w$ | $dp/dx$ | $t^* = \delta/u_\tau$ |
| ---------------- | -------- | -------- | ------- | --------------------- |
| 180.2            | 2.75122  | 3.5485   | -709.79 | 1.817e-3              |

Two levels of refinement, targeted on the EB, are employed in order to sufficiently resolve the boundary layer. The mesh characteristics are summarized in the table below. The $y^+$ value is that of the cell center of the first full cell (uncut by the IB).

**Mesh characteristics**

| $\text{Re}_\tau$ | $\Delta z^+$ (L0) | $\Delta z^+$ (L2) | $z^+$ | Cells count |
| ---- | ---- | ---- | ---- | ---- |
| 180.2 |  | |  |  M |

Simulations are carried out for 20 eddy turn over time $t^*$ to reach statistically steady conditions and data are then spatially averaged in the periodic directions and averaged in time over 10 $t^*$ to get the velocity statistics in the direction normal to the wall.
