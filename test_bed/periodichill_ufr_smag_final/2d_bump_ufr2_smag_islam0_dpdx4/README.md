### UFR 3-30: Flow Over Periodic Hills 

This is an infinitely wide and long "periodic hill" test case, based on the geometry defined by Mellen et al. (2000). A row of 2D hills (height h=0.028m) sits in a channel with total height 3.036h, hill spacing (period) 9h, hill length 3.857h, defined by six polynomial segments. The flow separates at the hill crest, forms a recirculation zone in the lee of the hill, and reattaches naturally around x/h ≈ 4.5 before accelerating over the next hill. This is a canonical case for testing separation from a curved surface.

The case as implemented uses the following parameters:

| Parameter | Value |
| --- | --- |
| h | 0.028 m |
| $Re_h$ | 10595 |
| $dp/\rho dx$ | 0.067 |
| $U_{bulk}$ | 0.3785 m/s |

In post-processing, $U_{bulk}$ is estimated as $0.85\bar{U}_{max}$. Though not exact, $dp/\rho dx$ was modified until an approximate $U_{bulk}$ and therefore $Re_h$ was achieved. 

One flow-through time is approximately 0.6658 seconds. Results in the output file labeled with "h_mean" are averaged both on flow-through time and the channel span direction. Results in the output file labeled with "center_slice" are only averaged on flow-through time.