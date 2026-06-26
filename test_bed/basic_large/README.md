# Basic large test bed

## Setup and design choices

This is a test bed for a large section of the Rosario Strait. It is basic because it uses the channel builder for the bathymetry and velocity setup as opposed to FVCOM or field data. It is based on snapshots of FVCOM data, including the velocity as a function of depth and the bathymetry, as seen below:

![Velocity as a function of depth](figures/FVCOM_image_velocity_cross_sections.png) 

![Bathymetry overlaid with preliminary segments for channel builder](figures/FVCOM_image_initial_segments.png)

The bathymetry above is overlaid with preliminary segments for the channel builder setup. These 8 segments were modified so that they could overlap, and two more segments (not shown) were added to better approximate portions of the domain:

![Modified segments](figures/FVCOM_image_modified_segments.png)

Using these segments, the computational domain looks like this (on a coarse mesh), reasonably mimicking the actual bathymetry from the FVCOM data.

![Coarse top view](figures/top_view_coarse_bathymetry.png)

Due to the large lateral extent of the domain, the overall domain has a high aspect ratio (looks like a pancake), as seen below in the 3D image and 2D slice:

![3D image at an angle](figures/3D_pancake.png)

![2D slice in the middle](figures/2D_pancake.png)

## Computational observations

With no additional refinements (using amr.max_level = 0), the projections converged easily without much parameter modification. Adding mesh refinements proved to be more difficult for the solver. Using field refinement alone to add cells at the terrain boundary, both amr.max_level = 1 and amr.max_level = 2 simulations could be brought to converge, but the nodal projection would stall out (use the maximum number of iterations) later on. Introducing a uniform mesh resolution across the interface fixed this, but this does increase the total number of cells.

These modifications helped the amr.max_level = 1 case run, but the amr.max_level = 2 case blew up immediately. Will continue to investigate. The image below shows some of the mesh from the successful run, indicating the nature of the mesh refinements:

![Mesh slices with max_level = 1](figures/mesh_slices.png)

## Flow observations

From the longest run performed, which shows the velocity at hub height (z = -17 m), the flow appears to develop well, but the free surface motion is significant and impactful at this depth, which is a concern. Because it is not in our interests to fully resolve the surface waves, we will need to add some functionality to handle this spurious phenomenon, which will likely be a type of relaxation zone for the free surface.

![Hub height velocity video, amr.max_level = 1](https://github.com/user-attachments/assets/387e5052-a72e-4d83-a3e7-73131a9d50be)
