## Matlab Scripts

* `bicontinuous_rho_r.m`: 

  Calculates the real space two point correlation functions of the bicontinuous structure (inset of *Fig.3*).

  * `bicontinuous_spectrum.m`, `bicontinuous_spectrum_time.m`:  

  Calculates the reciprocal space two point correlation functions of the bicontinuous structure.

  * `bicontinuous_spectrum_data.m` saves the calculated results

* `C_kt_example.m`: Draws an example of spatial-temporal correlation of the Gaussion form.

  `curvature*.m`: Calculates the curvature of the bicontinuous interface and. Calls the `curvatures.m` function.
  
* `isosurface_bicontinuous*.m`: Renders the 3D bicontinuous interface (*Fig.1d*).

* `isosurface_bicontinuous_2D*.m`: Plots the 2D bicontinuous interface (*Fig.1a-c*).

* `plot_curvature*.m`: Plots the curvature of the bicontinuous interface (*Fig.2*).

* `plot_heterogeneity.m`: Plots the 2D contour of the local orientation correlation.

* `plot_histogram.m`: Plots the probability distribution of the local orientation correlation.

* `plot_rho_r.m`: Plots the real space two point correlation functions of the bicontinuous structure (inset of *Fig.3*).

* `plot_temp_prop*.m`: Plots tau_M and the lengthscale of local orientation correlation distributions as a function of temperature (*Fig. 3*).

* `plot_temp_prop_Arial_color_inset.m`: Plots tau_M as a function of with colored marks.

* `susceptibility.m`: Calculates the four-point susceptibility.

### Matlab Functions

* `apply_ambient_occlusion.m`: Called by `isosurface_bicontinuous*.m`.
* `calculate_shape.m`: Calculates the cluster shape and orientations from the local gyration tensor (`rg_*.mat`).
* `curvatures.m`: Given the surface coordinates, calculates the local curvature.

### Image Outputs

* `BW`: The bicontinuous interface colored in black-and-white style.
