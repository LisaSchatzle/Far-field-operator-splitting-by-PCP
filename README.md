# Far-field-operator-splitting-by-PCP

## License

Copyright (c) 2025, Roland Griesmaier, Lisa Schätzle

This software is released under GNU General Public License, Version 3.
The full license text is provided in the file LICENSE included in this repository
or can be obtained from http://www.gnu.org/licenses/

## Content of this project

This is a guide to generate the figures and tables that have been used in the work

**"Far field operator splitting by principal component pursuit"**

by Roland Griesmaier and Lisa Schätzle.

The following versions of the paper are available:

- [ ] 

You find all needed Matlab files to generate the figures and tables.

## Requirements

The following additional software is necessary to run the code in this project:

- [ ] a recent version of Matlab

## An overview

- [ ] *addNoise.m* adds p% complex valued uniformly distributed additive error to a matrix.
- [ ] *elementwise_soft_thresholding.m* applies elementwise soft thresholding to a complex matrix.
- [ ] *evaluateFarfieldNystrom.m* evaluates the far field patterns for n2 incident and observation directions on an equidistant grid on the unit sphere for each configuration, i.e. simulates the far field operator, by using a Nyström method.
- [ ] *example_4_1.m* tests RPCP reconstruction for two scatterers and synthetic data, study on optimal choice of coupling parameter depending on number of illumination and observation directions.
- [ ] *example_4_2.m* tests RPCP reconstruction for two scatterers and synthetic data depending on amount of equally distributed added random noise to the data.
- [ ] *example_4_3_1.m* tests RPCP reconstruction for two scatterers and synthetic data depending on distance between scatterer's coponents.
- [ ] *example_4_3_2.m* tests RPCP reconstruction for two scatterers and synthetic data depending on sizes of scatterer's coponents.
- [ ] *figure_2_1.m* provides plots of the absolute values of the expansion coefficients of the two components in a sparse plus low-rank far field operator split as well as of the geometry of the corresponding scatterers.
- [ ] *figure_2_2.m* provides plots of the incoherence measure beta as well as of its upper bound depending on t=|c_1-c_2| and for two different fixed values of N_1.
- [ ] *kurve.m* provides curve values for different shapes of scatterers.
- [ ] *MyCMap.mat* provides colormap for the plots of expansion coefficients of far field operators.
- [ ] *RPCP.m* solves the relaxed principal component pursuit problem for far field operator splitting  
by using an accelerated proximal gradient method.
- [ ] *translOp.m* applies generalized translation operator to a far field matrix.

## Generating the files

For generating the Figures 2.1 and 2.2 from the work, run

* *figure_2_1.m* -  **Figure 2.1** (Example 2.1, geometry of kite-shaped and nut-shaped scatterer, expansion coefficients of components in sparse plus low-rank far field operator split)
* *figure_2_2.m* -  **Figure 2.2** (Remark 2.5, plots of the incoherence measure beta as well as of its upper bound depending on t=|c_1-c_2| and for two different fixed values of N_1)

For generating the Figures 4.1--4.5, run

* *example_4_1.m* - **Figure 4.1** (Example 4.1, splitting for varying number of illumination and observation directions and different selection strategies of coupling parameter)
* *example_4_2.m* - **Figure 4.2** (Example 4.2, splitting for varying noise level)
* *example_4_3_1.m* - **Figure 4.3** (Example 4.3, splitting for varying distance between scatterer's components)
* *example_4_3_2.m* - **Figure 4.4**, **Figure 4.5** (Example 4.3, splitting for varying sizes of scatterer's components)
