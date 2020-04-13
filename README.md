# layered-magnetic-structure-modeling

This project includes some approaches to model the magnetic structure for MnBi<sub>2</sub>Te<sub>4</sub>, as a layered magnetic material. Here we mainly discuss the **mean-field approach** (MF) applied to intralayer interactions, while the interlayer interactions are considered completely. Self-consistent methods are performed to get the thermodynamical value of the magnatization for each layer, and sometimes Monte-Carlo methods (Metropolis algorithm) are used for numerial integration. 

The setup of the system is fully described by the directions of the magnetic moment in each layer, and the class `Status` provides a direct description of a setup for the system. Since the absolute value of the magnetic moment is fixed, *2N* coordinates are enough to fully describe the whole setup. (*N* is the number of layers) The energy of one setup is computed, if the average values of the layer magnetism is obtained. So a self-consistant calculation process (iteration) must be applied to ensure that the thermodynamics expectation value of layer magnetism obtained from the energy with such average magnetization produces the same result as the average value itself used. Such a process is realized in `sys_run`. For high-dimension system, the thermodynamics expectation value is not easily obtained by normal integration methods, so a Monte Carlo simulation is needed to calculate the integral, realized in the class `M_C`.

The details of mathematical derivations will come soon, in the form of a refrence to the paper where this program is used.
