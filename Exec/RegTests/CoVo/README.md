Convected Vortex (CoVo)
================================

Description
-----------

The CoVo case is a classical test case for computational fluid dynamics solvers. It allows to directly evaluate the numerical characteristics of the solver by comparison with analytical solution. It consists in convecting an isentropic vortex across a periodic 2D box. The analytic solution of the isentropic vortex is given by the stream function:

<a href="https://www.codecogs.com/eqnedit.php?latex=\large&space;\psi(x,y)&space;=&space;\psi_0&space;e^{&space;-&space;{(x-x_c)^2&plus;(y-y_c)^2&space;\over&space;2&space;r_c^2}&space;}=\psi_0&space;e^{-r^2/2r_c^2}" target="_blank"><img src="https://latex.codecogs.com/svg.latex?\large&space;\psi(x,y)&space;=&space;\psi_0&space;e^{&space;-&space;{(x-x_c)^2&plus;(y-y_c)^2&space;\over&space;2&space;r_c^2}&space;}=\psi_0&space;e^{-r^2/2r_c^2}" title="\large \psi(x,y) = \psi_0 e^{ - {(x-x_c)^2+(y-y_c)^2 \over 2 r_c^2} }=\psi_0 e^{-r^2/2r_c^2}" /></a>


User inputs
-----------

In the present implementation, the vortex position, size and strength are controlled in the probin file using the following keywords:
 * Position: `xvort`, `yvort` [m]
 * Characteristic size: `rvort` [m]
 * Strength: `forcevort` []

Additionally, the underlying uniform velocity magnitude and direction can be adjusted from the probin file using:
 * Velocity magnitude: `meanFlowMag` [m/s]
 * Streamwise direction: `meanFlowDir` 1 -> X, 2 -> Y, 3 -> 45 angle between X/Y 


PeleLM testing
--------------

The CoVo is used as part of PeleLM testing suite:
 * It is used along the in Pele-convergence.ini along with the multirun.py and ppConvOrder.py scripts to evaluate the convergence order of the PeleLM convection scheme.
 * It is used as part of the regression testing suite (Pele-regtest.ini), changing the viscosity from 0.0, to a constant value of 2.00e-5 and finally to a local value computed using EGLiB in order to test gradually various sections of the code.
