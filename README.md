# Computing Flow Around Obstacles - Uisng Lattice Boltzmann Methods

Lattice Boltzmann methods (LBM), originated from the lattice gas 
automata (LGA) method (Hardy-Pomeau-Pazzis and Frisch-Hasslacher-Pomeau 
models), is a class of computational fluid dynamics (CFD) methods for
fluid simulation. Instead of solving the Navier–Stokes equations directly,
a fluid density on a lattice is simulated with streaming and collision 
(relaxation) processes. The method is versatile as the model fluid 
can straightforwardly be made to mimic common fluid behaviour like 
vapour/liquid coexistence, and so fluid systems such as liquid droplets 
can be simulated. 

More Information: https://en.wikipedia.org/wiki/Lattice_Boltzmann_methods

Based on lectures: 
Simulation and modeling of natural processes
Week 5: Lattice Boltzmann modeling of fluid flow
Authoor: Jonas Latt
Universite de Geneve, Switzerland

Course link: 
https://www.coursera.org/learn/modeling-simulation-natural-processes

Modified by: Harivinay V
Github: https://www.github.com/M87K452b

--- Code restructured only a bit. Contains additional comments
--- Added functionality to export a video directly from code

Dependencies: 
    ffmpeg - install using pip install ffmpeg-python
    graphviz - install using pip install graphviz

*  Inputs can me set up as desired in the flow conditons sections
*  The obstacle geometry can be changed by entering the equation of a desired shaped in obstacle_func

## Demo
[![Flow around an ellipse](https://img.youtube.com/vi/HxTHIJGJsYY/0.jpg)](https://youtu.be/HxTHIJGJsYY)

**Click on the above image to see the simulation.(redirected to youtube)**

