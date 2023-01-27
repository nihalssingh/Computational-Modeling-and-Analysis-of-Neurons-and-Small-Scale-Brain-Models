# Computational-Modeling-and-Analysis-of-Neurons-and-Small-Scale-Brain-Models

Taking insights from cutting-edge neuroscience
research, this work proposes to develop a generic brain
simulation interface. 
In this work, we have developed a scalable
interface that can be used to model different regions of the brain
by specifying the distribution and dynamics of different kinds of
neurons specific to the region of interest. 
We have used MATLAB
to implement and test our simulator. We have used raster plots to
benchmark our circuits. 
As an illustration, we have also modeled
the third layer of a cat's visual cortex and obtained the raster
plots to study its spiking sequence.

We have uploaded 3 files excluding this README file.

Brain_Simulator_generic.m - This is a generic Brain Simulator where we have 
randomized synaptic connections and weights for a random mix of neurons in 
randomized percentages. (This code is best to view the impact of synaptic 
connections, and see the general  behaviour of different neurons across a 
transient simulation.

CatCortex_Sim.m - While the previous code is all randomized for no particular 
case, this code is modelled based on the neuron models and corresponding 
percentages found in Level 3 of a Cat's Visual Cortex. Running this code will 
generate the spike raster plot discussed in the report/slides.

User_input_BSim.m - This code allows users to customize exactly which 
circuits/layers to simulate. The user-interface allows one to input the number
of neurons, percentages of each neuron type, simulation time step, transient 
time, number of synapses for noise modelling, and specific explicit synaptical 
connections.
