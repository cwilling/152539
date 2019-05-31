Recurrent random network with excitatory and inhibitory Leaky Integrate-and-Fire neurons with current-based 
synapses from the paper “Comparison of the dynamics of neural interactions between current-based and 
conductance-based integrate-and-fire recurrent networks” written by S.Cavallari, S.Panzeri and A.Mazzoni and 
published in Frontiers in Neural Circuits (2014), 8:12. doi:10.3389/fncir.2014.00012. The paper compares the 
activity of this current-based network (i.e. code_CUBN.c) with the activity of a comparable network of LIF 
neurons with conductance-based synapses (whose source code is in the “LIF_COBN” folder). 

The function code_CUBN.c is a mex source code. You have to compile this routine in Matlab to generate the 
mex file (e.g.: code_CUBN.mexw64). Note that you have to include the functions ran1.c and gasdev.c in the 
compiling instruction in the Matlab workspace, in the following way: 
mex code_CUBN.c ran1.c gasdev.c 

After you compiled the function, you can call it as specified in the help.
For more information use the help of the function (and see the example below): 
help code_CUBN 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
In the following an easy example for setting the arguments to generate the data used in figures 4I; 7 (i.e. LFP) 
and 6A (i.e. average firing rate) once you compiled the mex file.

In the Matlab workspace:

1.	parameters_CUBN; % to generate the structure net_CUBN with all the parameters of the network
2.	simulation_length = 4500; % units: (ms)
3.	M = simulation_length/net_CUBN.Dt; % length of the simulation in time steps
4.	external_signal_intensity = 2; % units; (spikes/ms)/cell 
5.	external_signal = ones(M,1) * external_signal_intensity * net_CUBN.Dt;
6.	SEED_OU = 1; % positive integer number
7.	external_noise = OU_process(M, net_CUBN.Dt, 16, 0.16*net_CUBN.Dt, SEED_OU);
8.	INPUT2E = external_signal + external_noise;
9.	INPUT2I = INPUT2E;
10.	SEED_connections = 2; % positive integer number
11.	SEED_poisson = 3; % positive integer number
12.	[E2EI,I2EI,eFR,iFR] = code_CUBN(net_CUBN, INPUT2E, INPUT2I, SEED_connections, SEED_poisson);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Question on how to use the model should be addressed to:
ste.cavallari@gmail.com

Please cite the paper if you use the code.
