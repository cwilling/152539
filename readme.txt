This is the readme for the two neural network models (mex file) associated with the following paper:

S.Cavallari, S. Panzeri and A.Mazzoni (2014) Comparison of the
dynamics of neural interactions between current-based and conductance-based 
integrate-and-fire recurrent networks, Frontiers in Neural Circuits
8:12. doi: 10.3389/fncir.2014.00012

Recurrent networks, each of two populations (excitatory and inhibitory) of
randomly connected Leaky Integrate-and-Fire (LIF) neurons with either
conductance-based synapses (COBN) or current-based synapses (CUBN)
were studied. The activity of the LIF COBN model were compared with the
activity of the associated model LIF CUBN. 
Instructions are provided in the ReadMe files in each model associated sub folder, LIF_COBN and
LIF_CUBN.

If you have any questions about the implementation of these matlab
models, which require compilation with mex, contact:
ste.cavallari@gmail.com

Please cite the paper if you use the codes.


20140709 Comments in LIF_COBN/code_COBN.c, LIF_COBN/code_COBN.m,
LIF_CUBN/code_CUBN.c, LIF_CUBN/code_CUBN.m were enhanced.

20150722 The code to generate the OU process has been added (LIF_COBN/OU_process.m
and LIF_CUBN/OU_process.m) together with the instructions
(LIF_COBN/ReadMe_COBN and LIF_CUBN/ReadMe_CUBN)
to set the arguments to generate the data used in some figures of the paper.
