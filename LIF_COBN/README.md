#### Overview.

The original COBN code is extended here to allow multiple, successive, calls from Matlab driver code.

This facilitates a sequence of short time simulations of one LIF group, within a longer time sequence.
Thus two, or multiple, interacting LIF groups can be simulated.

The original code has been divided into two pieces (_connections\_COBN.c_ & _processing\_COBN.c_) each with its own driver code
(_driver\_connections.m_ & _driver\_processing.m_). The first piece of code generates a connectivity matrix which is used,
along with various network parameters, by the second piece to solve the DE's in sequential time windows.
The extended code preserves a single connectivity matrix, along with multiple internal variables which capture the current LIF state.
These are recycled between the Matlab driver and C codes for successive time windows.

The new code can be used to simulate interacting groups of LIF neurons. Simulation time windows down to 16 time steps work OK
and reproduce the results of longer simulations (> 1sec) quite closely.

#### Installation.
Either clone this repository or unpack a release tarball somewhere e.g.
```
  mkdir ~/src
  cd ~/src
  git clone https://github.com/cwilling/152539
```
Then run Matlab, changing directory to the newly created ~/src/152539/LIF_COBN directory.

#### Usage.

If installed as above, under Matlab run _`driver_processing`_

Basic characteristics of the neuron's connections matrix
have been set in _driver\_connections.m_ which runs the C code in connections_COBN.c to actually generate the connections
matrix. The _driver\_processing.m_ file sets various network parameters (in the net_COBN structure), in particular _sample_width_
which determines the number of time steps to be processed in each iteration of the C code solver ( _processing\_COBN.c_).
Change the parameters in _driver\_connections.m_ and _driver\_processing.m_ to test alternate scenarios.

#### Authors (2025 extensions).
Project Lead: Prof. Bernard Pailthorpe, Physics, University of Sydney, Australia
e: bernard.pailthorpe@sydney.edu.au

Coding: Christoph Willing
[https://github.com/cwilling/152539/issues](https://github.com/cwilling/152539/issues)

