## Overview

The original COBN code is extended here to allow multiple, successive, calls from Matlab driver code.

This facilitates a sequence of short time simulations of one LIF group, within a longer time sequence.
Thus two, or multiple, interacting LIF groups can be simulated.

The original code has been divided into two pieces (_connections\_COBN.c_ & _process\_COBN.c_) each with its own driver code
(_driver\_connections.m_ & _driver\_process.m_). The first piece of code generates a connectivity matrix which is used,
along with various network parameters, by the second piece to solve the DE's in sequential time windows.
The extended code preserves a single connectivity matrix, along with multiple internal variables which capture the current LIF state.
These are recycled between the Matlab driver and C codes for successive time windows.

The new code can be used to simulate interacting groups of LIF neurons. Simulation time windows down to 16 time steps work OK
and reproduce the results of longer simulations (> 1sec) quite closely.


#### Project Lead
Prof. Bernard Pailthorpe, Physics, University of Sydney, Australia
e: bernard.pailthorpe@sydney.edu.au

#### Coding
Christoph Willing
