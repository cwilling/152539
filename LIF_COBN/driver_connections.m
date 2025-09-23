% driver code for CW's connections_COBN.c

%% 0.01 - mex the c code
  mex connections_COBN.c ran1.c gasdev.c % nb. "helper" codes

%% basic set up
%addpath('/Users/bap/bap_working/MatLabfiles/LIFmodels'); % for other codes
%addpath('/Users/bap/bap_working/MatLabfiles/LIFmodels/ChrisCodes'); % for model codes

% NETWORK PROPERTIES -----------------------------------------------------
  fprintf('\nCOBN sim: set neuron parameters \n')
% addpath('/bap_working/MatLabfiles/MatlabFiles/MarmosetBrain/Models/LIF')

% Ensure clean start
clear net_COBN;

% Number of excitatory and inhibitory neurons, eNnrn and iNnrn
  % # of neurons, totNnrn:  use 4:1 ratio
  net_COBN.eNnrn = 80; %800; %400; %40; %40; %400; % 240
  net_COBN.iNnrn = 20; %200;  %100; % 60

% Connectivity, i.e., connection probability, p:
  net_COBN.p = 0.2 ;  % Careful ! Not too large, if n small
 
  net_COBN.SEED_connections = 2; % positive integer number; const  for reproducible ran samples
  net_COBN.SEED_poisson = 3;     % positive integer number

% Generate A
  A = connections_COBN(net_COBN);

spy(A); title('A at end of driver_connections')
