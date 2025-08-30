

function loop_COBN (sample_width)
  arguments
	sample_width = 16;	% time steps of input window
  end
  restart = 0;

  fprintf("\nloop_COBN() Iteration = 1\n");

% MATLAB code to set up call to code_COBN.c      (nb. % for comments)

net_COBN.restart = restart;
net_COBN.sample_width = sample_width;

%% 0.1 cf. Sim codes, time vec
  x0=0; xf= 20.0; % t range (sec)
  h = 0.0001;   % t step (1 ms; [& 0.5, 0.25 ms, for debug])
  x = [x0:h:xf];  % t domain array
  x=x(2:end); % match vector length, if needed
  simL_sec =max(x)
  clear simL* 

%% 1.0 set up THE NETWORK struct (parameters)  { %% indicates a nee code section }
% Time resolution [ms]
  net_COBN.Dt = 0.1; %0.05;
 
% NETWORK PROPERTIES -----------------------------------------------------
  fprintf('\n COBN sim: set neuron parameters \n')
% addpath('/bap_working/MatLabfiles/MatlabFiles/MarmosetBrain/Models/LIF')

% Number of excitatory and inhibitory neurons, eNnrn and iNnrn, and total
  % # of neurons, totNnrn:  use 4:1 ratio
  net_COBN.eNnrn = 800; %800; %400; %40; %40; %400; % 240
  net_COBN.iNnrn = 200; %200;  %100; % 60
 
% Connectivity, i.e., connection probability, p:
  net_COBN.p = 0.4;  %net_COBN.p = 0.2;   % Careful ! Not too large, if n small
  frnLinks= net_COBN.p % check
  clear frnLinks
% >>>>>>>>>>>>>>>>
% Membrane time constant, [ms]
  net_COBN.eTm = 20; 
  net_COBN.iTm = 10; 
 
% Leak membrane potential, [mV]
  net_COBN.V_leaky = -70.; 
 
% Membrane potential firing threshold, [mV]
  net_COBN.Vthr = -52; 
% the neuron fires following the following sequence:
 % 1. the neuron potential is reset to a value Vres, [mV]: 
  net_COBN.eVres = -59; 
  net_COBN.iVres = -59; 
 % 2. the neuron cannot fire again for a refractory period, Trp, equal to
 %    2ms for E neurons and 1ms for I neurons, [ms]:
  net_COBN.eTrp = 2; 
  net_COBN.iTrp = 1; 
 % 3. all post-synaptic neurons receive a spike with a delay, Tl, equal to
 %    of 1ms after the time of threshold crossing, [ms]:
  net_COBN.eTl = 1;
  net_COBN.iTl = 1; 
 
% Rise and deacy times, Tr and Td [ms]:
% - of E => E synaptic currents:
  net_COBN.e2eTr = 0.4;
  net_COBN.e2eTd = 2.; 
% - of E => I synaptic currents:
  net_COBN.e2iTr = 0.2; 
  net_COBN.e2iTd = 1.; 
% - of I synaptic currents: rise and decay times are assumed idependent of
%   the type of neuron they act on (excitatory or inhibitory)
  net_COBN.iTr = 0.25; 
  net_COBN.iTd = 5.; 
 
% Synaptic reversal potentials, [mV]
  net_COBN.VsynAMPA = 0; 
  net_COBN.VsynGABA = -80;
 
% Synaptic conductances [nS]
% - on inhibitory neurons:
  net_COBN.gi2i = 2.698602679456193;
  net_COBN.ge2i = 0.233373613159943; 
  net_COBN.gx2i = 0.316721332145637; 
% - on excitatory neurons:
  net_COBN.gi2e = 2.008771996214003; 
  net_COBN.ge2e = 0.178441556321102; 
  net_COBN.gx2e = 0.233673466610967; 
 
% Membrane resistances, [GOhm]
  net_COBN.eRm = 0.04; 
  net_COBN.iRm = 0.05; 
 
% > >   >   >

%% 1.1 other sim parameters
  fprintf('\n COBN sim: set basic sim parameters \n')
 % need to execute this In the Matlab workspace:
 % parameters_COBN; % above: generates the struct net_COBN with all the parameters of the network
  %simulation_length = 20000; % units: (msec). ie 10, 20, 40 sec
  simulation_length = xf * 1000; % units: (msec). ie 10, 20, 40 sec
  M = simulation_length/net_COBN.Dt; % length of the simulation in time steps (Dt)
 
  external_signal_intensity = 2 %default % units; (spikes/ms/cell)
 % >>>>>>   Default: const stim           >>>>>>>>>>>
  external_signal = ones(M,1) * external_signal_intensity * net_COBN.Dt; % const (200 Hz default)
  SEED_OU = 1; % positive integer number
  % tmp=clock;  % go to #1.1a  vary seed via clock 
  external_noise = OU_process(M, net_COBN.Dt, 16, 0.16*net_COBN.Dt, SEED_OU);
 %INPUT2E = external_signal + external_noise; % const + OU noise
  INPUT2E = external_signal; % const ext stim only
   % ALSO: apply extra ext stim, set up at #1.1c, above
    %extra_signal = waveStim*net_COBN.Dt; % modulated by wave {from # 1.1c
    %INPUT2E = INPUT2E -min(INPUT2E); % ensure >=0
  INPUT2I = INPUT2E; % default - stim to i is same as basic stim for e [local]
  %INPUT2E = INPUT2E + extra_signal; % add n extra (modulated, etc) stim
   %INPUT2E = INPUT2E -min(INPUT2E); % ensure >=0 rate, if neg net stim present?
  SEED_connections = 2; % positive integer number; const  for reproducible ran samples
  SEED_poisson = 3; % positive integer number
 
 % & variables for other codes: nb. output vectorss @ 1ms bins!
  h = net_COBN.Dt*1e-3; %(now in msec; nb Dt in ms)
  %x=[0: 1e-4: simulation_length*1e-3]'; x=x(2:end); % t axis (sec) [for fft]
  %TotTime = max(x) % nb. need to match lengths, so omit t=0 
  Ntot = net_COBN.eNnrn + net_COBN.iNnrn % tot # neurons
 
  figure; plot(x, external_noise); title('ext OU noise stimulus'); xlabel('t (sec)'); % Debug
  figure; plot(x, external_signal); title('ext const. stimulus'); xlabel('t (sec)');
   % figure; plot(x, extra_signal); title('ext modulated stimulus');
   % xlabel('t (sec)'); % if present?
  figure; plot(x, INPUT2E); xlabel('t (sec)'); title('tot ext stimulus: const + OU noise');
   title('tot ext stimulus: const + theta modulated + OU noise');
  maxInput = max(INPUT2E)
  minInput = min(INPUT2E)
  clear TotTime maxInput minInput

%% 1.1b other ext stim
  % eg. sin wave
  fprintf('\n COBN sim: add ext wave \n')
  Aampl = 0.6; % (spikes/ms) : amplitude of ext wave
  wfreq = 7.7 % {Hz) : ~2-15-150 used in orig paper
  phi=0; % phase (rad)
 % time vector, x set above: steps of Dt ~1e-4, from "1" to 10, 20s [200k elements)
  %tline = x*10000.0; % convert ms scale to sec, to match pulses & DE solver
  %max(tline)
  waveStim = Aampl*sin(2*pi*wfreq*x -phi); % base case;  row vec
  figure; plot(x, waveStim); ylim([-1 1]); title('ext wave stimulus')
%%%% INPUT2E = INPUT2E+ waveStim;
  figure; plot(x, INPUT2E); %ylim([-1 1]); 
  title('total ext stimulus')
  clear Aampl wfreq phi  

% "Standard" initialisation done

  %  Initialise a vector for saved membrane potentials
  %V = zeros(Ntot,1) + net_COBN.V_leaky;
  V = ones(1, Ntot) * net_COBN.V_leaky;

  %  Initialise a vector for time of last spike
  %tLastSP = zeros(Ntot,1) - M;
  tLastSP = - ones(1, Ntot) * simulation_length;

  % Initialise array of neuron connections
  A = zeros(Ntot);

  fprintf("Setup done!\n");
  pause(1);


  % iterations should be set to total sample length divided by our sample_width
  iterations = M / sample_width;
  fprintf("iterations set by: %d/%d = %d\n", M, sample_width, iterations);

  [E2EI,I2EI,eFR,iFR] = code_COBN(net_COBN, INPUT2E, INPUT2I, SEED_connections, SEED_poisson);

  for i = 2:iterations
    fprintf("\nloop_COBN() Iteration = %d of %d\n", i,iterations);

    net_COBN.restart = net_COBN.restart + sample_width;
    [E2EI,I2EI,eFR,iFR] = code_COBN(net_COBN, INPUT2E, INPUT2I, SEED_connections, SEED_poisson);
  end

  figure;
  tiledlayout(4,1)
  ax1 = nexttile; plot(E2EI(1:end), 'Color', 'black');	title('E2EI'); xlabel('time');
  ax2 = nexttile; plot(I2EI(1:end), 'Color', 'blue');	title('I2EI'); xlabel('time');
  ax3 = nexttile; plot(eFR(1:end), 'Color', 'yellow');	title('eFR');  xlabel('time');
  ax4 = nexttile; plot(iFR(1:end), 'Color', 'magenta');	title('iFR');  xlabel('time');
  linkaxes([ax1 ax2 ax3 ax4],'x');

end


