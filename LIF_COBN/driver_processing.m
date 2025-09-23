
% To provide a choice betwen connection matrix generators
% (the default being connections_COBN.c called by driver_connections.m)
% let's assume we have an alternate generator named "Aconvert",
% whose output is Aconn. If Aconn is empty (Aconvert wasn't used)
% we rename A (output of connections_COBN.c) as Aconn to be passed
% as input parameter to process_COBN

Aconn = [];
%Aconvert(); % leaves new connection matrix as Aconn

% In case Aconvert wasn't used (so Aconn is empty)
% run driver_connections to obtain A and assign Aconn to A
if isempty(Aconn) 
  % Generate connections, initial net_COBN and A
  driver_connections();
  Aconn = A;
end

% Display which A will be used
figure; spy(Aconn); title('Using Aconn as A')


% Compile processor
mex process_COBN.c ran1.c gasdev.c


% !!! Be careful not to overwrite !!!
% !!! existing net_COBN settings made above !!!

net_COBN.restart = 0;
net_COBN.sample_width = 128;

% Additions to net_COBN for processing
%  time vec
FactorInner = 20 ; %50 % 200 ms  / modulo 20 sec 
x0=0; xf= 20.0/FactorInner ;%1, 20.0 (s); % t range (sec)  % 10 x 2sec
h = 0.0001;   % t step (1 ms; [& 0.5, 0.25 ms, for debug])
x = [x0:h:xf];  % t domain array
x=x(2:end); % match vector length, if needed
simL_sec =max(x) ;
clear simL*  FactorInner

% Time resolution [ms]
net_COBN.Dt = 0.1; %0.05; % 0.1 ms t step

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

%% 1.1 other sim parameters; ext noise vecor[* t length]
fprintf('\n COBN sim: set basic sim parameters; ext noise \n')
% need to execute this In the Matlab workspace:
% parameters_COBN; % above: generates the struct net_COBN with all the parameters of the network
%simulation_length = 20000; % units: (msec). ie 10, 20, 40 sec
simulation_length = xf * 1000; % units: (msec). ie 10, 20, 40 sec
M = simulation_length/net_COBN.Dt; % length of the simulation in time steps (Dt)

external_signal_intensity = 2 ; %default % units; (spikes/ms/cell)
% >>>>>>   Default: const stim           >>>>>>>>>>>
external_signal = ones(M,1) * external_signal_intensity * net_COBN.Dt; % const (200 Hz default)
SEED_OU = 1; % positive integer number
% tmp=clock;  % go to #1.1a  vary seed via clock 
external_noise = OU_process(M, net_COBN.Dt, 16, 0.16*net_COBN.Dt, SEED_OU);
INPUT2E = external_signal + external_noise; % const + OU noise
%INPUT2E = external_signal; % const ext stim only
% ALSO: apply extra ext stim, set up at #1.1c, above
%extra_signal = waveStim*net_COBN.Dt; % modulated by wave {from # 1.1c
%INPUT2E = INPUT2E -min(INPUT2E); % ensure >=0

INPUT2I = INPUT2E; % default - stim to i is same as basic stim for e [local]
%INPUT2E = INPUT2E + extra_signal; % add n extra (modulated, etc) stim
%INPUT2E = INPUT2E -min(INPUT2E); % ensure >=0 rate, if neg net stim present?

% & variables for other codes: nb. output vectorss @ 1ms bins!
h = net_COBN.Dt*1e-3; %(now in msec; nb Dt in ms)
%x=[0: 1e-4: simulation_length*1e-3]'; x=x(2:end); % t axis (sec) [for fft]
%TotTime = max(x) % nb. need to match lengths, so omit t=0 
Ntot = net_COBN.eNnrn + net_COBN.iNnrn; % tot # neurons

%figure; plot(x, external_noise); title('ext OU noise stimulus'); xlabel('t (sec)'); % Debug
%figure; plot(x, external_signal); title('ext const. stimulus'); xlabel('t (sec)');
% figure; plot(x, extra_signal); title('ext modulated stimulus');
% xlabel('t (sec)'); % if present?
figure; plot(x, INPUT2E); xlabel('t (sec)'); title('tot ext stimulus: const + OU noise');
title('tot ext stimulus: const  + OU noise');
maxInput = max(INPUT2E);
minInput = min(INPUT2E);
clear TotTime maxInput minInput

% iterations should be set to total sample length divided by our sample_width
iterations = M / net_COBN.sample_width;
fprintf("iterations set by: %d/%d = %f\n", M, net_COBN.sample_width, iterations);


% Dummy inputs for first iteration, where arrays will be allocated
tLastSP = [];
aX_rec  = [];
aX_ext  = [];
gX      = [];
aS_rec  = [];
aS_ext  = [];
gS      = [];
V       = [];
E2EI    = [];
I2EI    = [];
 eFR    = [];
 iFR    = [];
eSPC    = [];
iSPC    = [];


% iterate forward in time
for i = 1:iterations
    fprintf("\nprocess_COBN() Iteration = %d of %d\n", i, iterations);

    [E2EI,I2EI,eFR,iFR,tLastSP,aX_rec,aX_ext,gX,aS_rec,aS_ext,gS,V,eSPC,iSPC] = process_COBN(net_COBN, INPUT2E, INPUT2I, tLastSP, aX_rec, aX_ext, gX, aS_rec, aS_ext, gS, V, Aconn, E2EI, I2EI, eFR, iFR, eSPC, iSPC);

    net_COBN.restart = net_COBN.restart + net_COBN.sample_width;
end


figure;  tiledlayout(4,1)
ax1 = nexttile; plot(E2EI(1:end), 'Color', 'blue');   
title('E2EI'); xlabel('time steps (0.1 ms)');
ax2 = nexttile; plot(I2EI(1:end), 'Color', 'red'); title('I2EI(t) 2 sec @ 128  steps'); xlabel('time');
ax3 = nexttile; plot(eFR(1:end), 'Color', 'cyan');    title('eFR');  xlabel('time');
ax4 = nexttile; plot(iFR(1:end), 'Color', 'magenta'); title('iFR');  xlabel('time');
linkaxes([ax1 ax2 ax3 ax4],'x');

