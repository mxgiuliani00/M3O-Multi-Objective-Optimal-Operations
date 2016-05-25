%% Multi-Objective Optimal Operation (M3O) Toolbox
%
% This toolbox purpose is to design the optimal operations of multipurpose 
% water reservoir system, in an attempt to close the gap between research 
% studies and real-world applications. The toolbox, called Multi-Objective 
% Optimal Operations (M3O), allows users to design Pareto optimal (or 
% approximate) operating policies through several alternative state-of-the-art
% methods. The application of all these techniques on the same case study 
% contributes a step-forward with respect to traditional literature review 
% papers as the availability of the source code, along with the possibility 
% of cross-comparing the results on the same problem, allows a better 
% understanding of pros and cons of each approach. At the same time, the
% modular structure of M3O allows experienced users to easily customize, and
% possibly further develop, the implemented code according to their specific
% requirements.


clear all;
clc

global sys_param;

addpath(genpath('sim'))
addpath(genpath('lib'))

%% Configure general system parameters

sys_param.simulation.q     = load('inflow.txt','-ascii') ;
sys_param.simulation.h_in  = 0.6 ;
sys_param.simulation.w     = 370 ;    % irrigation demand
sys_param.simulation.hFLO  = 0.8 ;    % m the flooding threshold
sys_param.simulation.h0    = -0.5;
sys_param.simulation.A     = 145900000;
sys_param.simulation.r_min = 0;
sys_param.simulation.r_max = 518;
sys_param.simulation.delta = 60*60*24;

[sys_param.simulation.vv, sys_param.simulation.VV] = deal(0); 

%% --- Run DDP (Diterministic Dynamic Programming) ---
clc
addpath('./DP')

% Configure the parameters 
load 'grids.mat';

sys_param.algorithm = grids;
sys_param.algorithm.name = 'ddp';
sys_param.algorithm.Hend = 0 ; % penalty set to 0

[vv, VV] = construct_rel_matrices(grids); % compute daily minimum/maximum release matrixes
sys_param.algorithm.min_rel = vv;
sys_param.algorithm.max_rel = VV;


% weights for aggregation of objectives
lambda = [1 0; .75 .25; .5 .5 ; .35 .65; .2 .8; .1 .9;  0 1];
Nalt   = size(lambda,1);
JJ_ddp = nan(Nalt,2);
Hddp   = cell(Nalt,1);

for i = 1: Nalt
  sys_param.algorithm.weights = lambda(i,:);
  [JJ_ddp(i,:), Hddp{i}] = run_ddp() ;
end

figure; plot( JJ_ddp(:,1), JJ_ddp(:,2), 'o' );
xlabel('flooding'); ylabel('irrigation');

%% --- Run SDP (Stochastic Dynamic Programming) ---
clc
addpath('./SDP')

% Configure the parameters 
load 'grids.mat';

sys_param.algorithm = grids;
sys_param.algorithm.name = 'sdp';
sys_param.algorithm.Hend = 0 ; % penalty set to 0
sys_param.algorithm.T = 1 ;    % the period is equal 1 as we assume stationary conditions

[vv, VV] = construct_rel_matrices(grids); % compute daily minimum/maximum release matrixes
sys_param.algorithm.min_rel = vv;
sys_param.algorithm.max_rel = VV;

% Estimate inflow probability density function assuming log-normal
% distribution fitting
sys_param.algorithm.q_stat = lognfit(sys_param.simulation.q);
sys_param.algorithm.gamma = 1; % set future discount factor

lambda = [1 0; .75 .25; .5 .5 ; .35 .65; .2 .8; .1 .9;  0 1];
Nalt   = size(lambda,1);
JJ_sdp = nan(Nalt,2);
Hsdp   = cell(Nalt, 1);

for i = 1: Nalt
  sys_param.algorithm.weights = lambda(i,:);
  [JJ_sdp(i,:), Hsdp{i}] = run_sdp();
end

figure; plot( JJ_sdp(:,1), JJ_sdp(:,2), 'o' );
xlabel('flooding'); ylabel('irrigation');

%% --- Run EMODPS (Evolutionary Multi-Objective Direct Policy Search) ---
clc
addpath('./EMODPS')

% Define the parameterized class for the policy (i.e., standard operating
% policy)
sys_param.algorithm.name = 'emodps' ;
pClass = 'stdOP'; 

% Define MOEA and its setting (i.e., NSGAII)
moea_param.name = 'NSGAII';
moea_param.pop  = 40; % number of individuals
moea_param.gen  = 50; % number of generation

[JJ_emodps, Popt] = run_emodps(pClass, moea_param) ;

figure; plot( JJ_emodps(:,1), JJ_emodps(:,2), 'o' );
xlabel('flooding'); ylabel('irrigation');

%% --- Run FQI (Fitted Q-Iteration) --- 
clc
addpath('./FQI')

load 'grids.mat';
sys_param.algorithm = grids;

% Construction of sample dataset by Monte Carlo simulation 
sys_param.algorithm.name = 'doe';
[F, G] = run_doe(50) ; % create 50 samples of tuples 

% Aggregation of the immediate costs
lambda = [1 0; .75 .25; .5 .5 ; .35 .65; .2 .8; .1 .9;  0 1];
Nalt   = size(lambda,1);
JJ_fqi = nan(Nalt, 2);
Qfqi   = cell(Nalt, 1);

% Define regressor and its parameters (i.e. extra-trees). 
reg_param.name    = 'ET';
reg_param.M       = 200;  % number of trees
reg_param.nmin    = 25;   % minimum number of points per leaf
reg_param.maxIter = 40;

sys_param.algorithm.name = 'fqi';
sys_param.algorithm.gamma = 0.99;

for i = 1: Nalt
  sys_param.algorithm.weights = lambda(i,:);
  [JJ_fqi(i,:), Qfqi{i}] = run_fqi(F, G, reg_param);
end

figure; plot( JJ_fqi(:,1), JJ_fqi(:,2), 'o' );
xlabel('flooding'); ylabel('irrigation');


%% --- Run MPC (Model Predictive Control) --- 
clc;
addpath('./MPC')

% MPC settings
sys_param.algorithm.name = 'mpc';
sys_param.algorithm.P = 3;            % Length of the moving prediction horizon
sys_param.algorithm.mi_e    = mean(sys_param.simulation.q);
sys_param.algorithm.sigma_e = std(sys_param.simulation.q);

% mpc_input  = sys_param.simulation.q;  % Candidate disturbance variable to be predicted for MPC
errorLevel = 0;                       % Disturbance prediction error [%] 

% Weights value for aggregation of objectives
lambda = [1 0; .75 .25; .5 .5 ; .35 .65; .2 .8; .1 .9;  0 1];
Nalt   = size(lambda,1);
JJ_mpc = nan(Nalt, 2);
Ompc   = cell(Nalt, 1);

for i = 1: Nalt
  sys_param.algorithm.weights = lambda(i,:);
  [JJ_mpc(i,:), Ompc{i}] = run_mpc(errorLevel);
end

figure; plot( JJ_mpc(:,1), JJ_mpc(:,2), 'o' );
xlabel('flooding'); ylabel('irrigation');

%% --- Run ISO ( Implicit Stochastic Optimization) ---
clc;
addpath('./DP')

% Configure the parameters 
load 'grids.mat';

sys_param.algorithm = grids;
sys_param.algorithm.name = 'iso';

[vv, VV] = construct_rel_matrices(grids); % compute daily minimum/maximum release matrixes
sys_param.algorithm.min_rel = vv;
sys_param.algorithm.max_rel = VV;
sys_param.algorithm.Hend    = 0 ; % penalty set to 0

% Define regression method
regressor = 'linear_spline';
sys_param.algorithm.regressorName = regressor;

% weights value for aggregation of objectives
lambda = [1 0; .75 .25; .5 .5 ; .35 .65; .2 .8; .1 .9;  0 1];
Nalt = size(lambda,1);

[JJ_iso, err_perc] = deal(nan(Nalt,2));
policy = cell(Nalt,1);

for i = 1: Nalt
  sys_param.algorithm.weights = lambda(i,:);
  [JJ_iso(i,:), policy{i}, err_perc(i,:)] = run_iso(regressor);
end

% plot
figure; plot( JJ_iso(:,1), JJ_iso(:,2), 'o' );
xlabel('flooding'); ylabel('irrigation');

%% --- Run SSDP (Sampling Stochastic Dynamic Programming) --- 
clc;
addpath('./SSDP')

% Configure the parameters 
load 'grids.mat';

sys_param.algorithm = grids;
sys_param.algorithm.name = 'ssdp';
sys_param.algorithm.Hend = 0 ; % penalty set to 0

[vv, VV] = construct_rel_matrices(grids); % compute daily minimum/maximum release matrixes
sys_param.algorithm.min_rel    = vv;
sys_param.algorithm.max_rel    = VV;
sys_param.algorithm.T          = 365;
sys_param.algorithm.interp_foo = @interp1qr; % default interpolator 
sys_param.algorithm.gamma      = 1;
sys_param.algorithm.Pr_mode    = 2; 
sys_param.algorithm.cycle_T    = 300;
sys_param.algorithm.forecast_T = 60;

% estimate inflow pdf (log-normal distribution) and create the test samples
% of total 25 ensembles with 365 days each
q_stat = lognfit(sys_param.simulation.q) ;

esp_sample = lognrnd( q_stat(1), q_stat(2), 365, 25 );
sys_param.algorithm.esp_sample = esp_sample ;

% weights value for aggregation of objectives
lambda = [1 0; .75 .25; .5 .5 ; .35 .65; .2 .8; .1 .9; 0 1];
Nalt = size(lambda,1);

JJ_ssdp  = nan(Nalt, 2);
H_policy = cell(Nalt, 1);

for i = 1: Nalt
  sys_param.algorithm.weights = lambda(i,:);
  [JJ_ssdp(i,:), H_policy{i}] = run_ssdp() ;
end

figure; plot( JJ_ssdp(:,1), JJ_ssdp(:,2), 'o' );
xlabel('flooding'); ylabel('irrigation');
