function E = set_enkf_inputs
%% set_enkf_inputs.m  
%
% Set the input parameters for the Lorenz 1963 model, run within an Ensemble Kalman Filter
% All the inputs are put together in a matlab structure, E, wich is returned by this code. 
% E then becomes the input for EnKF_l63.m.
%
% Code by Lisa Neef, 1 August 2015
%---------------------------------------------------------------------------------

%-------------modify this part----------------------------------------------
%% MODEL PARAMETERS 
sigma 	= 10;
rho	= 28;
beta	= 8/3;
dt 	= 0.01;

%% ASSIMILATION PARAMETERS
run_filter 	= 0;	% set to 1 to run the EnKF
sig0 		= 0.5;	% initial forecast error variance
sig_obs		= 0.5;	% observation errror covariance
N		= 10;	% ensemble size
Tend		= 30;	% total integration time
tobs 		= 3.0;	% observation interval 
localize 	= 0;	% set to 1 to use only the diagonals of the covariance matrix

%% OBSERVATION SETTINGS
obsx 	= 1;	% set to 1 to observe variable x
obsy 	= 1;	% set to 1 to observe variable y
obsz 	= 1;	% set to 1 to observe variable z

obs_meanxy	= 0;	% set to 1 to observe the mean of x and y
obs_meanyz	= 0;	% set to 1 to observe the mean of y and z
obs_meanxz	= 0;	% set to 1 to observe the mean of x and z

%% INITIAL CONDITIONS
xt0 = zeros(3,1);
xt0(1) = 10*randn(1);	% initial true value for x
xt0(2) = 10*randn(1);	% initial true value for y
xt0(3) = 10*randn(1);	% initial true value for z

xf0 = zeros(3,1);
xf0 = xt0 + sig0*randn(3,1);	% initial vector for forecast


%-------------do not modify this part----------------------------------------------
%% put it all together as a structure
E = struct('sigma',sigma,...
		'rho',rho,...
		'beta',beta,...
		'dt',dt,...
		'obsx',obsx,...
		'obsy',obsy,...
		'obsz',obsz,...
		'obs_meanxy',obs_meanxy,...
		'obs_meanyz',obs_meanyz,...
		'obs_meanxz',obs_meanxz,...
		'xt0',xt0,...
		'xf0',xf0,...
		'run_filter',run_filter,...
		'sig0',sig0,...
		'sig_obs',sig_obs,...
		'N',N,...
		'Tend',Tend,...
		'localize',localize,...
		'tobs',tobs);
		
