function enkf_l63(E)
%% EnKF_l63.m
%
% Run the Lorenz 1963 model and 
% assimilate observations using the Ensemble Kalman Filter

% Lisa Neef; started 17 June 2013
%
% INPUT:
%   E: a Matlab structure that holds all the model and assimilation
%	parameters. It is generated using E = set_enkf_inputs
%-----------------------------------------------------------


%% extract the individual assimilation parameters from E
sigma 	= E.sigma;
rho 	= E.rho;
beta 	= E.beta;
dt 	= E.dt;
xt0 	= E.xt0;
xf0 	= E.xf0;
sig0 	= E.sig0;
sig_obs = E.sig_obs;
N 	= E.N;
Tend 	= E.Tend;
tobs 	= E.tobs;

%% Initialize arrays to hold everything
t 	= 1:dt:Tend;
nT 	= length(t);
XT 	= zeros(3,nT)+NaN;	% array to hold the true state in time
XENS 	= zeros(N,3,nT)+NaN;	% array to hold the ensemble
S  	= zeros(3,nT)+NaN;	% array to hold analysis error variance in time
YOBS 	= zeros(3,nT)+NaN;

%% initial conditions
XT(:,1) = xt0;			% initial truth
xf 	= xt0+sig0*rand(3,1);	% initial forecast
for iens = 1:N
  XENS(iens,:,1) = xf0 + sig0*randn(3,1);
end

%% Define other matrices needed in the assimilation
H = eye(3);
R = sig_obs*eye(3);

%% Loop in time
xfens = squeeze(XENS(:,:,1));
for k = 1:nT-1

  % observations?
  if (mod(t(k),tobs) == 0)
    YOBS(:,k) = XT(:,k)+sig_obs*rand(3,1);

    % get the forecast error covariance matrix from the ensemble
    D = zeros(3,N);
    for iens = 1:N
      D(:,iens) = xfens(iens,:) - mean(xfens,1);
    end
    Pf = (1/N)*D*D';

    % update the entire ensemble with the observations
    K_bottom = H*Pf*H' + R;
    K_top = H*Pf';
    K = K_top * inv(K_bottom);
    for iens = 1:N
      XENS(iens,:,k) = xfens(iens,:)' + K*(YOBS(:,k) - H*xfens(iens,:)');
    end

  else
    % if no observation, then the forecast becomes the analysis
    XENS(:,:,k) = xfens;
  end

  % regardless of whether there's been an observation, the analysis error covariance matrix comes from the ensemble
    D = zeros(3,N);
    for iens = 1:N
      D(:,iens) = XENS(iens,:,k) - mean(XENS(:,:,k),1);
    end
    Pa = (1/N)*D*D';
 

  % save the diagonals of the analysis error covariance matrix -- these are the variances
  S(:,k) = diag(Pa);


  % evolve the truth forward
  XT(:,k+1) = lorenz63(XT(:,k), sigma, rho, beta, dt);

  % evolve the analysis ensemble forward to become the next forecast
  xfens = zeros(N,3);
  for iens = 1:N
    xfens(iens,:)  = lorenz63(XENS(iens,:,k), sigma, rho, beta, dt);
  end

end

%---------------PLOTTING----------------------------------

%% Compute a few other output quantities
XA = squeeze(mean(XENS,1));

%% Produce Plots!
YL = {'x','y','z'};

%% definte some plot settings
LW = 2;		% line width
tcol = [0,0,0];
acol = [217,95,2]/256.0;
acol = [102,166,30]/256.0;
ocol = [241,41,138]/256.0;
ecol = .7*ones(1,3);


% plot the state analysis versis truth
figure(1),clf
h = zeros(1,4);  % legend handle
T = ones(N,1)*t;
for ic = 1:3
  subplot(3,1,ic)
    ens = transpose(squeeze(XENS(:,ic,:)));
    dum  = plot(transpose(T),ens,'Color',ecol,'LineWidth',1);
    hold on
    h(1) = dum(1);
    h(2) = plot(t,XT(ic,:),'Color',tcol,'LineWidth',LW);
    h(3) = plot(t,XA(ic,:),'Color',acol,'LineWidth',LW);
    h(4) = plot(t,YOBS(ic,:),'o','Color',ocol,'MarkerSize',5,'LineWidth',LW);
    if ic == 1
	    title('Lorenz 1963 Model - State Variables')
    end
    if ic == 3
	    xlabel('time')
	    legend(h, 'ensemble','truth','analysis','obs','Location','SouthOutside','Orientation','Horizontal')
    end
    ylabel(YL(ic))
end

% plot the estimated and true errors
ET = sqrt((XT-XA).^2);
EA = sqrt(S);
figure(2),clf
for ic = 1:3
  subplot(3,1,ic)
    h = zeros(1,2);
    h(1) = semilogy(t,ET(ic,:),'Color',tcol,'LineWidth',LW);
    hold on
    h(2) = semilogy(t,EA(ic,:),'Color',acol,'LineWidth',LW);
    ylabel(YL(ic))
    if ic == 1
	    title('Lorenz 1963 Model - RMSE')
    end
    if ic == 3
	    xlabel('time')
	    legend(h, 'true RMSE','estimated RMSE','Location','SouthOutside','Orientation','Horizontal')
    end
end

% average the errors for each variable over the integration time, and print to the screent
ETave = nanmean(ET,2);
EAave = nanmean(EA,2);
names = {'x','y','z'};
disp('+++++++Assimilation Run Average Errors++++++++++')
for ic = 1:3
	str_out = strcat(names(ic),':	True Error = ',num2str(ET(ic),3),'  Estimated = ',num2str(EA(ic),3));
	disp(str_out)
end
