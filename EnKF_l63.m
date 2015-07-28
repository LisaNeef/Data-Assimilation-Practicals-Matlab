function [XT,XA,t,ET] = enkf_l63(tobs,Tend,N)
%% EnKF_l63.m
%
% Run the Lorenz 1963 model and 
% assimilate observations using the Ensemble Kalman Filterq

% Lisa Neef, 17 June 2013
%
% INPUTS:
%   tobs: the observation interval
%   Tend: the model integration time
%   N: ensemble size
%
% OUTPUTS:
%   XT: the 3-variable true state
%   XA: the 3-variable analysis state (the mean of the ensemble)
%   t: the time array
%   ET: the true error
%-----------------------------------------------------------

%% Other assimilation Parameters:

sig0 	= 0.5;	% initial forecast error variance
sig_obs	= 0.01;	% observation errror covariance

%% Model Parameters

sigma 	= 10;
rho	= 28;
beta	= 8/3;
dt 	= 0.01;

%% Initial conditions

% initial true state
xt0 = zeros(3,1);
xt0(1) =  1.508870;
xt0(2) = -1.531271;
xt0(3) = 25.46091 ;

xf0 = xt0+sig0*rand(3,1);

%% Initialize arrays to hold everything
t = 1:dt:Tend;
nT = length(t);
XT = zeros(3,nT)+NaN;	% array to hold the true state in time
XENS = zeros(N,3,nT)+NaN;	% array to hold the ensemble
S  = zeros(3,nT)+NaN;	% array to hold analysis error variance in time
YOBS = zeros(3,nT)+NaN;

%% initial conditions
XT(:,1) = xt0;		% initial truth
xf = xf0;		% initial forecast
Pf = sig0*eye(3);	%

%% generate the initial ensemble 
for iens = 1:N
  XENS(iens,:,1) = xf0 + sig0*randn(3,1);
end

% this is the initial forecast ensemble
xfens = squeeze(XENS(:,:,1));

%% Define other matrices needed in the assimilation
H = eye(3);
R = sig_obs*eye(3);


%% Loop in time

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


%% Compute a few other output quantities
XA = squeeze(mean(XENS,1));

%% Produce Plots!

YL = {'x','y','z'};


% plot the state analysis versis truth
figure(1),clf
h = zeros(1,4);  % legend handle
T = ones(N,1)*t;
for ic = 1:3
  subplot(3,1,ic)
    dum  = plot(T,squeeze(XENS(:,ic,:)),'Color',0.7*ones(1,3));
    hold on
    h(1) = dum(1);
    h(2) = plot(t,XT(ic,:),'k');
    h(3) = plot(t,XA(ic,:),'b');
    h(4) = plot(t,YOBS(ic,:),'ro');
    xlabel('time')
    ylabel(YL(ic))
    legend(h, 'ensemble','truth','analysis','obs')
    title('Lorenz 1963 Model')
end

% plot the estimated and true errors
ET = sqrt((XT-XA).^2);
EA = sqrt(S);
figure(2),clf
for ic = 1:3
  subplot(3,1,ic)
    semilogy(t,ET(ic,:),'k',t,EA(ic,:),'b')
    xlabel('time')
    ylabel(YL(ic))
    legend('true error','analysis error')
    title('Errors')
end

fig_name_1 = ['lorenz_EnKF_tobs',num2str(tobs),'_N',num2str(N),'.png'];
fig_name_2 = ['lorenz_EnKF_tobs',num2str(tobs),'_N',num2str(N),'_error.png'];
pw1 = 10;
ph1 = 10;


%% Export Plots

exportfig(1,fig_name_1,'width',pw1,'height',ph1,'format','png','color','cmyk','FontSize',1)
exportfig(2,fig_name_2,'width',pw1,'height',ph1,'format','png','color','cmyk','FontSize',1)



