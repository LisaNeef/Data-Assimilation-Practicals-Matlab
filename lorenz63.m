%--------------------------------------------------------------
function xout = lorenz63(xin, sigma, rho, beta, dt)
%     
%     Integrate Lorenz 1963 model by one time-step.
%     Time integration using 4th order Runge-Kutta scheme.  
%     See Numerical Recipes ch. 15.1
%     Lisa Neef, 14 June 2013
%     based on code by S. Polavarapu 24 March 2002
%
% INPUTS:
%  xin: the model state to be advanced
%  rho, sigma, b: model parameters
%  dt: the time step
%
% OUTPUTS:
%  xout: the model state advanced forward 1 timestep


%-----------------
nv = size(xin,1);
xx = zeros(nv,1);
x1 = zeros(nv,1);
x2 = zeros(nv,1);
x3 = zeros(nv,1);
x4 = zeros(nv,1);
fp = zeros(nv,1);
w1 = 1.0/6.0;
w2 = 1.0/3.0;
w3 = 1.0/3.0;
w4 = 1.0/6.0;

%     Mean trajectory calculation
xx1 = xin;
fp  = lorenz63_rhs(xx1,sigma,rho,beta);
x1  = dt*fp;
xx2 = xin+ 0.5*x1;
fp  = lorenz63_rhs(xx2,sigma,rho,beta);
x2  = dt*fp;
xx3 = xin+ 0.5*x2;
fp  = lorenz63_rhs(xx3,sigma,rho,beta);
x3  = dt*fp;
xx4 = xin+ x3;

xout = xin + w1*x1 + w2*x2 + w3*x3 + w4*x4;
