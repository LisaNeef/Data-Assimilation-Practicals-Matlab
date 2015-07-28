%--------------------------------------------------------------
function f = lorenz63_rhs (vv,sigma,rho,beta)

%  Gets the right hand side of Lorenz's chaos butterfly model
%
%  Lisa Neef, 14 June 2013

nv = size(vv,1);
f = zeros(nv,1);

% --- Rename input variables

x = vv(1);
y = vv(2);
z = vv(3);

% --- RHS of eq. (5) of Lorenz (1986) p1550

f(1) = sigma*(y-x);
f(2) = rho*x-y-x*z;
f(3) = x*y-beta*z;
