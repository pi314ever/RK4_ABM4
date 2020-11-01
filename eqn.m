function [dx] = eqn(t,x,P)
%eqn    Equations of motion for an N-pendulum 

% System parameters
g = P.g; % Gravity
m = P.m; % Mass vector
L = P.l; % Massless rod length vector
n = length(m); % Number of masses
if n ~= length(L) % Number of masses must equal number of rods
    error('Mismatched vector lengths for m and L')
end

% Initialize % Definitions
dx = zeros(2*n,1); % Output 
M = zeros(n,1); % Cumulative mass vector
A = zeros(n); % Coefficient matrix
f_vec = zeros(n,1); % f vector function
g_vec = zeros(n,1); % g vector function

for ii = 1:n
    M(ii) = sum(m(ii:n));
end

for ii = 1:n
    for jj = 1:n
        A(ii,jj) = L(ii)*L(jj)*cos(x(ii)-x(jj))*M(max([ii jj]));
    end
end

for k = 1:n
    for jj = k:n
        for ii = 1:jj
            f_vec(k) = f_vec(k)+m(jj)*L(k)*L(ii)*x(n+ii)*sin(x(k)-x(ii))*...
                (x(n+ii)-x(n+k));
            g_vec(k) = g_vec(k)+L(k)*x(n+k)*m(jj)*L(ii)*x(n+ii)*...
                sin(x(ii)-x(k));
        end
    end
end
g_fun_PE = g*L.*sin(x(1:n)).*M;
g_vec = g_vec - g_fun_PE; % Adjust for PE partial derivative

h_vec = A\(g_vec-f_vec); % Solve for thetadd using linear algebra 

% Outputs
dx(1:n) = x(n+1:2*n); % d/dt theta = thetad
dx(n+1:2*n) = h_vec; % Import remaining dx values
end

