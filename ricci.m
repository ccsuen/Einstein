function [] = ricci(m, e, c, a, N)
% Program that numerically propogates the ricci soliton equations
% Created by Cissy Suen for Dr. Mckenzie Wang, Summer 2015
%   Creates a Taylor series around the origin
%   Produced initial values taken as initial values of Runge Kutta
%   Runge Kutta propogated to create plots of behaviour of the space
%   m is the parameter experimentally set
%   e, or epsilon,
%   c
%   a
%   N is the number of terms desired for the initial Taylor polynomial

% Set constants dependent on m
d1 = 2;
d2 = 4*m;
A2 = 2*m*(m+2);
A3 = -m/2;

% Creating polynomial equations
% z1 is odd
% z3, z5 are even
syms t;
a = sym('a',[1,N]);
b = sym('b',[1,N]);
c = sym('c',[1,N]);
z1 = sum (a.*(t).^[1:N])
isa (z1)
z2 = diff(z1);
z3 = sum (b.*(t).^[1:N])+b0

% Convert to scalar from vector or find another polynomial generator
% Matlab confused about array vs. scalar
z4 = diff(z3);
z5 = sum (c.*(t).^[1:N]);
z6 = diff(z5);
e1 = diff(z2);
e2 = diff(z4);
e3 = diff(z6);

% Collecting coefficients
% Setting formula bounds
E1 = collect(-2*(d1-1)*z3^4*d1*z2^2 - 2*d2*z2*z4*z3^3*z1*d1 + 2*z2*z6*z1*z3^4*d1 + 2*(d1-1)*z3^4*d1 + 2*A3*z1^4 + e*z1^2*z3^4*d1 - 2*e1*z1*z3^4*d1, t)
E2 = -2*d1*z2*z4*z3^3*d2 - 2*(d2 - 1)*z4^2*z1*z3^2*d2 + 2*z4*z6*z1*z3^3*d2 + 2*A2*z1*z3^2 - 4*A3*z1^3 + e*z3^4*z1*d2 - 2*e2*z1*z3^3*d2
E3 = -z6*(d1*z2*z3+d2*z4*z1)+z6^2+c+e*z5-e3*z1*z3

% Gathering coefficients
c1 = coeffs(E1, t)
c2 = coeffs(E2, t)
c3 = coeffs(E3, t)

% Setting initial conditions
initcond1 = strcat(char(subs(z1,t,0)),'=0');
initcond2 = strcat(char(subs(z2,t,0)),'=1');
initcond3 = strcat(char(subs(z3,t,0)),'=a');
initcond4 = strcat(char(subs(z4,t,0)),'=0');
initcond5 = strcat(char(subs(z5,t,0)),'=0');
initcond6 = strcat(char(subs(z6,t,0)),'=0');
initcond7 = strcat(char(subs(e3,t,0)),'=-1/3');

% Solving coefficients
% [a5, b5, c5] = solve(initcond1, initcond2, initcond3, initcond4, initcond5, initcond6, initcond7, e1, e2, e3)

end

