% This is an example script showcasing how CaSOS can be used for polynomial fitting


%% 1D example: Generate nonlinear function
x = linspace(-1, 11, 100);
y = 3 * sin(.2.*x) + 0.03 * randn(size(x)); % Nonlinear function with noise



xP = casos.PD('x', 1);
P = casos.PD(monomials(xP, [0:3]).to_vector);
c = casos.PS.sym('c', size(P));
PFun = to_function(P);
Pval = PFun(x);

% Minimizer
sos.f = (c'*Pval-y)*(c'*Pval-y)';
sos.x = c
S = casos.sossol('S', 'mosek', sos)
sol = S();

% Back to function
fSol = to_function(sol.x'*P);
fSolVal = full(fSol(x))


fig = figure(1);
clf;
hold on;
plot(x,y);
plot(x,fSolVal);
xlabel('$x$', Interpreter='Latex');
ylabel('$f(x)$', Interpreter='Latex');
hold off;

