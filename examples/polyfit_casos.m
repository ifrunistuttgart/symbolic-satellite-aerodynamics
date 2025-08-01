% This is an example script showcasing how CaSOS can be used for polynomial fitting


%% 1D example: Generate nonlinear function
x = linspace(-1, 11, 100);
y = 3 * sin(.2.*x) + 0.03 * randn(size(x)); % Nonlinear function with noise



xP = casos.Indeterminates('x', 1);
P = casos.PD(monomials(xP, [0:4]).to_vector);
c = casos.PS.sym('c', size(P));
PFun = to_function(P);
Pval = PFun(x);

% Define as sos problem
sos.f = (c'*Pval-y)*(c'*Pval-y)';
sos.x = c;
S = casos.sossol('S', 'mosek', sos);
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


%% Solve using conic interface
[X1, X2] = meshgrid(linspace(-1,1,50),linspace(-1,1,50));
Y = 3 * sin(.2.*X1) + cos(1.*X2);

% User defines polynomial for fitting
xP = casos.Indeterminates('x', 2);
P = casos.PS.sym('c', monomials(xP, [0:1]));

% Inside a function: extract casadi symbolic variables
PFun = to_function(casos.PD(P.monomials.to_vector));
c = poly2basis(P)


% Evaluate function
Pval = PFun(X1(:)', X2(:)');


sos.x = c;
sos.f = (c'*Pval-Y(:)')*((c'*Pval-Y(:)'))';

S = casos.sossol('S', 'mosek', sos);
sol = S();

YP = reshape(full(sol.x'*Pval), size(Y));
figure(2);
clf;
hold on;
surface(X1,X2,Y, EdgeAlpha=0, MeshStyle="both");
surface(X1,X2,YP, EdgeAlpha=0);
axis equal
hold off;


%% Call fit function
x = casos.PD('x', 2);
p = casos.PS.sym('c',monomials(x, [0:2]));
[S, sol] = cfit.fit([X1(:), X2(:)], Y(:), p);