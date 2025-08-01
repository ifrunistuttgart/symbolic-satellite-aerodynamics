function [Problem, sol] = fit(x,y,p)
% Inputs:
%  - x: kxn with n: number of inputs variables and k: number of
%  measurements
%  - y: kx1 Measurement vector
%  - p: casos.PS symbolic polynomial. Number of variables needs to match k
arguments
    x (:,:) {mustBeNumeric} % allow n inputs
    y (:,1) {mustBeNumeric} % only fit functions with scalar output
    p (1,1) {isa(p, 'casos.PS')} % only fit scalar output polynomials
end

nP  = size(x,1); % number of evaluation points
nin = size(x,2); % number of input variables

% Check dimension of y
if ~(nP==size(y,1))
    error('fit:DimensionMismatch: Number of columns in x (%d) must',...
        'match number of elements in y (%d).', nP, size(y,2));
end

% Check Polynomial dimension
nIndets = numel(p.indeterminates);
if ~(nIndets==nin)
    error('fit:DimensionMismatch: Number of indets in p (%d) must',...
        'match number of rows in x (%d).', nIndets, nin);
end


% Grab polynomial coefficients
c = poly2basis(p); % 1xc_n vector

% Create vectorial monomial function
% f:(x_1,...,x_n)->(poly[number of monomials]) SXFunction
pFun = to_function(casos.PD(p.monomials.to_vector));

% Evaluate polynomial function without coefficients
xCell = num2cell(x',2);
Pval = pFun(xCell{:});


% Define fitting problem using casos
sos.x = c;

% Minimizer
sos.f = (c'*Pval-y(:)')*((c'*Pval-y(:)'))';

Problem = casos.sossol('S', 'mosek', sos);
sol = Problem();
end