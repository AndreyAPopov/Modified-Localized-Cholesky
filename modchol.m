%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% (c) 2020 Andrey A Popov
%
% Licensed under the MIT license
%
% This function takes the UNSCALED anomalies
%
%    A = X - mu
%
% a distance function and a radius of influence
% and outputs T, Di, such that
%
%    B^-1 = T.'*Di*T
%
% sparse, such that vector products with B^-1
% can be done in an efficient manner.
%
% If the factorization has to be done multiple
% times in a row, with the same distance function
% then the sparsity pattern can  be pre-computed
% and returned, then subsequently supplied to a
% new call of the function
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [T, Di, sparcity] = modchol(A, dist, radius, sparcity)

if nargin < 3
    error('This function requires 3-4 input arguments');
end

P = @(i, js) dist(i, js) <= radius;

[n, N] = size(A);

if nargin < 4
    sparcity = calcsparcity(P, n);
end

nel = nnz(sparcity) + n;

Tis = zeros(1, nel);
Tjs = zeros(1, nel);
Tvs = zeros(1, nel);

Tis(1:n) = 1:n;
Tjs(1:n) = 1:n;
Tvs(1:n) = ones(1, n);

ind = n + 1;

E = zeros(n, N);
E(1, :) = A(1, :).';
for i = 2:n
    js = find(sparcity(:, i));
    Zi = A(js, :);
    bi = (Zi*Zi.')\(Zi*(A(i, :).'));
    
    E(i, :) = A(i, :).' - Zi.'*bi;
    
    njs = numel(js);
    
    Tis(ind:(ind + njs - 1)) = i*ones(1, njs);
    Tjs(ind:(ind + njs - 1)) = js;
    Tvs(ind:(ind + njs - 1)) = -bi.';
        
    ind = ind + njs;
end

T = sparse(Tis, Tjs, Tvs, n, n);

d = sum(E.*E, 2)/(N - 1);

Di = spdiags(1./d, 0, n, n);

end

function sparcity = calcsparcity(P, n)

% I don't think that there is an easy way to pre-allocate this.
% Fortunately, this is a function that should only be ran when 
% either the distances between the states change, or the number
% of states changes.

is = [];
js = [];
vs = [];

for i = 2:n
    ne = find(P(i, 1:(i - 1)));
    is = [is, ne];
    js = [js, i*ones(1, numel(ne))];
    vs = [vs, ones(1, numel(ne))];
end

sparcity = sparse(is, js, vs, n, n);

end

