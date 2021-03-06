# Modified-Localized-Cholesky
 
Usage

```matlab
% given X which is (n x N) sampled from a distribution N(mu, B)
% a vectorized distance function 'dist' s.t. dist(i, js) is a 
% row vector of distances from each j to i, and a radius of 
% influence, calculate T, Di s. t. T.'*Di*T is an approximation 
% to B^-1.

% calculate the anomalies
muE = mean(X, 2);
A = X - muE;

[T, Di, sparcity] = modchol(A, dist, radius);

% modify X somehow
X = modifyX(X, T, Di);

% recalculate the anomalies
muE = mean(X, 2);
A = X - muE;

% re-use the sparcity developed in the previous function call
% to speed up the next calculation
[T, Di] = modchol(A, dist, radius, sparcity);
```
