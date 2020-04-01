function values = d2f(t, params)
% The hessian of the analytical function that is to be fitted to data,
% through parameters in params, as function of t. Here we use a powerlaw

a = params(1);
b = params(2);

d2fTemp = zeros(length(params), length(params), length(t));

% d2f / da da
% df2(1,1) = zeroes(size(t))

% d2f / da db
d2fTemp(1,2,:) = log(t) .* power(t, b);

% d2f / db da
d2fTemp(2,1,:) = d2fTemp(1,2);

%d2f / db db
d2fTemp(2,2,:) = a*(log(t).^2) .* power(t, b);

values = d2fTemp;
