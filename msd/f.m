function value = f(t, params)
% The analytical function that is to be fitted to data, through
% parameters in params, as function of t. Here we use a powerlaw

%value = params(1)*power(t, params(2));
value = params(1)*abs(1.0-exp(-t/params(2)));
