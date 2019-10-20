function value = fexp(t, params)
% The analytical function that is to be fitted to data, through
% parameters in params, as function of t. Here we use a powerlaw

value = params(1)*abs(1.0-exp(-t/params(2)));
%modelfun=@(b,x)(b(1))*(1-exp(-x/b(2)));