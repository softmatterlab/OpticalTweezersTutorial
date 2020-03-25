function values = df(t, params)
% The gradient of the analytical function that is to be fitted to data,
% through parameters in params, as function of t. Here we use a powerlaw

df_d_lamda = zeros(length(params), length(t));
df_d_lamda(1,:) = power(t, params(2));
df_d_lamda(2,:) = params(1)*log(t) .* power(t, params(2));
values = df_d_lamda;
