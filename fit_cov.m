function [params, chi2_min, C, sigma] = fit_cov(tau, msd_exp, guess)
% Perform the correlated corrected weighted least squares fit on input data.
% M is number of tranjectories, N number of sampling points,

[M, N] = size(msd_exp');
assert(N == length(tau));

%%% Make covariance matrix
msd_mean = mean(msd_exp, 2);
size(msd_mean)
Y = zeros(N,M);
for m = 1:M

    Y(:,m) = msd_exp(:,m)' - msd_mean';
end
C = Y * Y' / (M - 1.0);

C=C(2:end, 2:end);
C = C / M;
R = diag(1./diag(C), 0);
%R=1e-19*inv(C);
%chi2 = @(p, t, y, R) (y - fexp(t, p)')' * R * (y - fexp(t, p)');

% Do the parameter fitting
[params, chi2_min] = fminunc(@(par)(chiqu(par, tau(2:end), msd_mean(2:end), R)), guess)

gradient = df(tau(2:end), params);
fhessian = d2f(tau(2:end), params);

q = length(params);
first_term = zeros(q, q);
delta = f(tau(2:end), params) - msd_mean(2:end);

for a = 1:q
    for b = 1:q
        for i = 1:N-1
            for j = 1:N-1
                first_term(a,b) = first_term(a,b) + 2 * fhessian(a,b,i) * R(i,j) * delta(j);
            end
        end
    end
end

second_term = 2 * gradient * R * gradient';
hessian = first_term + second_term;

H_inv = inv(hessian);
RCR = R * C * R;

error = zeros(q, q);
for a = 1:q
    for b = 1:q
        % Should we not be able to move it to up here? Like so:
        % dfRCRdf = gradient * RCR * gradient';
        for c = 1:q
            for d = 1:q
                dfRCRdf = gradient(c,:) * RCR * gradient(d,:)';
                error(a,b) = error(a,b) + 4 * H_inv(a,c) * dfRCRdf * H_inv(d,b);
            end
        end
    end
end

sigma = sqrt(diag(error));


