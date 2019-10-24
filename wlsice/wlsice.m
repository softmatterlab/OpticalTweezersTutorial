function [params, sigma, chi2_min, C] = wlsice(time, trajectories, guess)
% Perform the correlated corrected weighted least squares fit on input data.
% M is number of tranjectories, N number of sampling points,
%trajectories=trajectories';

[M, N] = size(trajectories);
assert(N == length(time));


%%% Make covariance matrix
y_mean = mean(trajectories)';

Y = zeros(N,M);
for m = 1:M
    Y(:,m) = trajectories(m,:) - y_mean';
end
C = Y * Y' / (M - 1.0);
C=C(2:end, 2:end);
y_sigma = sqrt(diag(C));                 %  In case we want it


%%% WLS-ICE parameter
C = C / M;

%R = inv(diag(diag(C), 0));             % works in octave 4.0.2
R = diag(1./diag(C), 0);              % works in matlab and octave

% chi2-square function to be minimized:
chi2 = @(p, t, y, R_) (y - f_acf_lf(t, p)')' * R_ * (y - f_acf_lf(t, p)');

% Do the parameter fitting
[params, chi2_min] = fminunc(@(par)(chi2(par, time(2:end), y_mean(2:end), R)), guess);


%%% WLS-ICE sigma
% WLS-ICE error estimation, valid also for non-linear fitting. R could be R =
% inv(cov), or the diagonal of that, or some other symmetric matrix of our
% choosing. With N sampling points, and k parameters we have:
%
% df     is a  k x N      dimensional array
% d2f    is a  k x k x N  dimensional array
% delta  is a  N          dimensional array

gradient = df(time(2:end), params);
fhessian = d2f(time(2:end), params);

q = length(params);
first_term = zeros(q, q);
delta = f_acf_lf(time(2:end), params) - y_mean(2:end);


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
