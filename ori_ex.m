% Licensed under the terms of the GPL v3. See LICENSE for details.
clear; close all; clc

% Where the data is:
datafolder = '../data/';
list = dir(strcat(datafolder,'trajectory*'));


%%% Read in data to matrix
M = length(list);                       % number of trajectories
fobj = fopen(strcat(datafolder, list(1).name));
trajInput = textscan(fobj, '%f%f', 'HeaderLines',1, 'Delimiter','\n', 'CollectOutput',1);
fclose(fobj);
traj = cell2mat(trajInput);
N = length(traj);

trajectories = zeros(M,N);

for i = 1:length(list)
    fname = strcat(datafolder, list(i).name);
    fobj = fopen(fname);
    trajInput = textscan(fobj, '%f%f', 'HeaderLines',1, 'Delimiter','\n', 'CollectOutput',1);
    fclose(fobj);
    traj = cell2mat(trajInput);
    trajectories(i,:) = traj(:,2)';     % add M trajectories to matrix
     % fprintf( 'Read in file #%d, "%s"\n', i, fname);
end
time = traj(:,1);                       % our x-axis, (assume same in all files)



%%%% OPTIONAL: Find index of first time to include in fitting procedure
starttime = 200;                    % first time point to include in fitting
epsilon0 = time(2) - time(1);           % time step size
start_idx = int16(1 + starttime / epsilon0);

trajectories = trajectories(:,start_idx:end);
time = time(start_idx:end);
fprintf('# Removed the first %d sampling times, new N = %d\n', start_idx, N);
N = length(time);



%%% OPTIONAL: Reduce N
everyN = 1;
time = time(1:everyN:N);
trajectories = trajectories(:, 1:everyN:N);
N = length(time);



%%%% OPTIONAL: Reduce M (if we want to work with a smaller set of trajectories)
new_M = 0;                              % If 0, don't use this
if new_M > 0
    fprintf('# Found M = %d trajectories\n', M);
    start_traj = 0;          % can be changed to use a _different_ set of the full M
    stop_traj = start_traj + new_M;
    assert(stop_traj <= M);
    trajectories = trajectories(1+start_traj:stop_traj, :);
    M = new_M;
end
%%%% END

guess = [1,1];                          % initial parameter guess

% Perform the actual fit, (using fminunc, with default options)
[params, sigma, chi2_min] = wlsice(time, trajectories, guess);

% RESULT:
fprintf('Optimal param: %d\n', params);
fprintf('Sigma: %d\n', sigma);
fprintf('Chi^2 value: %d\n', chi2_min);


%%% Also get goodness-of-fit parameters
% note, by y_mean here we mean the average over N, not M!

y = mean(trajectories);

y_mean = sum(y) / N;

SS_tot = sum(square(y_mean - y));
SS_reg = sum(square(y_mean - f(times, params)));
SS_res = sum(square(y - f(times, params)));
coeff = (1- SS_res / SS_tot)
