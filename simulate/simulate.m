% Set simulation parameters as desired.
n_subjects = 1024; % number of subjects to simulate
n_time = 1024; % number of timepoints to simulate per subject
rho = 0.5; % correlation between brain and motion signals
transform_motion_time = @(x)(1 + x.*x); % multiplicative motion transform in time domain, default @(x)(1)
transform_motion_sub = @(x)(x); % transform of per-subject motion mixing proportion, default @(x)(x)
transform_brain_sub = @(x)(x); % same transforms for brain signal
motion_factor = 1; % ratio of motion to brain signal, default 1 = same intensity
smooth_n = 1; % integer for how much to temporally smooth motion data, default 1 = no smoothing
proper_mixing = true; % whether to use the variance correction strategy for proper timeseries mixing.

%% Perform Simulation

% Load brain images and motion images.
img1 = loadpng('brain.png'); % brain
img2 = loadpng('motion.png'); % motion
assert(isequal(size(img1), size(img2)));

% Infer dimensions from image files.
n_nodes = length(img1); % number of nodes
n_edges = (n_nodes*n_nodes - n_nodes) / 2; % number of unique edges

% Compute nearest positive semidefinite correlation matrices.
% (Higham's alternating projection algorithm.)
corrmat1 = nearcorr(img1);
corrmat2 = nearcorr(img2);
clear img1 img2

% Compute cholesky decompositions.
chol1 = chol(corrmat1);
chol2 = chol(corrmat2);

% Randomly generate per-subject proportions of brain and motion signal.
% 1st column img1 (brain), 2nd column img2 (motion)
X_full = normrnd(0, 1, n_subjects, 2);

% Pre-apply per-subject transforms.
X_t = [transform_brain_sub(X_full(:,1)), ...
      transform_motion_sub(motion_factor .* X_full(:,2))];

% Shift columns to be non-zero.
minX = abs(min([min(X_full); zeros(size(X_full,2))]));
X_full = X_full + minX;
minX_t = abs(min([min(X_t); zeros(size(X_t,2))]));
X_t = X_t + minX_t;

% Pre-compute smoothing kernel.
smooth_kernel = ones(smooth_n,1);
smooth_kernel = smooth_kernel ./ sum(smooth_kernel);

% Pre-allocate memory for simulated timeseries data.
data = zeros(n_time, n_nodes, n_subjects); % fMRI BOLD timeseries
data_nomotion = zeros(n_time, n_nodes, n_subjects); % without motion
motion = zeros(n_time, n_subjects); % motion timeseires, analogous to FD

% Generate motion-only simulated data first.
for j_sub=1:n_subjects
    % Generate random motion timeseries from standard normal distribution.
    data_j_sub = normrnd(0, 1, n_time, n_nodes);
    
    % Add desired correlation structure using Cholesky factorization.
    % Generate univariate FD timeseries from simulated motion data without
    % accounting for per-subject mixing proportions or offset noise.
    motion(:, j_sub) = var(data_j_sub .* sqrt(X_full(j_sub,2)) * chol2, 0, 2);
    
    % Apply desired "confounding" timeseries transforms.
    % The per-subject transform adds variance to all time points.
    % The time transform adds variance as a function of the univariate
    % timeseries above.
    % The direct transform operates on the timeseries data directly.
    data_j_sub = data_j_sub .* ... % direct transform
        sqrt(X_t(j_sub,2)) .* ... % per-subject transform
        transform_motion_time(motion(:,j_sub)); % time-dependent transform
    
    % Temporally smooth the motion data.
    % Restore original variance after smoothing.
    var_factor = sqrt(sum(data_j_sub.^2)./n_time);
    for j_node=1:n_nodes
        data_j_sub(:,j_node) = conv( ...
            data_j_sub(:,j_node), ...
            smooth_kernel, 'same');
    end
    sumsq = sum(data(:,:,j_sub).^2);
    if sumsq > 0
        data(:,:,j_sub) = data(:,:,j_sub) ...
            .* var_factor ...
            ./ sqrt(sumsq./n_time);
    end
    
    % Update the no-motion data.
    data_nomotion(:,:,j_sub) = data_j_sub;
    
    % Use Cholesky factorization to add desired correlation structure.
    data_j_sub = data_j_sub * chol2;
    
    % Update full simulation data.
    data(:,:,j_sub) = data_j_sub;
        
    % Generate true-but-unobservable FD timeseries after transform and
    % compute true mixing proportion from it (needed for simulation only).
    X_t(j_sub,2) = mean(var(data_j_sub, 0, 2));
end
clear var_factor sumsq

% Unshift columns.
X_full = X_full - minX;
X_t = X_t - minX_t;
clear minX

% Introduce desired correlation between brain and motion proportions.
X_full(:,1) = gencorr(X_t(:,2), rho, X_full(:,1));
X_t(:,1) = gencorr(X_t(:,2), rho, X_t(:,1));

% Rescale to 1 and then apply motion factor.
X_t(:,1) = X_t(:,1)./std(X_t(:,1)).*std(X_t(:,2));
X_full(:,1) = X_full(:,1)./std(X_full(:,1)).*std(X_full(:,2));
X_t(:,1) = X_t(:,1) ./ motion_factor;
X_full(:,1) = X_full(:,1) ./ motion_factor;

% Shift all columns to be non-zero.
minX_t = abs(min([min(X_t); zeros(size(X_t,2))]));
X_t = X_t + minX_t;

% Find global variance scale term so that we can faithfully simulate
% correlation, not just covariance.
max_var = max(sum(X_t,2));

% Complete simulation by adding brain signal.
for j_sub=1:n_subjects
    % Generate random brain timeseries.
    brain_j_sub = normrnd(0, sqrt(X_t(j_sub,1)), n_time, n_nodes) * chol1;
    
    % Add to full simulated data.
    data_j_sub = data(:,:,j_sub) + brain_j_sub;
    
    % Now we must add some additional random noise to "correct" the
    % variance so that our simulation satisfies the assumptions for doing
    % regression on correlation matrices (as opposed to covariance
    % matrices, which would be much more straightforward).
    if proper_mixing
        data_j_sub = data_j_sub + correct_var_(data_j_sub, max_var);
    end
    
    % Add in a small amount of offset noise so that no subject can
    % accidentally have zero timeseries signal.
    offset_j_sub = normrnd(0, 1, n_time, n_nodes);
    data(:,:,j_sub) = data_j_sub + offset_j_sub;
    
    % Do same for no-motion simulated data.
    data_j_sub = data_nomotion(:,:,j_sub) + brain_j_sub;
    if proper_mixing
        data_j_sub = data_j_sub + correct_var_(data_j_sub, max_var);
    end
    data_nomotion(:,:,j_sub) = data_j_sub + offset_j_sub;
end
clear data_j_sub brain_j_sub offset_j_sub

% As a negative control, uncomment to simulate some brain-only data
% without any motion in it.
% for j_sub=1:n_subjects
%     data_nomotion(:,:,j_sub) = ...
%         normrnd(0, sqrt(X_t(j_sub,1)), n_time, n_nodes) * chol1 + ...
%         normrnd(0, sqrt(max_var - X_t(j_sub,1)), n_time, n_nodes) + ...
%         normrnd(0, 1, n_time, n_nodes);
% end

% Unshift brain columns (not strictly necessary).
X_t = X_t - minX_t;
clear minX_t

%% Write Simulated Data to Disk

for i = 1:n_subjects
    sub_i = struct;
    sub_i.data = data(:,:,i);
    sub_i.fd = motion(:,i);
    sub_i.trait = X_full(i,1);
    save(['sub' num2str(i, '%04d') '.mat'], '-struct', 'sub_i');
end
clear sub_i

%% Variance Correction

function var_correction_ = correct_var_(data_, max_var_)
    % Generate some random noise that, when added to data, "corrects" the
    % variance so that our simulation satisfies the assumptions for doing
    % regression on correlation matrices (as opposed to covariance
    % matrices, which would be much more straightforward).
    % 
    % We need to add noise such that:
    % 1. The difference in edgewise variance at each timepoint is minimized
    % (i.e. the total signal variance is constant over time).
    % 2. The average variance for this subject is equal to the "maximum"
    % variance of the highest-variance subject (i.e. the total signal
    % variance is constant across subjects).
    
    % First figure out for how many timepoints the edgewise variance
    % temporarily exceeds the maximum variance.
    data_var_ = var(data_, 0, 2);
    excess_var_ = find(data_var_ > max_var_);
    
    % Find a target variance for the remaining timepoints such such that
    % the average variance across time will equal the maximum variance.
    target_var_ = ...
        ... % total variance at maximum
        ( max_var_ * size(data_, 1) - ...
        ... % total excessive variance
        sum(data_var_(excess_var_) - max_var_) ) / ...
        ... % number of timepoints
        size(data_, 1);
    
    % Generate a random variance correction timeseries.
    var_correction_ = normrnd(0, 1, size(data_,1), size(data_,2));
    
    % So that we don't have to keep track of nested covariance, make each
    % row of the variance correction timeseries orthogonal to the
    % corresponding row in the timeseries data.
    if any(data_var_)
        for j_=1:size(data_,1)
            var_correction_(j_,:) = orthogonalize(var_correction_(j_,:)',data_(j_,:)')';
        end
    end
    
    % Scale each row so that it adds the needed amount of variance to reach
    % the target variance.  Zero out rows that already have excess
    % variance.
    var_correction_ = ...
        var_correction_ ./ std(var_correction_, 0, 2) .* ...
        sqrt(max([zeros(size(data_,1),1), target_var_ - data_var_], [], 2));
end