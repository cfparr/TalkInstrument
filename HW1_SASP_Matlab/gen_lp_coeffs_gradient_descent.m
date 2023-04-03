function a = gen_lp_coeffs_gradient_descent(x,y, M)
    % To efficiently compute the LPC coefficients, we use the Levinson-Durbin recursion 
    % (fast algorithm exploiting the Toeplitz structure of the auto-correlation matrix)
    % returns a_0, a_1, ... a_M for a signal x
    % length of the source signal x
N_x = length(x);
N_y = length(y);
%% Wiener-Hopf solution
% Compute the autocorrelation matrix

% Compute an estimate of the autocorrelation of the input signal
[r_ac, r_lags] = xcorr(x);
    
% Build auto-correlation matrix
r_ac = r_ac(r_lags >= 0); % Take positive lags 0,...,p
r_ac = 1/N_x * r_ac(1:M); % normalizing with the length of the source signal 

R = toeplitz(r_ac); % Create Toeplitz matrix from vector
%A Toeplitz matrix is a diagonal-constant matrix, which means all elements along a diagonal have the same value

% Estimation of the cross-correlation between the input signal and the desired response y
[p_cc, p_lags] = xcorr(x,y);
p_cc = flipud(p_cc);
% Build cross-correlation matrix
p_cc = p_cc(p_lags >= 0); % Take positive lags 0,...,p

% compute the cross-correlation vector
p = 1/N_y * p_cc(1:M); % normalizing with the length of the desired signal

% compute the Wiener-Hopf solution
w_o = R\p; 
  %% Steepest Descent
% Determine the step-size parameter 
num_steps = 2000;   %Number of steps of the iterative algorithm
eigenvalues = eig(R);   %Eigenvalues of the autocorrelation matrix
factor = 0.95;               

mu = factor * (2/max(eigenvalues));  %Update step

% determine the global time constant of the Steepest Descent algorithm
tau = 1 / (2 * mu * min(eigenvalues));

% initialize the vector of filter taps
w = zeros(M,1); %Initial guess of the filter coefficients: we don't have any prior knowledge so we set it to null vector

for n = 1:num_steps
% perform the iterative update
    w = w + mu*(p - R*w); 
end
    
    a = [1; w];
end
