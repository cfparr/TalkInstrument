function [lp_coeffs] = gen_lp_coeffs(x, M)
% Efficient computation of LPC coefficients using Levinson-Durbin recursion
% Returns a_0, a_1, ..., a_M for a signal x

    % Compute autocorrelation coefficients using built-in function
    rx = xcorr(x, 'biased');
    
    % Extract autocorrelation coefficients up to order M
    rx = rx(numel(x):(numel(x)+M));
    
    % Generate Toeplitz matrix using toeplitz() function
    toeplitz_matrix = toeplitz(rx(1:end-1));
    
    % Solve the system of linear equations using backslash operator
    lp_coeffs = [1; -toeplitz_matrix \ rx(2:end)];
end
