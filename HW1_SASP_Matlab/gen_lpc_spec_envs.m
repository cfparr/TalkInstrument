function spec_envs = gen_lpc_spec_envs(windowed_signal, M, nfft)
% This function returns a matrix of spectral envelopes, where column m
% is spectral envelope for m'th signal frame.
% Parameters of the function:
% windowed_input: In this matrix, each column is a windowed signal.
% M: Defines the order of linear predictor.
% nfft: Defines the fft size.

num_frames = size(windowed_signal, 2);
spec_envs = zeros(nfft, num_frames);

for m = 1:num_frames
    xm = windowed_signal(:, m); % get mth column
    coeffs =gen_lp_coeffs(xm, M); %gen_lp_coeffs_gradient_descent(x,y,M);% % compute LP coefficients
    spec_env = 1 ./ abs(fft(coeffs, nfft)); % compute spectral envelope
    spec_envs(:, m) = spec_env; % assign to m-th column of spec_envs
end

end
