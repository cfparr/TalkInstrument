function signal = get_istft(stft, R)
    % Performs Overlap-Add reconstruction of original signal
    nfft = size(stft, 1); % size of the FFT
    num_frames = size(stft, 2); % number of FFT windows
    signal = zeros(R * (num_frames - 1) + nfft, 1);
    for m = 1:num_frames
        idx = (m-1)*R + 1;
        windowed_signal = real(ifft(stft(:, m)));
        signal(idx:idx+nfft-1) = signal(idx:idx+nfft-1) + windowed_signal;
    end
    % Normalize by the number of FFT points
    signal = signal ./ nfft;
end
