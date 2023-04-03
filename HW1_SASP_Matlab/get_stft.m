function ffts = get_stft(windowed_signal, nfft)
%GET_STFT computes the STFT of a windowed signal
%   windowed_signal: windowed signal matrix
%   nfft: the size of the FFT (optional, defaults to the size of the window)
%
%   Returns a 2-D matrix where each column is a single FFT frame

    if nargin < 2
        nfft = size(windowed_signal, 1);
    end
    ms = size(windowed_signal, 2);
    ffts = zeros(nfft, ms);
    for m = 1:ms
        xm = windowed_signal(:, m);
        freq_window = fft(xm, nfft);
        ffts(:, m) = freq_window;
    end
end
