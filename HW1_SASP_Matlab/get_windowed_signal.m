function windowed_signal = get_windowed_signal(signal, window_size, hop_size, window_type)
    % Generate a windowed signal by dividing x into overlapping frames
    % and applying a window function to each frame.
    %
    % Parameters:
    % signal: input signal (vector or matrix)
    % window_size: length of window function
    % hop_size: hop size (distance between adjacent frames)
    % window_type: window function type ('rect', 'hann', or 'hamming')
    %
    % Returns:
    % windowed_signal: matrix of overlapping frames, where each column is
    % a windowed signal frame.
    
    % handle input arguments
    if nargin < 4
        window_type = 'hann';
    end
    
    % generate window function
    switch window_type
        case 'rect'
            window_fn = ones(window_size, 1);
        case 'hann'
            window_fn = hann(window_size);
        case 'hamming'
            window_fn = hamming(window_size);
        otherwise
            error('Invalid window type');
    end
    
    % partition signal into windows and apply window function
    num_windows = ceil((size(signal, 1) - window_size) / hop_size) + 1;
    windowed_signal = zeros(window_size, num_windows);
    for n = 1:num_windows
        start_idx = (n-1) * hop_size + 1;
        end_idx = start_idx + window_size - 1;
        if end_idx > size(signal, 1)
            % zero-pad the last window if necessary
            windowed_signal(:, n) = [signal(start_idx:end, :); zeros(end_idx - size(signal, 1), size(signal, 2))];
        else
            windowed_signal(:, n) = signal(start_idx:end_idx, :);
        end
        windowed_signal(:, n) = windowed_signal(:, n) .* window_fn;
    end
end
