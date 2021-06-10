function [windows] = active_sess_split(windows, max_window_length, window_overlap)
    % Function to split large active session windows in a binary fashion.
    % Each active session that is larger than the max_window_length will be
    % continually divided until it meets this cutoff. Split up sessions
    % will be overlapping (# of overlap samples = window_overlap) the
    % adjacent chunks in order to detect the packets that are sliced in the
    % middle
    % Note: max_window_length is expected to be larger than window_overlap,
    % undefined bahavior will occur otherwise.

    i = 0;
    while i < length(windows)
        i = i + 1;
        if(abs(windows(i, 2) - windows(i, 1)) > max_window_length)  % need to split window in half
            len = abs(windows(i, 2) - windows(i, 1)) / 2;               % get bisected length
            win1_end = ceil(windows(i, 1) + len + (window_overlap / 2));        % add overlap to each half
            win2_start = floor(windows(i, 2) - len - (window_overlap / 2));
            
            windows = [windows(1:i - 1,:); [windows(i,1), win1_end]; windows(i:end,:)]; % bisect the window
            windows(i+1, 1) = win2_start;
            i = i - 1;  % see if the first bisected half is within the max window length
        end
    end
end