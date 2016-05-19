function out = tfilt(in,window)
% This function filters the in vector with a triangle filter of
% (odd) window size window.

    sz = size(in);
    if (sz(1) == 1)
        in = in';
    end
    
    out = conv2(in,triang(window)/sum(triang(window)),'same');
    out(1:((window-1)/2)) = NaN;
    out(end-((window-1)/2-1):end) = NaN;
    if (sz(1) == 1)
        out = out';
    end
end
