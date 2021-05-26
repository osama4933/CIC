function [pnts] = get_bounded_max(arr,up_thresh,low_thresh)
%GET_BOUNDED_MAX, this function returns all points in arr that are >
%low_thresh and < up_thresh

        pnts_up = find(arr < up_thresh);
        pnts_low = find(arr > low_thresh);
        pnts = intersect(pnts_up,pnts_low);
end

