function [out] = get_max(arr,threshold,num_pnts)
%GET_MAX, this function returns uptil num_pnts many maximums above
% 'threshold'

    out = [];
    for i = 1:num_pnts
        [a b] = max(arr);
        if(a < threshold)
            return;
        else
            out = [out b];
            arr(b) = 0;
        end
    end

end

