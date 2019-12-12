function [matches] = lf_select_points(matches,mask)

    matches.matches = matches.matches(mask,:);
    matches.rays1 = matches.rays1(mask);
    matches.rays2 = matches.rays2(mask);
    matches.disparity1 = matches.disparity1(mask);
    matches.disparity2 = matches.disparity2(mask);
    matches.M1s = matches.M1s(mask);
    matches.M2s = matches.M2s(mask);
    matches.npoints = sum(mask);
end