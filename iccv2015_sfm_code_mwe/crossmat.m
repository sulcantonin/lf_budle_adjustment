%% Compute matrix corresponding to a vector cross product
%
% (c) Bastian Goldluecke 11/2014, University of Konstanz
% License: Creative Commons BY-SA 4.0
%          Please cite our paper if you use our code in your work.

function M = crossmat( v )

    M = zeros(3,3);
    M(1,2) = -v(3);
    M(1,3) = v(2);
    M(2,3) = -v(1);
    M(3,2) = v(1);
    M(3,1) = -v(2);
    M(2,1) = v(3);
    
end
