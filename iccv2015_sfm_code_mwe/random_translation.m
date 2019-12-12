%% Create a random translation vector up to a maximum length
%
% (c) Bastian Goldluecke 11/2014, University of Konstanz
% License: Creative Commons BY-SA 4.0
%          Please cite our paper if you use our code in your work.
%
function [t] = random_translation( max_length )

    t = rand( 3,1 );
    n = norm(t);
    if n>0
        t = t / n;
    end
    
    t = (max_length*rand()) * t;
end
