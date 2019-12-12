%% Compute errors for rotation and translation estimates
%% See example.m for a verbose description of returned errors.
%
% (c) Bastian Goldluecke 11/2014, University of Konstanz
% License: Creative Commons BY-SA 4.0
%          Please cite our paper if you use our code in your work.
%
function errors = lf_project_absolute( R, t, Ra, ta )

    cos_angle = 0.5 * ( trace( Ra * R' ) - 1.0 );
    if cos_angle < -1.01 || cos_angle > 1.01
        fprintf( 'cosa violation R\n' );
        cat( 2,R,Ra )
    end
    cos_angle = max( -1.0, min( 1.0, cos_angle ));
    R_angular = 180.0 / pi * acos( cos_angle );
    Rdiff = Ra - R;
    R_absolute = norm( Rdiff(:) );

    N = max( 0.0001, norm(Ra(:)) );
    R_relative = R_absolute / N;
    
    cos_angle2 = dot( t, ta ) / (norm( t ) * norm( ta ));
    if cos_angle2 < -1.01 || cos_angle2 > 1.01
        fprintf( 'cosa violation t\n' );
        cat( 2,t,ta )
    end
    cos_angle2 = max( -1.0, min( 1.0, cos_angle2 ));
    t_angular = 180.0 / pi * acos( cos_angle2 );
    t_absolute = norm ( t-ta );
    
    N = max( 0.0001, norm(ta) );
    t_relative = t_absolute / N;

    errors = [ R_angular R_absolute R_relative t_angular t_absolute t_relative ];
    
end
