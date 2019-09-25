function test_sys( sys, type, param1, param2, param3 )
    
    if ( nargin == 1 )
        F_s = 2e3;
        type = 'sine';
    elseif ( nargin < 3 && ~strcmp( type, 'custom' ) )
        F_s = 2e3;
    elseif ( ~strcmp( type, 'custom' ) )
        F_s = param1;
    else
        F_s = param3;
    end
    
    %% Generate Input Signal
    switch ( type )
       
        case 'sine'
            [input, t] = gensig( type, 1/F_s, 5/F_s, 1/( 50 * F_s ) );
            F_s = 50 * F_s;
            
        case 'square'
            [input, t] = gensig( type, 1/F_s, 5/F_s, 1/( 50 * F_s ) );
            F_s = 50 * F_s;
            
        case 'sawtooth'
            t = 0 : 1/( 50 * F_s ) : 5/F_s;
            input = 0.5 + 0.5 * sawtooth( 2 * pi * F_s * t );   % between [0,1]
            F_s = 50 * F_s;
            
        case 'custom'
            t = param1;
            input = param2;
        
    end
    
    %% Simulate Output Signal
    output_sine = lsim( sys, input, t );
    
    %% Plots
    figure
    subplot(2, 1, 1)
    plot( t, input );

    subplot(2, 1, 2)
    plot_spectrum( input, F_s );

    figure
    
    subplot(2, 1, 1)
    plot( t, output_sine );

    subplot(2, 1, 2)
    plot_spectrum( output_sine, F_s );

end

