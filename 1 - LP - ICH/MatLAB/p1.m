clear plot_transfer_function_nofig
clc
AEM = [9 0 2 6];

%% Design Requirements ( Inverse Chebysev BE )
a_max = 0.6125;         % db
a_min = 22.75;
LF_Gain_Req_DB = 0;

f_p = 4.4 * 1e3;        % Hz
f_s = 1.9 * f_p;

omega_p = 2*pi * f_p;   % rad/sec
omega_s = 2*pi * f_s;

% normalize acc. to omega_s
Omega_p = omega_p / omega_s;
Omega_s = 1;


%% Initial Calculations
%   - filter degree
n = acosh(...
        sqrt( (10 ^ (0.1 * a_min) - 1) / (10 ^ (0.1 * a_max) - 1) ) ...
    ) / acosh(1 / Omega_p);
n = ceil(n);

%   - epsilon param ( for Inverse Chebysev )
epsilon = 1 / sqrt(10 ^ (0.1 * a_min) - 1);

%   - half power frequency
Omega_hp = 1 / cosh( (1 / n) * acosh(1 / epsilon) );
assert( Omega_hp < 1 );
omega_hp = Omega_hp * omega_s;


%% Poles / Zeros of Inverse Chebysev LPF
%   - number of pole-pairs
n_pairs = floor( n / 2 ) + mod(n, 2);

%   - find Butterworth angles
psi = zeros(n_pairs, 1);
if ( mod( n, 2 ) == 0 )
    psi( 1 ) = 90 / n;
end
for k = 2 : n_pairs
    psi( k ) = psi( k - 1 ) + 180 / n;
end

%   - find Chebysev poles ( Guillemin Algorithm )
alpha = (1 / n) * asinh(1 / epsilon);
poles_ch = Pole.empty( n_pairs, 0 );
for k = 1 : n_pairs
   poles_ch(k) = Pole.guillemin( alpha, psi( k ) );
end

%   - find Inverse Chebysev poles
poles_ich = Pole.empty( n_pairs, 0 );
for k = 1 : n_pairs
   poles_ich(k) = poles_ch( k ).inverse;
end

%   - find Inverse Chebysev zeros
Omega_z = Inf( n_pairs, 1 );
zi = 1;
for k = 1 : 2 : ( n - 1 )
    Omega_z( zi ) = sec( ( k * pi ) / ( 2 * n ) );
    zi = zi + 1;
end


%% Utilize Units
%   - init units holder
n_units = n_pairs;
units = FilterUnit( n_units, 0 );   % max degree utilized is 2 ( biquads )

%   - units' parameters
for k = 1 : n_units
    
    units( k ) = FilterUnit( ...
        omega_s * poles_ich( k ).Omega0, ...
        poles_ich( k ).Q(), ...
        omega_s * Omega_z( k ) ...
    );
    
    if ( ~units( k ).is_bilinear )
        
        units( k ) = boctor_lpn( units( k ) );
%         units( k ) = fried_lpn( units( k ) );
        
    else
        
        % TODO: passive bilinear realization
        
    end

end


%% Combine sub-units
LF_Gain = 1;
A = units(1).TF;

% Plot tf of each sub-unit
for k = 1 : n_units
   
    plot_transfer_function( ...
        units(k).TF, ...
        ( 0.5 / pi ) * [omega_s, omega_p, omega_hp] ...
    );

    set(gcf, 'name', ['Unit #' num2str(k) ' | ' units(1, k).name], ...
        'numbertitle','off' );
    
    % Calculate Gain at HF
    LF_Gain = LF_Gain * units(k).k_lf;

    if ( k > 1 )
        
        A = series( A, units(k).TF );
        
    end
    
end

% Compensate gain ( gain @ DC should be $LF_Gain_Req_DB dB )
A = ( 10^( LF_Gain_Req_DB / 20 ) ) * ( 1 / LF_Gain ) * A;

% Plot Amplitude
plot_transfer_function( A, ( 0.5 / pi ) * [omega_s, omega_p, omega_hp] );
set(gcf, 'name', 'Total Response | Amplitude', 'numbertitle','off' );

% % Plot Attenuation
% a = inv(A);
% plot_transfer_function( a, ( 0.5 / pi ) * [omega_s, omega_p, omega_hp] );
% set(gcf, 'name', 'Total Response | Attenuation', 'numbertitle','off' );


%% Test Resulting System
test_sys( A, 'sawtooth', 2e3 );

