clear, clc
AEM = [9 0 2 6];

%% Design Requirements
a_max = 0.54 + AEM(4) / 36; % db
a_min = 24 + AEM(3) * 6/9;
HF_Gain_Req_DB = 5;

f_p = 4 * 1e3;              % Hz
f_s = f_p / 2.6;

omega_s = 2*pi * f_s;       % rad/sec
omega_p = 2*pi * f_p;


%% Prototype ( LP ) Filter Parameters
Omega_p = 1;
Omega_s = omega_p / omega_s;

%   - filter degree
n = log10( ( 10^(0.1 * a_min) - 1 ) / ( 10^(0.1 * a_max) - 1 ) );
n = n / ( 2 * log10( Omega_s / Omega_p ) );
n = ceil(n);

%   - half power frequency
Omega_0 = Omega_p / ( 10^(0.1 * a_max) - 1 ) ^ ( 1 / ( 2 * n ) );
assert( Omega_0 > 1 );
omega_0 = omega_p / Omega_0;


%% Poles / Zeros of Prototype Butterworth LPF
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

%   - find Butterworth poles
poles = Pole.empty( n_pairs, 0 );
for k = 1 : n_pairs
   poles(k) = Pole.fromOmega0AndQ( 1, 1 / ( cosd( psi( k ) ) * 2 ) );
end


%% Frequency Transformation ( LP -> HP )
%   - transform poles
%     Since Butterworth HPF is being designed, no pole changes will occur.
%     Frequency transformation simply adds $n_pairs zeros at origin.

%   - scale poles
hp_poles = Pole.empty( n_pairs, 0 );
for k = 1 : n_pairs
   hp_poles(k) = poles(k).scaleOmega0( omega_0 );
end

%   - add zeros
hp_zeros = zeros( n_pairs, 1 );

% Pole-Zero Plot
[pole_plus, pole_minus] = hp_poles(1).sigmaOmega;
P = [pole_plus, pole_minus];
for k = 2 : length( hp_poles )
    [pole_plus, pole_minus] = hp_poles(k).sigmaOmega;
    P = cat( 2, P, [pole_plus, pole_minus] );
end

Z = cat( 1, 1i * hp_zeros, -1i * hp_zeros );
pzplot( zpk( Z, P', 1 ) )


%% Utilize Units
%   - init units holder
n_units = n_pairs;
units = FilterUnit( n_units, 0 );

%   - units' parameters
for k = 1 : n_units
    
    units( k ) = FilterUnit( ...
        hp_poles( k ).Omega0, ...
        hp_poles( k ).Q, ...
        hp_zeros( k ) ...
    );

    units( k ) = sallenkey_hpf( units( k ) );

end


%% Combine sub-units
HF_Gain = 1;
A = units(1).TF;

% Plot tf of each sub-unit
for k = 1 : n_units
   
    plot_transfer_function( ...
        units(k).TF, ...
        ( 0.5 / pi ) * [omega_s, omega_p, omega_0] ...
    );

    set(gcf, 'name', ['Unit #' num2str(k) ' | ' units(1, k).name], ...
        'numbertitle','off' );
    
    % Calculate Gain at HF
    HF_Gain = HF_Gain * units(k).k_hf;

    if ( k > 1 )
        
        A = series( A, units(k).TF );
        
    end
    
end

% Compensate gain ( gain @ omega_0 should be 10dB )
A = ( 10^( HF_Gain_Req_DB / 20 ) ) * ( 1 / HF_Gain ) * A;

% Plot Amplitude
plot_transfer_function( A, ( 0.5 / pi ) * [omega_s, omega_p, omega_0] );
set(gcf, 'name', 'Total Response | Amplitude', 'numbertitle','off' );

% Plot Attenuation
a = inv(A);
plot_transfer_function( a, ( 0.5 / pi ) * [omega_s, omega_p, omega_0] );
set(gcf, 'name', 'Total Response | Attenuation', 'numbertitle','off' );


%% Test Resulting System
Fs = 5e4;
t = 0 : 1/Fs : 1 - 1/Fs;
input = cos( 0.5 * omega_s * t ) + 0.6 * cos( 0.8 * omega_s * t ) + ...
    cos( 1.2 * omega_p * t ) + 0.8 * cos( 2.4 * omega_p * t ) + ...
    0.4 * cos( 3.5 * omega_p * t );
test_sys( A, 'custom', t, input, Fs );













