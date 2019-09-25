clear, clc
AEM = [9 0 2 6];

%% Design Requirements ( Inverse Chebysev BP )
a_max = 0.4 + AEM(4) / 36;  % db
a_min = 34 + AEM(3) * 5/9;
Omega_0_Gain_Req_DB = 10;

f_0 = 1.65 * 1e3;           % Hz
f_1 = 1400 + 25 * AEM(4);
f_2 = f_0^2 / f_1;
D = 2.5 * ( f_0^2 - f_1^2 ) / f_1;
f_3 = 0.5 * ( -D + sqrt( D^2 + 4 * f_0^2 ) );
f_4 = f_0^2 / f_3;

omega_1 = 2*pi * f_1;       % rad/sec
omega_2 = 2*pi * f_2;
omega_3 = 2*pi * f_3;
omega_4 = 2*pi * f_4;

bw = omega_2 - omega_1;
omega_0 = sqrt( omega_1 * omega_2 );
q_c = omega_0 / bw;


%% Prototype ( LP ) Filter Parameters
Omega_p = 1;
Omega_s = ( omega_4 - omega_3 ) / bw;

% Inverse Chebysev requires Omega_s = 1
Omega_p = Omega_p / Omega_s;
Omega_s = 1;

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


%% Poles / Zeros of Prototype Inverse Chebysev LPF
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

% Inverse Chebysev requires Omega_s = 1 but LP -> BP requires Omega_p = 1
% Therefore, we require Omega_p' = 1
%   - scale poles
for k = 1 : n_pairs
   poles_ich(k) = poles_ich(k).scaleOmega0( 1 / Omega_p );
end
%   - scale zeros
for zi = 1 : n_pairs
    Omega_z( zi ) = Omega_z( zi ) * ( 1 / Omega_p );
end
%   - scale frequencies
Omega_p = 1;
Omega_s = Omega_s / Omega_p;


%% Frequency Transformation ( Geffe Algorithm )
%   - transform poles
bp_poles = Pole.empty( 2 * n_pairs, 0 );
bpi = 1;
for k = 1 : n_pairs
    
    [bp_poles(bpi), bp_poles(bpi + 1)] = ...
        Pole.geffe( poles_ich(k), omega_0, bw );
    bpi = bpi + 2;
    
    % PLUS TWO ZEROS @ omega = 0

end

%   - transform zeros
bp_zeros = Inf( 2 * n_pairs, 1 );
bzi = 1;
for zi = 1 : n_pairs
    
    if ( Omega_z( zi ) < Inf )
        
        [bp_zeros(bzi), bp_zeros(bzi + 1)] = ...
            Pole.geffez( Omega_z( zi ), omega_0, bw );
        bzi = bzi + 2;

        % PLUS TWO POLES @ omega = 0
        
    end
    
end


%% Utilize Units
%   - init units holder
n_units = 2 * n_pairs;
units = FilterUnit( n_units, 0 );

%   - units' parameters
for k = 1 : n_units
    
    units( k ) = FilterUnit( ...
        bp_poles( k ).Omega0, ...
        bp_poles( k ).Q, ...
        bp_zeros( k ) ...
    );

    if ( bp_poles( k ).Omega0 < bp_zeros( k ) )
        
        units( k ) = fried_lpn( units( k ) );
%         units( k ) = boctor_lpn( units( k ) );
        
    else
        
        units( k ) = fried_hpn( units( k ) );
        
    end

end

%% Combine sub-units
Omega_0_Gain = 1;
A = units(1).TF;

% Plot tf of each sub-unit
for k = 1 : n_units
   
    plot_transfer_function( ...
        units(k).TF, ...
        ( 0.5 / pi ) * [omega_1, omega_2, omega_3, omega_4, omega_0] ...
    );

    set(gcf, 'name', ['Unit #' num2str(k) ' | ' units(1, k).name], ...
        'numbertitle','off' );
    
    % Calculate Gain at omega_0
    Omega_0_Gain = Omega_0_Gain * squeeze( bode(units(k).TF, omega_0) );

    if ( k > 1 )
        
        A = series( A, units(k).TF );
        
    end
    
end

% Compensate gain ( gain @ omega_0 should be 10dB )
A = ( 10^( Omega_0_Gain_Req_DB / 20 ) ) * ( 1 / Omega_0_Gain ) * A;

% Plot Amplitude
plot_transfer_function( A, ( 0.5 / pi ) * [omega_1, omega_2, omega_3, omega_4, omega_0] );
set(gcf, 'name', 'Total Response | Amplitude', 'numbertitle','off' );
% 
% % Plot Attenuation
% a = inv(A);
% plot_transfer_function( a, ( 0.5 / pi ) * [omega_1, omega_2, omega_3, omega_4, omega_0] );
% set(gcf, 'name', 'Total Response | Attenuation', 'numbertitle','off' );


% %% Test Resulting System
% t = 0 : 1/40000 : 1/50 - 1/40000;
% input = cos( ( omega_0 - ( omega_0 - omega_1 ) / 3 ) * t ) + ...
%     0.6 * cos( ( omega_0 + ( omega_0 + omega_1 ) / 4 ) * t ) + ...
%     0.7 * cos( 0.5 * omega_3 * t ) + 0.8 * cos( 2.4 * omega_4 * t ) + ...
%     0.6 * cos( 3 * omega_4 * t );
% test_sys( A, 'custom', t, input, 40000 );















