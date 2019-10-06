clear, clc
AEM = [9 0 2 6];

%% Design Requirements ( Chebysev BE )
a_max = 0.4 + AEM(4) / 36;  % db
a_min = 27 + AEM(3) * 5/9;
LF_Gain_Req_DB = 10;

f_0 = 2.5 * 1e3;           % Hz
f_1 = 1700 + 50 * AEM(3);
f_2 = f_0^2 / f_1;
D = ( f_0^2 - f_1^2 ) / ( 2.1 * f_1 );
f_3 = 0.5 * ( -D + sqrt( D^2 + 4 * f_0^2 ) );
f_4 = f_0^2 / f_3;

omega_1 = 2*pi * f_1;       % rad/sec
omega_2 = 2*pi * f_2;
omega_3 = 2*pi * f_3;
omega_4 = 2*pi * f_4;

bw = omega_2 - omega_1;
omega_0 = sqrt( omega_1 * omega_2 );
q_c = omega_0 / bw;


%% Prototype ( LP Chebysev ) Filter Parameters
Omega_p = 1;
Omega_s = bw / ( omega_4 - omega_3 );

%   - filter degree
n = acosh(...
        sqrt( (10 ^ (0.1 * a_min) - 1) / (10 ^ (0.1 * a_max) - 1) ) ...
    ) / acosh(Omega_s);
n = ceil(n);

%   - epsilon param ( for Inverse Chebysev )
epsilon = sqrt(10 ^ (0.1 * a_max) - 1);

%   - half power frequency
Omega_hp = cosh( (1 / n) * acosh( 1 / epsilon ) );
assert( Omega_hp > 1 );


%% Poles / Zeros of Prototype Chebysev LPF
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


%% Frequency Transformation LP -> HP
%   - transform poles
for k = 1 : n_pairs
   poles_ch(k) = poles_ch(k).inverse;
end

%   - add zeros
Omega_z = zeros( n_pairs, 1 );


%% Frequency Transformation HP -> BE ( Geffe Algorithm )
%   - transform poles
be_poles = Pole.empty( 2 * n_pairs, 0 );
bei = 1;
for k = 1 : n_pairs
    
    [be_poles(bei), be_poles(bei + 1)] = ...
        Pole.geffe( poles_ch(k), omega_0, bw );
    bei = bei + 2;
    
    % PLUS TWO ZEROS @ omega = omega_0

end

%   - transform zeros
be_zeros = Inf( 2 * n_pairs, 1 );
bzi = 1;
for zi = 1 : n_pairs
    
    if ( Omega_z( zi ) < Inf )
        
        [be_zeros(bzi), be_zeros(bzi + 1)] = ...
            Pole.geffez( Omega_z( zi ), omega_0, bw );
        bzi = bzi + 2;

        % PLUS TWO POLES @ omega = 0
        
    end
    
end

% Pole-Zero Plot
[pole_plus, pole_minus] = be_poles(1).sigmaOmega;
P = [pole_plus, pole_minus];
for k = 2 : length( be_poles )
    [pole_plus, pole_minus] = be_poles(k).sigmaOmega;
    P = cat( 2, P, [pole_plus, pole_minus] );
end

Z = cat( 1, 1i * be_zeros, -1i * be_zeros );
pzplot( zpk( Z, P', 1 ) )


%% Utilize Units
%   - init units holder
n_units = 2 * n_pairs;
units = FilterUnit( n_units, 0 );

%   - units' parameters
for k = 1 : n_units
    
    units( k ) = FilterUnit( ...
        be_poles( k ).Omega0, ...
        be_poles( k ).Q, ...
        be_zeros( k ) ...
    );

    if ( be_poles( k ).Omega0 < be_zeros( k ) )
        
        units( k ) = fried_lpn( units( k ) );
        
    else
        
        units( k ) = fried_hpn( units( k ) );
        
    end

end

%% Combine sub-units
LF_Gain = 1;
A = units(1).TF;

% Plot tf of each sub-unit
for k = 1 : n_units
   
%     plot_transfer_function( ...
%         units(k).TF, ...
%         ( 0.5 / pi ) * [omega_1, omega_2, omega_3, omega_4, omega_0] ...
%     );
% 
%     set(gcf, 'name', ['Unit #' num2str(k) ' | ' units(1, k).name], ...
%         'numbertitle','off' );
    
    % Calculate Gain at omega_0
    LF_Gain = LF_Gain * units(k).k_lf;

    if ( k > 1 )
        
        A = series( A, units(k).TF );
        
    end
    
end

% Compensate gain ( gain @ omega_0 should be 10dB )
A = ( 10^( LF_Gain_Req_DB / 20 ) ) * ( 1 / LF_Gain ) * A;

% Plot Amplitude
plot_transfer_function( A, ( 0.5 / pi ) * [omega_1, omega_2, omega_3, omega_4, omega_0] );
set(gcf, 'name', 'Total Response | Amplitude', 'numbertitle','off' );

% % Plot Attenuation
% a = inv(A);
% plot_transfer_function( a, ( 0.5 / pi ) * [omega_1, omega_2, omega_3, omega_4, omega_0] );
% set(gcf, 'name', 'Total Response | Attenuation', 'numbertitle','off' );


%% Test Resulting System
t = 0 : 1/40000 : 1/50 - 1/40000;
input = 0.8 * cos( ( omega_0 - ( omega_0 - omega_3 ) / 2 ) * t ) + ...
    cos( ( omega_0 + ( omega_0 + omega_3 ) / 2 ) * t ) + ...
    cos( 0.5 * omega_1 * t ) + 0.8 * cos( 2.4 * omega_2 * t ) + ...
    0.4 * cos( 3.5 * omega_2 * t );
test_sys( A, 'custom', t, input, 40000 );














