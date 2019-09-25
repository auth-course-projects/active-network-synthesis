function friedHPN = fried_hpn(unit)
%FRIED_LPN Summary of this function goes here
%   Detailed explanation goes here

    k_1 = ( 1 / unit.Omega_z )^2 - 1;
    k_2 = ( (2 + k_1) * unit.Q^2 ) / ( (2 + k_1) * unit.Q^2 + 1 );
    k = k_2 * ( 1 / unit.Omega_z )^2;

    
    %% Circuit Elements
    %   - resistors
    unit.R = zeros( 4, 1 );
    unit.R(1) = 1;
    unit.R(2) = unit.Q^2 * ( k_1 + 2 )^2;
    unit.R(3) = 1;
    unit.R(4) = unit.Q^2 * ( k_1 + 2 );
    
    %   - capacitors
    unit.C = zeros( 3, 1 );
    C = 1 / ( unit.Q * ( 2 + k_1 ) );
    unit.C(1) = k_1 * C;
    unit.C(2) = C;
    unit.C(3) = C;
    
    %   - gains
    % HF_gain is given in eq. 7.137 ( Prof. Theocharis' notes - Chapter 7 )
    unit.k_hf = k;
    % LF_gain = HF_gain * ( 1 / k_1 + k_1 * unit.Omega_z ^ 2 - 1 )
    % IF k_1 == ( omega_0 / omega_z ) ^ 2 OR k_1 == 1, 
    % THEN LF_GAIN = HF_gain * ( omega_z / omega_0 ) ^ 2
    % k_1 is assumed 1 for DC gain calculation
    unit.k_lf = unit.k_hf * unit.Omega_z ^ 2;
    
    
    %% Scaling
    % requirements
    C_1n = 1e-6;    % 1.0uF
    
    % scaling coefficients
    k_f = unit.omega_0;
    k_m = C / ( C_1n * k_f );
    
    % scale resistors
    unit.R = unit.R .* k_m;
    unit.C = unit.C ./ ( k_m * k_f );
    
    
    %% Transfer function
    num = zeros( 1, 3 );
    denom = zeros( 1, 3 );
    
    num(1) = 1;
    num(2) = ( k_2 - 1 ) / ( k_2 * ( k_1 + 1 ) * unit.R(1) * unit.C(2) ) + ...
        ( k_1 + 2 ) / ( ( k_1 + 1 ) * unit.R(2) * unit.C(2) );
    num(3) = 1 / ( ( k_1 + 1 ) * unit.R(1) * unit.R(2) * unit.C(2)^2 );
    
    denom(1) = 1;
    denom(2) = ( k_1 + 2 ) / ( unit.R(2) * unit.C(2) );
    denom(3) = 1 / ( unit.R(1) * unit.R(2) * unit.C(2) ^ 2 );
    
    unit.TF = unit.k_hf * tf( num, denom );  
    
    % return new unit
    friedHPN = unit;
end

