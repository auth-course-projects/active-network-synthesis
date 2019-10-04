function boctorLPN = boctor_lpn(unit)
%BOCTOR_LPN Summary of this function goes here
%   Detailed explanation goes here

    %   - define k1
    k_1 = ( ( unit.Omega_0 / unit.Omega_z ) ^ 2 + 1 ) / 2;
    
    
    %% Circuit Elements
    %   - resistors
    unit.R = zeros( 6, 1 );
    unit.R(1) = 2 / ( k_1 * unit.Omega_z ^ 2 - 1 );
    unit.R(2) = 1 / ( 1 - k_1 );
    unit.R(3) = k_1 / ( 2 * unit.Q ^ 2 ) + 1 / unit.R(1);
    unit.R(4) = 1 / k_1;
    unit.R(5) = 1;
    unit.R(6) = 1;
    
    %   - capacitors
    unit.C = zeros( 2, 1 );
    unit.C(1) = k_1 / ( 2 * unit.Q );
    unit.C(2) = 2 * unit.Q;
    
    %   - gains
    % HF_gain is given in eq. 7.158 ( Prof. Theocharis' notes - Chapter 7 )
    unit.k_hf = 1 / ( 1 + unit.R(3) );
    % LF_gain = HF_gain * ( 1 / k_1 + k_1 * unit.Omega_z ^ 2 - 1 )
    % IF k_1 == ( omega_0 / omega_z ) ^ 2 OR k_1 == 1, 
    % THEN LF_GAIN = HF_gain * ( omega_z / omega_0 ) ^ 2
    % k_1 is assumed 1 for DC gain calculation
    unit.k_lf = unit.k_hf * unit.Omega_z ^ 2;
    
    
    %% Scaling
    % requirements
    C_1n = 1e-8;    % 0.01uF
    
    % scaling coefficients
    k_f = unit.omega_0;
    k_m = unit.C(1) / ( C_1n * k_f );
    
    % scale resistors
    unit.R = unit.R .* k_m;
    unit.C = unit.C ./ ( k_m * k_f );
    
    
    %% Transfer function
    num = zeros( 1, 3 );
    denom = zeros( 1, 3 );
    
    num(1) = 1;
    R_24 = ( unit.R(2) * unit.R(4) ) / ( unit.R(2) + unit.R(4) );
    num(2) = 0;     % there's a mistake in TF given in eq. 7-156
    num(3) = ( unit.R(1) + R_24 + unit.R(6) ) / ...
        ( unit.R(1) * R_24 * unit.R(6) * unit.C(1) * unit.C(2) );
    
    denom(1) = 1;
    denom(2) = 1 / ( unit.R(6) * unit.C(2) ) + 1 / ( R_24 * unit.C(2) );
    denom(3) = 1 / ( unit.R(4) * unit.R(6) * unit.C(1) * unit.C(2) );
    
    unit.TF = unit.k_hf * tf( num, denom );    

    % return new unit
    boctorLPN = unit;
end

