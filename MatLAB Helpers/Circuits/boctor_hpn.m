function boctorHPN = boctor_hpn(unit)
%BOCTOR_LPN Summary of this function goes here
%   Detailed explanation goes here

    % Check Fundamental Condition
    if( unit.Q >= ( 1 / ( 1 - unit.Omega_z ^ 2 ) ) )
        
        boctorHPN = fried_hpn( unit );
        return
        
    end

    %   - define gamma, alpha
    gamma = 1 / ( unit.Q * ( 1 - unit.Omega_z ^ 2 ) );
    alpha = ( gamma^2  - 1 ) / ( 1 + gamma^2 * unit.Omega_z ^ 2 );
    H = 2;
    
    %% Circuit Elements
    %   - resistors
    unit.R = zeros( 6, 1 );
    unit.R(1) = ( ( 1 + alpha ) / ( alpha * gamma^2 * unit.Omega_z ) ) ^ 2;
    unit.R(2) = 1;
    unit.R(3) = ( H * (1+alpha)^2 ) / ( alpha * gamma^2 );
    unit.R(4) = ( 1 / alpha ) * ( H / ( H - 1 ) );
    unit.R(5) = H / alpha;
    unit.R(6) = unit.R(3) / ( H * ( 1 + alpha - alpha * unit.Omega_z^2 ) - 1 );
    
    %   - capacitors
    unit.C = ( ( alpha * gamma ) / ( 1 + alpha ) ) * ones( 2, 1 );
    
    %   - gains
    % HF_gain is DEFINED
    unit.k_hf = H;
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
    
    R_eq13 = (unit.R(1)*unit.R(3)/(unit.R(1)+unit.R(3)));
    R_eq1 = (R_eq13*unit.R(6)/(R_eq13+unit.R(6)));
    R_eq2 = unit.R(2) + (unit.R(4)*unit.R(5)/(unit.R(4)+unit.R(5)));
    
    num(1) = 1;
    num(2) = 0;
    num(3) = 1 / ( unit.R(1) * unit.R(2) * unit.C(1) * unit.C(2) );
    
    denom(1) = 1;
    denom(2) = ( 1 / ( R_eq1 * unit.C(1) ) ) * ( 1 - (R_eq1 * R_eq2) / ( unit.R(1) * unit.R(2) ) );
    denom(3) = 1 / ( R_eq1 * R_eq2 * unit.C(1) * unit.C(2) );
    
    unit.TF = unit.k_hf * tf( num, denom );    

    % return new unit
    boctorHPN = unit;
end

