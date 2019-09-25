function sallenkeyHPF = sallenkey_hpf(unit)
%SALLENKEY_HPF Summary of this function goes here
%   Detailed explanation goes here
    
    k = 3 - 1 / unit.Q;

    
    %% Circuit Elements
    %   - resistors
    unit.R = ones( 2, 1 );
    unit.r = ones( 2, 1 );
    unit.r(2) = k - 1;
    
    %   - capacitors
    unit.C = ones( 2, 1 );
    
    %   - gains
    unit.k_hf = k;
    unit.k_lf = unit.k_hf * unit.Omega_z ^ 2;
    
    
    %% Scaling
    % requirements
    C_1n = 1e-6;    % 1.0uF
    
    % scaling coefficients
    k_f = unit.omega_0;
    k_m = unit.C(1) / ( C_1n * k_f );
    
    % scale resistors
    unit.R = unit.R .* k_m;
    unit.r = unit.r .* k_m;
    unit.C = unit.C ./ ( k_m * k_f );
    
    %% Transfer function
    num = zeros( 1, 3 );
    denom = zeros( 1, 3 );
    
    num(1) = 1;
    
    denom(1) = 1;
    denom(2) = 1 / ( unit.R(1) * unit.C(1) ) + 1 / ( unit.R(1) * unit.C(2) ) ...
        + ( 1 - k ) / ( unit.R(2) * unit.C(1) );
    denom(3) = 1 / ( unit.R(1) * unit.R(2) * unit.C(1) * unit.C(2) );
    
    unit.TF = unit.k_hf * tf( num, denom );  
    
    % return new unit
    sallenkeyHPF = unit;
end
