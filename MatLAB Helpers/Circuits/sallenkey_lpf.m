function sallenkeyLPF = sallenkey_lpf(unit)
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
    unit.k_lf = k;
    unit.k_hf = 0;
    
    
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
    % TODO
    
    unit.TF = unit.k_lf * tf( num, denom );  
    
    % return new unit
    sallenkeyLPF = unit;
end
