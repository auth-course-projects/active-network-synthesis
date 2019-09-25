classdef FilterUnit
    %UNIT Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        name
        
        % design parameters
        Q
        omega_0
        omega_z
        
        % normalized
        Omega_0        
        Omega_z
        
        % flags
        has_zero    % non Inf
        is_bilinear
        
        % circuit elements
        R
        r   % for Sallen-Key circuits
        C
        
        % gains
        k_hf
        k_lf
        
        % when utilized
        TF
    end
    
    methods
        function obj = FilterUnit(omega_0, Q, omega_z)
            obj.omega_0 = omega_0;
            obj.Q = Q;
            
            if ( nargin > 2 )
                obj.omega_z = omega_z;
                obj.has_zero = true;
            else
                obj.omega_z = Inf;
                obj.has_zero = false;
            end
           
            % normalize
            obj.Omega_0 = 1.0;
            if ( nargin > 2 )
                obj.Omega_z = omega_z / omega_0;
            else
                obj.Omega_z = Inf;
            end
            
            % check Q
            obj.is_bilinear = ( Q == 0.5 );
            
            % set name
            if ( obj.is_bilinear )
                obj.name = 'Bilinear Passive';
            else
                if ( obj.omega_z ~= Inf )
                    if ( omega_z > omega_0 )
                        obj.name = 'Biquad LPN';
                    else
                        obj.name = 'Biquad HPN';                        
                    end
                else
                    obj.name = 'Biquad LP';
                end
            end
        end
    end
end

