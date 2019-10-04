classdef Pole
    %POLE Summary of this class goes here
    %   Detailed explanation goes here
    
    
    properties
        % p = -sigma + j*omega
        sigma
        omega
        
        % omega_0, Q
        omega_0 = NaN
        q = NaN
        
        is_real
    end
    
    
    methods
        
        function obj = Pole(sigma, omega)
            obj.sigma = sigma;
            obj.omega = omega;
            
            obj.is_real = ( omega == 0 );
            
            omega_0 = sqrt( sumsqr( [sigma; omega] ) );
            obj.omega_0 = omega_0;
            obj.q = omega_0 / ( 2 * sigma );
        end
        
        function omega_0 = Omega0(obj)
            if isnan( obj.omega_0 )
                omega_0 = sqrt( sumsqr( [obj.sigma; obj.omega] ) );
            else
                omega_0 = obj.omega_0;
            end
        end
        function q = Q(obj)
            if isnan( obj.omega_0 )
                q = obj.Omega0 / ( 2 * obj.sigma );
            else
                q = obj.q;
            end
        end
        
        function obj = scaleOmega0(obj, Omega_k)
            obj = Pole.fromOmega0AndQ( obj.Omega0 * Omega_k, obj.Q );
        end
        
        function inversePole = inverse(obj)
            inversePole = Pole.fromOmega0AndQ( 1 / obj.Omega0, obj.Q );
        end
        
        function [pole_plus, pole_minus] = sigmaOmega(obj)
            pole_plus = -obj.sigma + 1i * obj.omega;
            pole_minus = -obj.sigma - 1i * obj.omega;
        end
        
    end
    
    
    methods(Static)
        
        function pole = fromOmega0AndQ( omega_0, Q )
            % Denominator of TS: S^2 + ( OMEGA_0 / Q ) * S + OMEGA_0^2
            re = omega_0 / ( 2 * Q );
            im = sqrt( omega_0 ^ 2 - re ^ 2 );
            
            pole = Pole( re, im );
        end
        
        function pole = guillemin(alpha, angle)
            pole = Pole( ...
                sinh( alpha ) * cosd( angle ), ...
                cosh( alpha ) * sind( angle ) ...
            ); 
        end
        
        function [bpf_pole1, bpf_pole2] = geffe( lpf_pole, omega_0, bw )
            q_c = omega_0 / bw;
            
            Sigma = lpf_pole.sigma;
            Omega = lpf_pole.omega;
            
%             fprintf( "Pole: -%f +- j%f", Sigma, Omega );

            C = sumsqr( [Sigma, Omega] );
            D = 2 * Sigma / q_c;
            E = 4 + C / ( q_c^2 );
            G = sqrt( E^2 - 4 * D^2 );
            Q = ( 1/D ) * sqrt( 0.5 * ( E + G ) );
            M = ( Sigma * Q ) / q_c;
            W = M + sqrt( M^2 - 1 );

            omega_0_1 = W * omega_0;
            omega_0_2 = ( 1/W ) * omega_0;

            % first of the two new pole-pairs
            sigma = omega_0_1 / ( 2 * Q );
            omega = omega_0_1 * sin( acos( 1 / ( 2 * Q ) ) );
            bpf_pole1 = Pole( sigma, omega );

            % last of the two new pole-pairs
            sigma = omega_0_2 / ( 2 * Q );
            omega = omega_0_2 * sin( acos( 1 / ( 2 * Q ) ) );
            bpf_pole2 = Pole( sigma, omega );

            % PLUS TWO ZEROS @ origin
        end
        
        function [bpf_zero1, bpf_zero2] = geffez( lpf_zero, omega_0, bw )
            q_c = omega_0 / bw;
            
%             fprintf( "Zero: +- j%f", lpf_zero );
            
            K = 2 + ( lpf_zero / q_c )^2;
            x = 0.5 * ( K + sqrt( K^2 - 4 ) );

            assert( x > 1 );

            bpf_zero1 = omega_0 * sqrt(x);
            bpf_zero2 = omega_0 / sqrt(x);
            
            % PLUS TWO POLES @ origin
        end
        
    end
    
end

