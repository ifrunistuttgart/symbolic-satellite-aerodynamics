classdef Sentman < saero.aerodynamics.CalculationModel
    %SENTMAN implements Calculation model with Sentman calculation method
    %   Detailed explanation goes here
    
    properties
        parameters
    end
    
    properties (Access=private)
        cosdelta = sym('cosd', 'real');
        significant_digits = 5
    end
    
    methods
        function obj = Sentman(parameters)
            %SENTMAN Construct an instance of this class
            %   Detailed explanation goes here
            arguments
                parameters.alpha_E = 0.95; % Energy accommodation coeff
                parameters.si = 5; % Molecular speed ratio
                parameters.Tw = 300; % Wall temperature
                parameters.Vi = 7800; % scalar velocity
                parameters.rho = 1e-10; % density
            end

            obj.parameters = parameters;
            
            % Add constant values
            obj.parameters.kB = 1.380649e-23;  % Boltzmann Constant
            obj.parameters.mT = 16 * 1.6605390689252e-27;

            % Define symvar used internally
            % obj.cosdelta = sym('cosd', 'real');
        end
        
        % Calculate the force vector for each panel:
        function force_vectors = calculate_force(obj, wind_direction, ...
            satellite_geometry)

            % Satellite geometry
            N = satellite_geometry.get_normals;
            A = satellite_geometry.get_panel_areas;
            
            % Number of panels
            n = size(N, 2);
            
            % Incoming wind velocity in B frame (vectorial)
            vi_B = repmat(wind_direction, 1, n);
            
            % Helper terms
            H1 = obj.H1();
            H2 = obj.H2();
            
            % Substitute H1 and H2 for each panel
            H1vec = sym(zeros(1,n));
            H2vec = H1vec;

            for i=1:n
                H1vec(i) = subs(H1, obj.cosdelta, -N(:,i)'*vi_B(:,1));
                H2vec(i) = subs(H2, obj.cosdelta, -N(:,i)'*vi_B(:,1));
            end

            % Stacking and significant digits
            H1vec = vpa(repmat(H1vec,3,1),obj.significant_digits);
            H2vec = vpa(repmat(H2vec,3,1),obj.significant_digits);
            
            p = obj.parameters;
            
            % momentum flux and force
            p_B = p.rho/2.*p.Vi.^2.*(H1vec.*N + H2vec.*vi_B);
            force_vectors = real(repmat(A,3,1).*p_B); % Repmat required for Casadi Implementation
        end
    end

    methods (Access=private) % Private helper functions for Sentman
        function H1 = H1(obj)
            % Shorthand for parameters
            p = obj.parameters;

            % Temperature ratio (1xn):
            enum = p.si.*obj.cosdelta.*erfc(-p.si.*obj.cosdelta);
            denom = 1./sqrt(pi).*exp(-p.si.^2*obj.cosdelta.^2) + ...
                p.si.*obj.cosdelta.*erfc(-p.si.*obj.cosdelta);
            Trat = p.alpha_E.*(2.*p.kB.*p.Tw)./(p.mT.*p.Vi.^2).*p.si^2 + ...
                (1-p.alpha_E).*(...
                    1 + p.si.^2./2 + 1/4.* enum./denom...
                );

            % Helper function
            Gamma1 = @(x) (x ./ sqrt(pi)) .* exp(-x.^2) + (0.5 + x.^2) .* (1 + erf(x));
            Gamma2 = @(x) (1 ./ sqrt(pi)) .* exp(-x.^2) + x .* (1 + erf(x));


            % H1 as 1x1 symexpr
            H1 = 1./(p.si.^2).*(-Gamma1(p.si.*obj.cosdelta) ...
                - sqrt(pi)/2.*sqrt(Trat).*Gamma2(p.si.*obj.cosdelta) + ...
                p.si.*Gamma2(p.si.*obj.cosdelta).*obj.cosdelta);
            return
        end

        function H2 = H2(obj)
            Gamma2 = @(x) (1 ./ sqrt(pi)) .* exp(-x.^2) + x .* (1 + erf(x));
            H2 = 1./(obj.parameters.si)...
                .*Gamma2(obj.parameters.si.*obj.cosdelta);
        end
    end
end