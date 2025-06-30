classdef SchaafChambre < saero.aerodynamics.CalculationModel
    %SENTMAN implements Calculation model with Sentman calculation method
    %   Detailed explanation goes here
    
    properties
        parameters
    end
    
    properties (Access=private)
        % cosdelta = sym('cosd', 'real');
        significant_digits = 5
    end
    
    methods
        function obj = SchaafChambre(parameters)
            %SENTMAN Construct an instance of this class
            %   Detailed explanation goes here
            arguments
                parameters.sigma_n = 0.95; % normal energy accommodation
                parameters.sigma_t = 0.95; % tangential energy accommodation
                parameters.si = 5; % Molecular speed ratio
                parameters.T_inf = 700; % Ambient Temperature
                parameters.Tw = 300; % Wall temperature
                parameters.Vi = 7800; % scalar velocity
                parameters.rho = 1e-10; % density
            end

            obj.parameters = parameters;
            
            % % Add constant values
            % obj.parameters.kB = 1.380649e-23;  % Boltzmann Constant
            % obj.parameters.mT = 16 * 1.6605390689252e-27;

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

            sigmaN = obj.parameters.sigma_n;
            sigmaT = obj.parameters.sigma_t;
            s = obj.parameters.si;
            Tw = obj.parameters.Tw;
            Tinf = obj.parameters.T_inf;

            % Assuming length one in normals and vi_B
            cosdelta = dot(-vi_B, N);
            sindelta = norm(cross(vi_B, N));

            
            cp   = 1/s^2.*(((2-sigmaN).*s/sqrt(pi).*cosdelta+sigmaN/2*...
                (Tw/Tinf)^0.5).*exp(-s^2.*(cosdelta).^2)+...
            ((2-sigmaN).*(0.5+s^2.*(cosdelta).^2)+sigmaN/2.*...
            (Tw/Tinf)^0.5.*sqrt(pi).*s.*cosdelta).*(1+erf(s.*cosdelta)));
            
            ctau = sigmaT.*sindelta/(s*sqrt(pi)).*(exp(-s^2.*...
                (cosdelta).^2)+ s.*sqrt(pi).*cosdelta.*...
                (1+erf(s.*cosdelta)));

            % ADBSat n_v = -vi_B
            tau = cross(N, cross(vi_B, N)); % shear direction

            % Force
            p = obj.parameters;
            force_vectors = A.*(ctau.*tau - cp.*N).*p.rho/2.*p.Vi^2;
        end
    end
end