classdef Satellite < handle
    %SATELLITE Summary of this class goes here
    %   Detailed explanation goes here
    
    properties (SetAccess=private)
        calculation_model (1,1)
        geometry (1,1) saero.geometry.SatelliteGeometry
    end
    
    methods
        function obj = Satellite(options)
            %SATELLITE Construct a satellite instance
            arguments
                options.calculation_model (1,1) ...
                    {mustBeA(options.calculation_model, ...
                    'saero.aerodynamics.CalculationModel')} ...
                    = saero.aerodynamics.Sentman();
                options.satellite_geometry (1,1) ...
                    {mustBeA(options.satellite_geometry, ...
                    'saero.geometry.SatelliteGeometry')} ...
                    = saero.geometry.SatelliteGeometry();
            end
            obj.calculation_model = options.calculation_model;
            obj.geometry = options.satellite_geometry;
        end
        
        

        % Force of each panel
        function force_B = get_aerodynamic_force(obj, wind_direction__B)
            force_B = obj.calculation_model.calculate_force( ...
                wind_direction__B, obj.geometry);
        end
        
        % Force of all panels
        function total_force_B = get_total_aerodynamic_force(obj, ...
                wind_direction__B)
            force_B_vectors = obj.calculation_model.calculate_force( ...
                wind_direction__B, obj.geometry);
            total_force_B = sum(force_B_vectors, 2);
        end

        % Torque of each panel
        function torque_B = get_aeroynamic_torque(obj, ...
                wind_direction__B)
            force_vectors = obj.get_aerodynamic_force(wind_direction__B);
            torque_B = cross(obj.geometry.get_cop_positions, force_vectors);
        end

        % Torque of all panels
        function torque_B = get_total_aerodynamic_torque(obj, ...
                wind_direction__B)
            torque_B = sum(obj.get_aeroynamic_torque(wind_direction__B), 2);
        end

        % Compute aerodynamic equilibria
        function [equilibria_table] = ...
            get_aerodynamic_equilibria(obj, input_ranges)
            equilibria_table = table(fieldnames(input_ranges)); % nx2 matrix
            
            % Parametrize using alpha and beta
            alpha = sym("alpha");
            beta = sym("beta");
            
            % aerodynamic direction in body coordinates using alpha/beta
            vi_B = ssmu.dcm.aeroToBody(alpha, beta)*[-1;0;0];
            
            % get torque expression with alpha/beta parametrization
            torque_expr = obj.get_total_aerodynamic_torque(vi_B);

            % extract additional variables (without alpha and beta) in
            % alphabetical order
            extra_vars = reshape(symvar(torque_expr), 1, []);
            extra_vars(ismember(extra_vars, [alpha, beta])) = [];

            % make sure the input ranges struct is ordered alphabetically
            input_ranges = orderfields(input_ranges, ...
                sort(fieldnames(input_ranges)));

            % make sure user has provided input ranges for all symvars
            if ~isequal(...
                    reshape(string(fieldnames(input_ranges)), [], 1), ...
                    reshape(string(extra_vars),[],1))
            
                expected = strjoin(string(extra_vars), ", ");
                defined  = strjoin(string(fieldnames(input_ranges)), ", ");
            
                error("Mismatch in variable definitions.\n" + ...
                    "Expected (extra_vars): %s\n" + ...
                    "Defined (struct fields): %s", ...
                    expected, defined);
            end

            
        end
    end
end

