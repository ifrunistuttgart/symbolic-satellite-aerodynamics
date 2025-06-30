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
    end
end

