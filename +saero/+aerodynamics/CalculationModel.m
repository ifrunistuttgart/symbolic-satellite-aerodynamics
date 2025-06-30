classdef (Abstract) CalculationModel < handle
    %CALCULATIONMODEL Abstract class of an aerodynamic calculation model
    
    properties (Abstract)
        parameters
    end

    methods (Abstract)
        % Calculate the force vector for each panel:
        force_vectors = calculate_force(obj, wind_direction__B, ...
            satellite_geometry)

    end
end

