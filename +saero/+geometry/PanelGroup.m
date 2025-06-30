classdef PanelGroup < matlab.mixin.Heterogeneous
    %PANELGROUP group of one-sided panels defining a satellite
    
    properties (SetAccess = private)
        cop_positions (3,:)  % 3xN matrix of center-of-pressure positions
        panel_normals (3,:)  % 3xN matrix of unit normals
        panel_areas (1,:)    % 1xN vector of panel areas
    end
    
    methods
        function obj = PanelGroup(cop_positions, panel_normals, panel_areas)
            %PANELGROUP Construct an instance of this class
            %   A Panel group defines a group of panels that is a part of a
            %   satellite geometry. Currently, rotation and translation of
            %   panel groups is not implemented.
            obj.cop_positions = cop_positions;
            obj.panel_normals = panel_normals;
            obj.panel_areas = panel_areas;
        end
        
        %% Getters
        function normals = get_normals(obj)
            %get_normals returns array of normals
            normals = obj.panel_normals;
        end

        function cop_positions = get_cop_positions(obj)
            %get_cop_positions returns array of panel centers of pressure
            cop_positions = obj.cop_positions;
        end

        function areas = get_panel_areas(obj)
            %get_panel_areas returns array of panel areas
            areas = obj.panel_areas;
        end

        %% Setters (return updated object)
        function obj = set_cop_positions(obj, cop_positions)
            obj.cop_positions = cop_positions;
        end

        function obj = set_panel_normals(obj, panel_normals)
            obj.panel_normals = panel_normals;
        end

        function obj = set_panel_areas(obj, panel_areas)
            obj.panel_areas = panel_areas;
        end
    end
end
