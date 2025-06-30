classdef SatelliteGeometry < handle
    %SATELLITEGEOMETRY Abstract class defining a satellite geometry
    
    properties (SetAccess=private)
        panel_groups (:,1) saero.geometry.PanelGroup
    end
    
    methods
        function obj = SatelliteGeometry(panel_groups)
            %SATELLITEGEOMETRY Construct an instance of this class
            %   Holds information about the satellites geometry independent
            %   of the environment, thermals or any other satellite data
            if nargin > 0
                obj.panel_groups = panel_groups;
            end
        end

        % function disp(obj)
        %     fprintf('Object of class %s\n', class(obj));
        %     % Display selected private properties
        %     fprintf('Panels: %s\n', obj.panel_groups);
        % end

        function add_panel_group(obj, panel_groups)
            arguments
                obj
                panel_groups (:,1) saero.geometry.PanelGroup
            end

            obj.panel_groups = [obj.panel_groups; panel_groups];
        end
        
        function normals = get_normals(obj)
            normals = [];
            for i=1:numel(obj.panel_groups)
                normals = [normals, obj.panel_groups(i).get_normals()];
            end
        end

        function cop = get_cop_positions(obj)
            cop = [];
            for i=1:numel(obj.panel_groups)
                cop = [cop, obj.panel_groups(i).get_cop_positions()];
            end
        end
        
        function areas = get_panel_areas(obj)
            areas = [];
            for i=1:numel(obj.panel_groups)
                areas = [areas, obj.panel_groups(i).get_panel_areas()];
            end
        end
    end
end

