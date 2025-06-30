classdef Box < saero.geometry.PanelGroup
    %BOX Represents a 3D rectangular box composed of six flat panels

    properties (SetAccess = private)
        length (1,1)
        width  (1,1)
        height (1,1)
    end

    methods
        function obj = Box(length, width, height, center_of_mass)
            arguments
                length (1,1)
                width  (1,1)
                height (1,1)
                center_of_mass (3,1) 
            end

            % Store dimensions
            l = length;
            w = width;
            h = height;

            % Panel normals (unit vectors)
            normals = [  1,-1, 0, 0, 0, 0;
                         0, 0, 1, 0,-1, 0;
                         0, 0, 0,-1, 0, 1];

            % Panel areas
            areas = [h*w, h*w, h*l, w*l, h*l, w*l];

            % Center-of-pressure positions (midpoint of each panel)
            cop_pos = - center_of_mass + ...
                            [l/2, -l/2,  0,    0,   0,      0;
                            0,   0,     w/2,  0,   -w/2,   0;
                            0,   0,     0,    -h/2,   0,     h/2];

            % Call superclass constructor
            obj@saero.geometry.PanelGroup(cop_pos, normals, areas);

            % Store dimensions in the subclass
            obj.length = l;
            obj.width = w;
            obj.height = h;
        end
    end
end
