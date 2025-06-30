% Use Sentman model with standard parameters
aero = saero.aerodynamics.Sentman();

% Define satellite geometry
l = sym('l', 'real'); % symbolic bus length
w = sym('w', 'real'); % symbolic bus width
h = sym('h', 'real'); % symbolic bus height

% Define satellite bus using the predefined box shape
bus = saero.geometry.shapes.Box( ...
    l, w, h, [0;0;0]);


% Define center of pressure positions and normal vectors of the wings 
% (keep in mind: Wings have two sides -> 4 entries)
d = sym('d', 'real');
cop_wings = [0,  0, 0, 0;
            -d, -d, d, d;
             0,  0, 0, 0];
normals_wings = [0,  0,  0,  0;
                 0,  0,  0,  0;
                 -1, 1, -1,  1];
% Wing area
l_w = sym('l_w', 'real'); % symbolic wing length
w_w = sym('w_w', 'real'); % symbolic wing width
wing_areas = l_w*w_w.*ones(1,4);

wings = saero.geometry.PanelGroup(cop_wings, normals_wings, wing_areas);

% Full geometry (Bus + Wings)
sat_geometry = saero.geometry.SatelliteGeometry([bus; wings]);

%% Define full satellite model
sat = saero.Satellite( ...
    "calculation_model",aero, "satellite_geometry", sat_geometry);

% Turn density into symbolic variable aswell
sat.calculation_model.parameters.rho = sym('rho', 'real');

% symbolic incoming wind
vi = sym('vi', [3,1]);

% get forces in body frame
forceExpr = sat.get_total_aerodynamic_force(vi);

% Turn symexpr into matlab function handle
forceFun = matlabFunction(forceExpr);
