% Creates a symbolic box

% box geometry
h = sym('h', 'real');
w = sym('w', 'real');
l = sym('l', 'real');

com = [0;0;0]; % from the middle of the box

body = saero.geometry.shapes.Box(l,w,h,com);

geometry = saero.geometry.SatelliteGeometry([body]);


calculation_method = saero.aerodynamics.Sentman();

satellite = saero.Satellite( ...
    "satellite_geometry", geometry, ...
    "calculation_model", calculation_method);


% Forces and torques as symbolic expressions
f = satellite.get_total_aerodynamic_force([-1;0;0]);
tau = satellite.get_total_aerodynamic_torque([-cos(pi/2);0;sin(pi/2)]);