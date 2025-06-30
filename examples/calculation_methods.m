
% Two-sided panel
delta = sym('delta', 'real');
normals = [ sin(delta),  -sin(delta);
            0,           0;
            cos(delta),  -cos(delta)];
panel = saero.geometry.PanelGroup(zeros(3,2), normals, [1,1]);
panel = saero.geometry.SatelliteGeometry([panel]);

% Two different models
sentman = saero.aerodynamics.Sentman( ...
    "alpha_E", 0.95, ...
    "rho", 1e-11);
schaafChambre = saero.aerodynamics.SchaafChambre(...
    "sigma_n",1,...
    "rho", 1e-11);

% Two different satellites
sentmanSat = saero.Satellite("calculation_model", sentman, ...
    "satellite_geometry", panel);
schaafSat = saero.Satellite("calculation_model", schaafChambre, ...
    "satellite_geometry", panel);
vi_B = [-1;0;0];
sentmanForce = sentmanSat.get_total_aerodynamic_force(vi_B);
schaafForce = schaafSat.get_total_aerodynamic_force(vi_B);

%% Plot forces
delta = linspace(-1,1,100);
figure(1);
clf;
hold on;
sFX = matlabFunction(sentmanForce(1));
sFZ = matlabFunction(sentmanForce(3));
schFX = matlabFunction(schaafForce(1));
schFZ = matlabFunction(schaafForce(3));

plot(delta, sFX(delta), 'DisplayName', 'Sentman f_x');
plot(delta, sFZ(delta), 'DisplayName', 'Sentman f_z');
plot(delta, schFX(delta), 'DisplayName', 'Schaaf f_x');
plot(delta, schFZ(delta), 'DisplayName', 'Schaaf f_z');

legend

hold off