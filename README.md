# Symbolic Satellite Aerodynamics Toolbox
This toolbox aims to simplify the process quickly obtaining aerodynamic models for satellites with simple geometries using panel methods while allowing the user to define all possible variables as Matlab symbolic expressions.
The toolbox can calculate forces and torques of a user defined geometry using different aerodynamic models for Free Molecular Flow (e.g. Sentman).

This toolbox not only evaluates forces numerically for given inputs, but actually allows the user to obtain analytical aerodynamic models (e.g. calculation of forces) as Matlab Function handles.

### What you can do
Using symbolic variables allows for workflows such as:
- Defining a satellite bus with a symbolic length `l`
- Defining a satellite with rotatable panel geometries by defining rotation angles as symbolic variables
- Obtaining the torque of a satellite with symbolic variables as a matlab Function handle

### What this toolbox is CAN NOT DO
- Computing wether a panel is shadowed by other panels
- High fidelity computations

## Requirements
- Install Matlab with the [Matlab Symbolic Math Toolbox](https://mathworks.com/products/symbolic.html)
- Install [git](https://git-scm.com/downloads)

## Install
Open the terminal and run:
```bash
git submodule update --init --recursive
```

## Add to path
You can either make sure the `+saero` folder in available on the Matlab path or simply run
```matlab
openProject('.')
```
in the root project folder. 

Alternatively double click the `Symbolicvleoaerodynamics.prj` file to setup the project with correct path dependencies.

----------------------------

## Coordinate Systems
All coordinates are defined in the body reference frame. All normal vectors are assumed to be unit vectors. The `B`-frame is defined as a right hand system.

## Usage

### Satellite Geometries
Each satellite has a `SatelliteGeometry` which consists of (potentially several) `PanelGroup` objects.

As a user you need to build your satellites geometry by defining these panel groups.

### Panel group
You can simply use predefined shapes like

```matlab
l = 2;
w = 1;
h = 1;
center_of_mass = [0;0;0];
body = saero.geometry.shapes.Box(l, w, h, center_of_mass);
```

Alternatively you can simply manually define a satellite geometry by providing panel normals, positions of each panels center of pressure and panel areas:

```
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

body = saero.geometry.PanelGroup(cop_pos, normals, areas);
```

### Geometry
then you can simply define the geometry with
```matlab
satellite_geometry = saero.geometry.SatelliteGeometry([body])
```

## Calculation Methods
Each satellite has a calculation method to obtain aerodynamic forces (and torques).

Currently available:
- Sentman's method [1]
- Schaaf-Chambre [2]

The calculation method is initialized via
```matlab
sentman = saero.aerodynamics.Sentman();
```
where by default, a set of parameters is defined. Parameters can be optionally varied using the syntax
```matlab
sentman = saero.aerodynamics.Sentman( ...
    "alpha_E", 0.95, ...
    "rho", 1e-11);
```
where all other model parameters are kept at the default value.


## Defining the satellite
Each satellite is an instance of the class `saero.Satellite`. The following class diagram shows the hierarchy:
![class hierarchy](docs/images/Class_hierarchy.png)


## Full example
We will derive the aerodynamic model symbolically for the following satellite geometry with a variable size bus and two solar panels:
![satellite example](docs/images/satellite_example.png)

We start by defining all symbolic variables
```matlab
l = sym('l', 'real'); % symbolic bus length
w = sym('w', 'real'); % symbolic bus width
h = sym('h', 'real'); % symbolic bus height
d = sym('d', 'real'); % symbolic distance body cop -> wing cop
l_w = sym('l_w', 'real'); % symbolic wing length
w_w = sym('w_w', 'real'); % symbolic wing width
```

Define satellite bus using the predefined box shape

```matlab
bus = saero.geometry.shapes.Box(l, w, h, [0;0;0]); % satellite bus as Panel Group instance
```

Define satellite wings
```matlab
% Wing center of pressure in body coordinates
cop_wings = [0,  0, 0, 0;
            -d, -d, d, d;
             0,  0, 0, 0];

% Normal vectors of wings
normals_wings = [0,  0,  0,  0;
                 0,  0,  0,  0;
                 -1, 1, -1,  1];

% Symbolic wing area for each panel
wing_areas = l_w*w_w.*ones(1,4);

% Define wings as one panel group
wings = saero.geometry.PanelGroup(cop_wings, normals_wings, wing_areas);
```

and define the full geometry with
```matlab
% Full geometry (Bus + Wings)
sat_geometry = saero.geometry.SatelliteGeometry([bus; wings]);
```

we now have the object `sat_geometry` which clearly defines our satellites geometry. We can inspect the instance in the Matlab terminal

```matlab
>> sat_geometry

sat_geometry = 

  SatelliteGeometry with properties:

    panel_groups: [2×1 saero.geometry.PanelGroup]
```

Now we can define a calculation method to obtain aerodynamic forces and torques. Let's use Sentman's method

```matlab
aero = saero.aerodynamics.Sentman();
```

Every aerodynamic calculation method needs specific parameters. We can check them with
```matlab
>> aero.parameters

ans = 

  struct with fields:

    alpha_E: 0.9500
         si: 5
         Tw: 300
         Vi: 7800
        rho: 1.0000e-10
         kB: 1.3806e-23
         mT: 2.6569e-26
```

let's say we want to change the density to a symbolic value `rho`. We do this with
```matlab
aero.parameters.rho = sym('rho', 'real');
```

Finally, we use the geometry and the calculation model to define the `saero.Satellite` instance:
```matlab
% Full satellite
sat = saero.Satellite( ...
    "calculation_model",aero, "satellite_geometry", sat_geometry);
```

Now we can use the class methods to calculate forces and torques. For example we can obtain the symbolic expression of the total torque with
```matlab
% We have to define the incoming velocity vector of length 1:
% Here, we do it symbolically:
incoming_velocity = sym('vi', [3,1]);

% Obtain symbolic torque expression
torqueExpr = sat.get_total_aerodynamic_torque(incoming_velocity)

% Optional: turn into Matlab Function handle:
torqueFun = matlabFunction(torqueExpr)
```

This way of defining satellite models yields much freedom. For example we can use a transformation matrix to introduce Euler angles with the included `symbolic-space-math-utils` toolbox:

```matlab
% Define Aerodynamic angles
alpha = sym('a', 'real');
beta = sym('b', 'real');

% Transformation Matrix from aerodynamic frame to body frame
T_fa = ssmu.dcm.T2(alpha)*ssmu.dcm.T3(-beta);

% Transform wind vector based on current attitude w.r.t. wind
vi_B = T_fa*[1;0;0]

% Obtain torque as function of alpha/beta
torqueAB = sat.get_total_aerodynamic_torque(vi_B)
```

## References

1. H. Sentman, “FREE MOLECULE FLOW THEORY AND ITS APPLICATION TO THE DETERMINATION OF AERODYNAMIC FORCES”.
2. S. A. Schaaf and P. L. Chambre, Flow of Rarefied Gases. Princeton University Press, 1958.
