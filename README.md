# Symbolic Satellite Aerodynamics Toolbox
This toolbox simplifies the process of quickly obtaining aerodynamic models for satellites with simple geometries using panel methods. It allows users to define all variables as Matlab symbolic expressions.
The toolbox can calculate forces and torques for user-defined geometries using different aerodynamic models for Free Molecular Flow (e.g., Sentman).

Unlike purely numerical tools, this toolbox enables users to obtain analytical aerodynamic models (e.g., force expressions) as Matlab function handles.

### What you can do
Using symbolic variables allows for workflows such as:
- Defining a satellite bus with a symbolic length `l`
- Defining a satellite with rotatable panel geometries by introducing symbolic rotation angles
- Obtaining the torque of a satellite with symbolic variables as a Matlab function handle

### What this toolbox is CANNOT do
- Computing whether a panel is shadowed by other panels
- Provide high-fidelity computations

## Requirements
- Matlab **with** the [Matlab Symbolic Math Toolbox](https://mathworks.com/products/symbolic.html)
- [git](https://git-scm.com/downloads)

## Installation
You can either:
- Download the `.zip` file of this repository to your preferred location **or**
- Clone the repository using git:  
  ```bash
  git clone https://github.com/ifrunistuttgart/symbolic-satellite-aerodynamics.git

(Optional if you want to use the `+ssmu` functionality):
Open the terminal in the root folder of this project and run
```bash
git submodule update --init --recursive
```

## Add to path
To ensure the toolbox is on the Matlab path, choose one of the following:
- Each time you open this repo in Matlab, double-click `Symbolicvleoaerodynamics.prj` (recommended) **or**
- make sure the `+saero` folder in available on the Matlab path **or** 
- run
```matlab
openProject('.')
```
in the root project folder. 

----------------------------

## Coordinate Systems
All coordinates are defined in the body reference frame (fixed with respect to the satellite geometry). All normal vectors are assumed to be unit vectors. The `B`-frame is defined as a right hand system.

## Usage

### Satellite Geometries

Each satellite has a SatelliteGeometry, which consists of one or more PanelGroup objects.
A PanelGroup is a collection of single-sided panels defined by:
- Center of pressure (3×1)
- Outward-pointing normal vector (3×1)
- Area (1×1)

These properties are combined into matrices (for centers/normals - 3xn) or vectors (for areas - 1xn), which may contain symbolic expressions.

Users construct their satellite geometries by defining such panel groups.

### Panel group
You can simply use predefined shapes like

```matlab
l = 2;
w = 1;
h = 1;
center_of_mass = [0;0;0];
body = saero.geometry.shapes.Box(l, w, h, center_of_mass);
```

Alternatively you can simply manually define a satellite panel group by providing panel normals, positions of each panels center of pressure and panel areas:

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
A geometry is basically just a list of panel groups.

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
where all other model parameters are kept at the default value. you can also see the current parameter settings with 

```matlab console
>> sentman.parameters
```
in the Matlab console.

Here is an overview of all parameters:
| Variable  | Description                                   | Method(s) Used          | SI Unit           |
| --------- | --------------------------------------------- | ----------------------- | ----------------- |
| `alpha_E` | Energy accommodation coefficient              | Sentman                 | – (dimensionless) |
| `si`      | Molecular speed ratio                         | Sentman; Schaaf–Chambre | – (dimensionless) |
| `Tw`      | Wall temperature                              | Sentman; Schaaf–Chambre | K                 |
| `Vi`      | Norm of incoming particle velocity            | Sentman; Schaaf–Chambre | m/s               |
| `rho`     | Density                                       | Sentman; Schaaf–Chambre | kg/m³             |
| `kB`      | Boltzmann constant                            | Sentman                 | J/K               |
| `mT`      | Particle mass                                 | Sentman                 | kg                |
| `sigma_n` | Normal momentum accommodation coefficient     | Schaaf–Chambre          | – (dimensionless) |
| `sigma_t` | Tangential momentum accommodation coefficient | Schaaf–Chambre          | – (dimensionless) |
| `T_inf`   | Free-stream temperature                       | Schaaf–Chambre          | K                 |



## Obtaining Torque/Force Expressions

Assume you defined a symbolic satellite `sat`. Define the incoming velocity vector, e.g., if facing the flow directly:

$$
v_i = [-1, 0, 0]^\top
$$

```matlab
vi = [-1;0;0];
forceExpr = sat.getTotalAerodynamicForce(vi);
torqueExpr = sat.getTotalAerodynamicTorque(vi);
```

> `vi` can also be symbolic (3×1). Ensure it is normalized when evaluating.

Turn expressions into Matlab function handles:

```matlab
forceFun = matlabFunction(forceExpr);
torqueFun = matlabFunction(torqueExpr);
```

See [Matlab’s `matlabFunction` docs](https://de.mathworks.com/help/symbolic/sym.matlabfunction.html).

**Key advantage:** Every satellite parameter can be symbolic, giving function handles that depend on those variables. This enables fast parametric studies of aerodynamic behavior.

---



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



# Known issues
- When using Git bash on Windows you might encounter **issues due to long file paths**. The solution is to setup git with `git config --system core.longpaths true` and to run `git submodule update --init --recursive` again