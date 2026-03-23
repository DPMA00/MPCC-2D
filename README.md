# MPCC-2D (IN DEVELOPMENT)

A C++ library in development for **Model Predictive Contouring Control (MPCC)** of ground-based wheeled vehicles (2D).

The library will support:

* Differential drive vehicles
* Kinematic bicycle vehicles for Ackermann-style steering

Given a path or track description, the library allows construction of a **parametric path representation** that is passed to the MPCC. The controller then solves the resulting nonlinear optimal control problem using **Sequential Quadratic Programming (SQP)**.

## Example Workflow (future)

``test/`` will contain examples.

### Fixed track

1. Define a **ParametricSpline** object with desired spline type.
2. Load your path/track as a **waypoints** struct.
3. Construct an MPCC object with prediction horizon length **N**, sampling time **$T_s$** and the spline created above.
4. Select a vehicle model and integration method using **configure_dynamics()**.
5. Select a track projection method and configure its settings using the **configure_projection()**.
6. Pass the path into **update_path()** to construct the parametric path representation.
7. Start the MPC loop using your state feedback or state estimation pipeline.
8. At each control step, call **solve()** with the current state estimate.

### Dynamic Track (with a global planner)

1. Define a **ParametricSpline** object with desired spline type.
2. Construct an MPCC object with prediction horizon length **N**, sampling time **$T_s$** and the spline created above.
3. Select a vehicle model and integration method using **configure_dynamics()**.
4. Select a track projection method and configure its settings using the **configure_projection()**.
5. Start the global planner outer loop.
6. At each planning update, load the planned path as a **waypoints** struct and pass the path into **update_path()** to construct the parametric path representation.
1. Start the MPC loop using your state feedback or state estimation pipeline.
2. At each control step, call **solve()** with the current state estimate.



