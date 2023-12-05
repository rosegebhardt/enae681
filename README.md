# Optimization of a vortex position near a hydrofoil to maximize resulting thrust

parameters.m : Defines physical parameters of the system. Includes hydrofoil size and shape parameters, uniform flow speed and angle-of-attack, water pressure, and vortex strength.

dragFunction.m : Computes the net force in the in-stream direction (force = drag - thrust) given a vortex position.

penaltyFunction.m : Computes a penalty term to ensure a minimum is not found inside the foil (prevents a non-physical result).

outputFunction.m : Defines a function used in fminunc to return the state and cost function evaluation at each iteration of the algorithm.

main.m : Runs the optimization algorithm.

visuals.m : Creates figures that show the flow field, pressure field, and cost functions over a domain for a given configuration as well as the resulting search path and learning curve of a specific run of the optimization.
