## Reflections

### The Model

The model state includes 6 parts describing the vehicle:

1. x position
2. y position
3. orientation (ψ)
4. velocity
5. cross track error (cte)
6. orientation error (eψ)

The actuator outputs of the model are:

1. steering angle
2. acceleration (converted to a throttle value)

The model uses the current state and actuator values to calculate the next state using the update equations:

1. x[t+1] = x[t] + v[t] * cos(ψ[t]) * dt
2. y[t+1] = y[t] + v[t] * sin(ψ[t]) * dt
3. ψ[t+1] = ψ[t] - v[t] / Lf * δ[t] * dt
4. v_[t+1] = v[t] + a[t] * dt
5. cte[t+1] = f(x[t]) - y[t] + v[t] * sin(eψ[t]) * dt
6. eψ[t+1] = ψ[t] - ψdes[t] - v[t] * δ[t] / Lf * dt

In these equations, **δ** is the steering angle.

### Timestep Length and Elapsed Duration (N & dt)

I chose **N** to be as large as possible without slowing down the solution calculation. I chose **dt** to be as small as possible while keeping the product `N * dt` as large as possible.

### Polynomial Fitting and MPC Preprocessing

I converted the waypoints to be in the car's coordinate system.

This greatly simplified the initial vehicle state, as the position and orientation could all be zero.

Additionally, the cross track error becomes just the constant term of the polynomial that describes the waypoints, as x is 0. Furthermore, the derivative of the polynomial, which is used in calculating the orientation, simplifies to just the coefficient of the 1st order x term.

### Model Predictive Control with Latency

To account for latency, I took the advice of my reviewer and used the kinematic equations to calculate the initial state of the vehicle 100ms into the future. Then the actuators calculated would be as if there was a 100ms latency between calculation and the actual actuation.