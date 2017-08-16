#include "MPC.h"
#include <cppad/cppad.hpp>
#include <cppad/ipopt/solve.hpp>
#include "Eigen-3.3/Eigen/Core"
#include <float.h>

using CppAD::AD;

// TODO: Set the timestep length and duration
size_t N = 25;
double dt = 0.05;

// This value assumes the model presented in the classroom is used.
//
// It was obtained by measuring the radius formed by running the vehicle in the
// simulator around in a circle with a constant steering angle and velocity on a
// flat terrain.
//
// Lf was tuned until the the radius formed by the simulating the model
// presented in the classroom matched the previous radius.
//
// This is the length from front to CoG that has a similar radius.
const double Lf = 2.67;

const double target_v = 45.0;

size_t x_begin = 0;
size_t y_begin = x_begin + N;
size_t psi_begin = y_begin + N;
size_t v_begin = psi_begin + N;
size_t cte_begin = v_begin + N;
size_t epsi_begin = cte_begin + N;
size_t steer_begin = epsi_begin + N;
size_t throttle_begin = steer_begin + N - 1;

class FG_eval {
 public:
  // Fitted polynomial coefficients
  Eigen::VectorXd coeffs;
  FG_eval(Eigen::VectorXd coeffs) { this->coeffs = coeffs; }

  typedef CPPAD_TESTVECTOR(AD<double>) ADvector;
  void operator()(ADvector& fg, const ADvector& vars) {
    // TODO: implement MPC
    // `fg` a vector of the cost constraints, `vars` is a vector of variable values (state & actuators)
    // NOTE: You'll probably go back and forth between this function and
    // the Solver function below.

    fg[0] = 0.0;

    // Cost based on reference state
    for (int i = 0; i < N; i++) {
      // Minimize cross track error
      fg[0] += CppAD::pow(vars[cte_begin + i], 2);

      // Minimize orientation error
      fg[0] += CppAD::pow(vars[epsi_begin + i], 2);

      // Keep velocity as close to the target velocity as possible
      fg[0] += CppAD::pow(vars[v_begin + i] - target_v, 2);
    } 

    // Cost based on change in actuators to make a smoother ride
    for (int i = 1; i < N - 1; i++) {
      // Minimize changes in steering
      fg[0] += CppAD::pow(vars[steer_begin + i] - vars[steer_begin + i - 1], 2);

      // Minimize changes in throttle
      fg[0] += CppAD::pow(vars[throttle_begin + i] - vars[throttle_begin + i - 1], 2);
    }

    // Initial state
    fg[1 + x_begin] = vars[x_begin];   
    fg[1 + y_begin] = vars[y_begin];   
    fg[1 + psi_begin] = vars[psi_begin];   
    fg[1 + v_begin] = vars[v_begin];   
    fg[1 + cte_begin] = vars[cte_begin];   
    fg[1 + epsi_begin] = vars[epsi_begin];   

    int latency = round(0.1 / dt);
    // Setup constraints
    for (int i = 0; i < N - 1; i++) {
      // Current state
      AD<double> x0 = vars[x_begin + i];
      AD<double> y0 = vars[y_begin + i];
      AD<double> psi0 = vars[psi_begin + i];
      AD<double> v0 = vars[v_begin + i];
      AD<double> cte0 = vars[cte_begin + i];
      AD<double> epsi0 = vars[epsi_begin + i];

      AD<double> steer0 = vars[steer_begin + i];
      AD<double> throttle0 = vars[throttle_begin + i];

      if (i > (latency - 1)) {
        throttle0 = vars[throttle_begin + i - latency];
        steer0 = vars[steer_begin + i - latency];
      }

      // Next state
      AD<double> x1 = vars[x_begin + i + 1];
      AD<double> y1 = vars[y_begin + i + 1];
      AD<double> psi1 = vars[psi_begin + i + 1];
      AD<double> v1 = vars[v_begin + i + 1];
      AD<double> cte1 = vars[cte_begin + i + 1];
      AD<double> epsi1 = vars[epsi_begin + i + 1];

      // Calculate the desired position 
      AD<double> f_x = coeffs[0];

      // Calculate the derivative at the desired position
      AD<double> f_dx = 0.0;
      AD<double> x = 1.0;
      for (int j = 1; j < coeffs.size(); j++) {
        f_dx += j * coeffs[j] * x;
        x = x * x0;
        f_x += coeffs[j] * x;
      }

      // Calculate the desired orientation
      AD<double> psi_des = CppAD::atan(f_dx);

      AD<double> v0_dt = v0 * dt;
      AD<double> next_psi = psi0 - v0_dt * steer0 / Lf;

      fg[1 + x_begin + i + 1] = x1 - (x0 + CppAD::cos(psi0) * v0_dt);
      fg[1 + y_begin + i + 1] = y1 - (y0 + CppAD::sin(psi0) * v0_dt);
      fg[1 + psi_begin + i + 1] = psi1 - next_psi;
      fg[1 + v_begin + i + 1] = v1 - (v0 + throttle0 * dt);
      fg[1 + cte_begin + i + 1] = cte1 - (f_x - y0 + CppAD::sin(epsi0) * v0_dt);
      fg[1 + epsi_begin + i + 1] = epsi1 - (next_psi - psi_des);
    }
  }
};

//
// MPC class definition implementation.
//
MPC::MPC() {}
MPC::~MPC() {}

vector<double> MPC::Solve(Eigen::VectorXd state, Eigen::VectorXd coeffs) {
  bool ok = true;
  size_t i;
  typedef CPPAD_TESTVECTOR(double) Dvector;

  // TODO: Set the number of model variables (includes both states and inputs).
  // For example: If the state is a 4 element vector, the actuators is a 2
  // element vector and there are 10 timesteps. The number of variables is:
  //
  // 4 * 10 + 2 * 9
  size_t n_vars = 6 * N + 2 * (N - 1);
  // TODO: Set the number of constraints
  size_t n_constraints = 6 * N;

  double x = state[0];
  double y = state[1];
  double psi = state[2];
  double v = state[3];
  double cte = state[4];
  double epsi = state[5];

  // Initial value of the independent variables.
  // SHOULD BE 0 besides initial state.
  Dvector vars(n_vars);
  for (int i = 0; i < n_vars; i++) {
    vars[i] = 0;
  }

  vars[x_begin] = x;
  vars[y_begin] = y;
  vars[psi_begin] = psi;
  vars[v_begin] = v;
  vars[cte_begin] = cte;
  vars[epsi_begin] = epsi;

  Dvector vars_lowerbound(n_vars);
  Dvector vars_upperbound(n_vars);
  // TODO: Set lower and upper limits for variables.

  for (int i = 0; i < steer_begin; i++) {
    vars_lowerbound[i] = -1.0e19;
    vars_upperbound[i] = 1.0e19;
  }

  // Steering angle needs to in range [-25, 25] degrees, but should be in radians
  double radians25 = M_PI * 25.0 / 180.0;
  for (int i = steer_begin; i < throttle_begin; i++) {
    vars_lowerbound[i] = -radians25;
    vars_upperbound[i] = radians25;
  }

  // Acceleration should be between [-1, 1]
  for (int i = throttle_begin; i < n_vars; i++) {
    vars_lowerbound[i] = -1.0;
    vars_upperbound[i] = 1.0;
  }

  // Lower and upper limits for the constraints
  // Should be 0 besides initial state.
  Dvector constraints_lowerbound(n_constraints);
  Dvector constraints_upperbound(n_constraints);
  for (int i = 0; i < n_constraints; i++) {
    constraints_lowerbound[i] = 0;
    constraints_upperbound[i] = 0;
  }
;
  constraints_lowerbound[x_begin] = x;
  constraints_lowerbound[y_begin] = y;
  constraints_lowerbound[psi_begin] = psi;
  constraints_lowerbound[v_begin] = v;
  constraints_lowerbound[cte_begin] = cte;
  constraints_lowerbound[epsi_begin] = epsi;

  constraints_upperbound[x_begin] = x;
  constraints_upperbound[y_begin] = y;
  constraints_upperbound[psi_begin] = psi;
  constraints_upperbound[v_begin] = v;
  constraints_upperbound[cte_begin] = cte;
  constraints_upperbound[epsi_begin] = epsi;

  // object that computes objective and constraints
  FG_eval fg_eval(coeffs);

  //
  // NOTE: You don't have to worry about these options
  //
  // options for IPOPT solver
  std::string options;
  // Uncomment this if you'd like more print information
  options += "Integer print_level  0\n";
  // NOTE: Setting sparse to true allows the solver to take advantage
  // of sparse routines, this makes the computation MUCH FASTER. If you
  // can uncomment 1 of these and see if it makes a difference or not but
  // if you uncomment both the computation time should go up in orders of
  // magnitude.
  options += "Sparse  true        forward\n";
  options += "Sparse  true        reverse\n";
  // NOTE: Currently the solver has a maximum time limit of 0.5 seconds.
  // Change this as you see fit.
  options += "Numeric max_cpu_time          0.5\n";

  // place to return solution
  CppAD::ipopt::solve_result<Dvector> solution;

  // solve the problem
  CppAD::ipopt::solve<Dvector, FG_eval>(
      options, vars, vars_lowerbound, vars_upperbound, constraints_lowerbound,
      constraints_upperbound, fg_eval, solution);

  // Check some of the solution values
  ok &= solution.status == CppAD::ipopt::solve_result<Dvector>::success;

  // Cost
  auto cost = solution.obj_value;
  std::cout << "Cost " << cost << std::endl;

  for (int i = 0; i < N; i++) {
    std::cout << "x: " << solution.x[x_begin + i] << ", y: " << solution.x[y_begin + i] << ", psi: " << solution.x[psi_begin + i] << ", v: " << solution.x[v_begin] << ", steer: " << solution.x[steer_begin] << ", throttle: " << solution.x[throttle_begin] << std::endl;
  }

  // TODO: Return the first actuator values. The variables can be accessed with
  // `solution.x[i]`.
  //
  // {...} is shorthand for creating a vector, so auto x1 = {1.0,2.0}
  // creates a 2 element double vector.

  vector<double> ret;
  ret.push_back(solution.x[steer_begin]);
  ret.push_back(solution.x[throttle_begin]);

  for (int i = 1; i < N; i++) {
    ret.push_back(solution.x[x_begin + i]);
    ret.push_back(solution.x[y_begin + i]);
  }

  return ret;
}
