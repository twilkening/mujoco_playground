#include <stdio.h>
#include <math.h>
#include <string.h>

// Include the same parameter structure for consistency
typedef struct {
    // Simulation parameters
    double timestep;
    
    // Wheel parameters
    double wheel_radius;
    double wheel_thickness;
    double wheel_mass;
    double wheel_separation;
    
    // Pendulum parameters
    double pendulum_length;
    double pendulum_width;
    double pendulum_thickness;
    double pendulum_mass;
    
    // Cart parameters
    double cart_height;
    double cart_pitch_offset;
    
    // Motor parameters
    double motor_gear_ratio;
    double motor_force_limit;
    
    // Additional physical parameters for control
    double gravity;
    double friction_coefficient;
    double motor_inertia;
    
} CartPendulumParams;

// Function to load parameters (can be from file or hardcoded)
CartPendulumParams load_system_parameters() {
    CartPendulumParams params = {
        // Simulation
        .timestep = 0.001,
        
        // Wheels
        .wheel_radius = 0.05,           // 5cm radius
        .wheel_thickness = 0.005,       // 0.5cm thickness
        .wheel_mass = 0.01,             // 10g each
        .wheel_separation = 0.12,       // 12cm apart
        
        // Pendulum
        .pendulum_length = 0.25,        // 25cm length
        .pendulum_width = 0.02,         // 2cm width
        .pendulum_thickness = 0.1,      // 10cm thickness
        .pendulum_mass = 0.5,           // 500g
        
        // Cart
        .cart_height = 0.15,            // 15cm above ground
        .cart_pitch_offset = -0.1,      // pitch joint offset
        
        // Motors
        .motor_gear_ratio = 1.0,
        .motor_force_limit = 10.0,      // 10N max force
        
        // Physics
        .gravity = 9.81,                // m/s^2
        .friction_coefficient = 0.02,   // bearing friction
        .motor_inertia = 1e-4          // kg*m^2
    };
    return params;
}

// Function to calculate system matrices based on parameters
void calculate_system_matrices(CartPendulumParams params, 
                             double A[4][4], double B[4]) {
    
    double g = params.gravity;
    double m = params.pendulum_mass;
    double m_c = params.wheel_mass * 2; // total wheel mass
    double l = params.pendulum_length / 2; // half length to CoM
    double mu_c = params.friction_coefficient;
    double r = params.wheel_radius;
    
    // Simplified moment of inertia (rod about end)
    double J = (1.0/3.0) * m * l * l;
    
    // Linearized system matrices around vertical equilibrium
    // State: [x, x_dot, theta, theta_dot]
    // Input: u (motor force)
    
    double denominator = m_c * m * l * l + J * m + J * m_c;
    
    // A matrix (4x4)
    // dx/dt = x_dot
    A[0][0] = 0; A[0][1] = 1; A[0][2] = 0; A[0][3] = 0;
    
    // d(x_dot)/dt
    A[1][0] = 0;
    A[1][1] = -(g * mu_c * l * l * m * m + g * m_c * mu_c * l * l * m + 
                J * g * mu_c * m + J * g * m_c * mu_c) / (r * denominator);
    A[1][2] = -(g * l * l * m * m) / denominator;
    A[1][3] = 0;
    
    // dtheta/dt = theta_dot
    A[2][0] = 0; A[2][1] = 0; A[2][2] = 0; A[2][3] = 1;
    
    // d(theta_dot)/dt
    A[3][0] = 0;
    A[3][1] = (l * m * (g * m * mu_c + g * m_c * mu_c)) / (r * denominator);
    A[3][2] = (l * m * (g * m * r + g * m_c * r)) / (r * denominator);
    A[3][3] = 0;
    
    // B matrix (4x1)
    B[0] = 0;
    B[1] = (m * r * l * l + J * r) / (r * denominator);
    B[2] = 0;
    B[3] = -(l * m) / denominator;
}

// Function to calculate optimal PID gains based on system parameters
void calculate_pid_gains(CartPendulumParams params, 
                        double* kp, double* ki, double* kd) {
    
    // Simple pole placement approach
    double l = params.pendulum_length / 2;
    double m = params.pendulum_mass;
    double g = params.gravity;
    
    // Natural frequency of pendulum
    double omega_n = sqrt(g / l);
    
    // Desired closed-loop characteristics
    double zeta = 0.7;          // damping ratio
    double omega_d = omega_n * 2; // desired frequency (faster than natural)
    
    // Calculate gains based on desired pole placement
    *kp = omega_d * omega_d * m * l / params.motor_gear_ratio;
    *kd = 2 * zeta * omega_d * m * l / params.motor_gear_ratio;
    *ki = 0.1 * (*kp);  // Small integral gain to eliminate steady-state error
    
    printf("Calculated PID gains based on system parameters:\n");
    printf("  Kp = %.6f\n", *kp);
    printf("  Ki = %.6f\n", *ki);
    printf("  Kd = %.6f\n", *kd);
    printf("  System natural frequency: %.3f rad/s\n", omega_n);
    printf("  Desired frequency: %.3f rad/s\n", omega_d);
}

// Function to save parameters to file
void save_parameters_to_file(const char* filename, CartPendulumParams params) {
    FILE* file = fopen(filename, "w");
    if (!file) {
        printf("Error: Could not create parameter file %s\n", filename);
        return;
    }
    
    fprintf(file, "# Cart-Pendulum System Parameters\n");
    fprintf(file, "# Generated automatically\n\n");
    
    fprintf(file, "# Simulation\n");
    fprintf(file, "timestep = %.6f\n\n", params.timestep);
    
    fprintf(file, "# Wheel parameters\n");
    fprintf(file, "wheel_radius = %.6f\n", params.wheel_radius);
    fprintf(file, "wheel_thickness = %.6f\n", params.wheel_thickness);
    fprintf(file, "wheel_mass = %.6f\n", params.wheel_mass);
    fprintf(file, "wheel_separation = %.6f\n\n", params.wheel_separation);
    
    fprintf(file, "# Pendulum parameters\n");
    fprintf(file, "pendulum_length = %.6f\n", params.pendulum_length);
    fprintf(file, "pendulum_width = %.6f\n", params.pendulum_width);
    fprintf(file, "pendulum_thickness = %.6f\n", params.pendulum_thickness);
    fprintf(file, "pendulum_mass = %.6f\n\n", params.pendulum_mass);
    
    fprintf(file, "# Motor parameters\n");
    fprintf(file, "motor_gear_ratio = %.6f\n", params.motor_gear_ratio);
    fprintf(file, "motor_force_limit = %.6f\n\n", params.motor_force_limit);
    
    fprintf(file, "# Physics\n");
    fprintf(file, "gravity = %.6f\n", params.gravity);
    fprintf(file, "friction_coefficient = %.6f\n", params.friction_coefficient);
    
    fclose(file);
    printf("Saved parameters to: %s\n", filename);
}

int main() {
    printf("Cart-Pendulum Parameter Calculator\n");
    printf("==================================\n\n");
    
    // Load system parameters
    CartPendulumParams params = load_system_parameters();
    
    // Print current parameters
    printf("System Parameters:\n");
    printf("  Pendulum mass: %.3f kg\n", params.pendulum_mass);
    printf("  Pendulum length: %.3f m\n", params.pendulum_length);
    printf("  Wheel radius: %.3f m\n", params.wheel_radius);
    printf("  Wheel mass: %.3f kg (each)\n", params.wheel_mass);
    printf("  Motor force limit: %.1f N\n\n", params.motor_force_limit);
    
    // Calculate system matrices
    double A[4][4], B[4];
    calculate_system_matrices(params, A, B);
    
    printf("System A matrix (linearized):\n");
    for (int i = 0; i < 4; i++) {
        printf("  [");
        for (int j = 0; j < 4; j++) {
            printf("%8.4f ", A[i][j]);
        }
        printf("]\n");
    }
    
    printf("\nSystem B matrix:\n");
    printf("  [%8.4f %8.4f %8.4f %8.4f]^T\n\n", B[0], B[1], B[2], B[3]);
    
    // Calculate optimal PID gains
    double kp, ki, kd;
    calculate_pid_gains(params, &kp, &ki, &kd);
    
    // Save parameters to file
    save_parameters_to_file("system_parameters.txt", params);
    
    printf("\nYou can use these gains in your main.c file:\n");
    printf("  double kp = %.6f;\n", kp);
    printf("  double ki = %.6f;\n", ki);
    printf("  double kd = %.6f;\n", kd);
    
    return 0;
}
