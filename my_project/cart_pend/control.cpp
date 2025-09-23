#include "main.h"

/* define globals once here */
double simend = 20; 
    //related to writing data to a file
FILE* fid = nullptr;
int loop_index = 0;
const int data_frequency = 50;

char path[] = ""; // "../myproject/dbpendulum/";
char xmlfile[] = "cart_pend_meshes.xml";
char datafile[] = "data.csv";

// holders of one step history of time and position to calculate dertivatives
mjtNum position_history = 0;
mjtNum previous_time = 0;

// controller related variables
double ctrl_update_freq = 100.0;

// noise injection parameters
double ctrl_noise_std = 0.01;
double ctrl_noise_rate = 50.0;

// scale by r to convert from force input to torque
const double r = 0.04;

// gains for state-space controller
// mjtNum K_ssp[4] = {-1.2957*r, -2.7642*r, -6.7969*r, -0.5284*r}; // scale by r to convert from force input to torque
mjtNum K_ssp[4] = {-0.0079*r, -1.6837*r, -3.6745*r, -0.1261*r};

// gains for state-space LQR controller with desired state
mjtNum K_ssp_desired[4] = {-0.0088*r, -3.1188*r, -8.7115*r, -0.7069*r}; // Q = diag([1 1 10 100]); R = 1e4;
// mjtNum K_ssp_desired[4] = {-0.0760*r, -2.9516*r, -9.1199*r, -1.0618*r}; // Q = diag([1 1 10 100]); R = 1e2;
// mjtNum K_ssp_desired[4] = {-0.0881*r, -3.1477*r, -8.7474*r, -0.7102*r}; // Q = diag([100 1 10 100]); R = 1e4;
// mjtNum K_ssp_desired[4] = {-1.8645*r, -3.4969*r, -11.9430*r, -2.1176*r}; // Q = diag([100 1 10 100]); R = 1; // unstable
// mjtNum K_ssp_desired[4] = {-0.0088*r, -3.1188*r, -8.7122*r, -0.7069*r}; // Q = diag([1 1 100 100]); R = 1e4;
// mjtNum K_ssp_desired[4] = {-0.0088*r, -3.1217*r, -8.7108*r, -0.7014*r}; // Q = diag([1 1 1000 10]); R = 1e4;
// mjtNum K_ssp_desired[4] = {-0.0859*r, -3.1184*r, -9.4975*r, -0.7817*r}; // Q = diag([1 1 1000 10]); R = 1e2;
// mjtNum K_ssp_desired[4] = {-0.1883*r, -2.2166*r, -8.1008*r, -1.9945*r}; // Q = diag([1 1 10 100]); R = 1; // unstable
// mjtNum K_ssp_desired[4] = {-0.0760*r, -2.9516*r, -9.1199*r, -1.0618*r}; // Q = diag([1 1 10 100]); R = 1e2; // not as good
// mjtNum K_ssp_desired[4] = {-0.0124*r, -3.1168*r, -8.7223*r, -0.7134*r}; // Q = diag([1 1 10 100]); R = 5e3;
// mjtNum K_ssp_desired[4] = {-0.0195*r, -3.1099*r, -8.7514*r, -0.7318*r}; // Q = diag([1 1 10 100]); R = 2e3;
// mjtNum K_ssp_desired[4] = {-0.0273*r, -3.0981*r, -8.7937*r, -0.7602*r}; // Q = diag([1 1 10 100]); R = 1e3;

// gains for state-space controller with integral action
mjtNum K_sspi[5] = {-1.9153*r, -3.7396*r, -9.4674*r, -0.7756*r, 0.0088*r}; // Q = diag([1 1 10 100 1]); R = 1e4;


//**************************
// Helper function for Quaternion to Euler angles
// Source: https://en.wikipedia.org/wiki/Conversion_between_quaternions_and_Euler_angles
void QuaternionToEuler(const mjtNum* quat, mjtNum* euler) {
    // roll (x-axis rotation)
    double sinr_cosp = 2 * (quat[0] * quat[1] + quat[2] * quat[3]);
    double cosr_cosp = 1 - 2 * (quat[1] * quat[1] + quat[2] * quat[2]);
    euler[0] = atan2(sinr_cosp, cosr_cosp); // roll angle   

    // pitch (y-axis rotation)
    double sinp = sqrt(1 + 2 * (quat[0] * quat[2] - quat[1] * quat[3]));
    double cosp = sqrt(1 - 2 * (quat[0] * quat[2] - quat[1] * quat[3]));
    euler[1] = 2 * atan2(sinp, cosp) - M_PI / 2; // pitch angle

    // yaw (z-axis rotation)
    double siny_cosp = 2 * (quat[0] * quat[3] + quat[1] * quat[2]);
    double cosy_cosp = 1 - 2 * (quat[2] * quat[2] + quat[3] * quat[3]);
    euler[2] = atan2(siny_cosp, cosy_cosp); // yaw angle
}

//****************************
//This function is called once and is used to get the headers
void init_save_data()
{
    //write name of the variable here (header)
    fprintf(fid,"t, ");
    // fprintf(fid,"KineticEnergy, PotentialEnergy, TotalEnergy");
    fprintf(fid,"quat_w, quat_x, quat_y, quat_z, ");
    fprintf(fid,"pos_x, pos_y, pos_z, ");
    fprintf(fid,"gyro_x, gyro_y, gyro_z, ");
    fprintf(fid,"accel_x, accel_y, accel_z, ");
    fprintf(fid,"LeftControl, RightControl,");
    fprintf(fid,"Pitch, PitchVel, ");
    //Don't remove the newline
    fprintf(fid,"\n");
}

//***************************
//This function is called at a set frequency, put data here
void save_data(const mjModel* m, mjData* d)
{
    //data here should correspond to headers in init_save_data()
    //seperate data by a space %f followed by space
    fprintf(fid,"%f, ",d->time);
    // fprintf(fid,"%f, %f, %f, ",
    //     d->energy[0], d->energy[1], d->energy[0] + d->energy[1]);
    // sensor data from cart_quaternion (w,x,y,z)
    fprintf(fid,"%f, %f, %f, %f, ", d->sensordata[0], d->sensordata[1], d->sensordata[2], d->sensordata[3]);
    // sensor data from cart_position (x,y,z)
    fprintf(fid,"%f, %f, %f, ", d->sensordata[4], d->sensordata[5], d->sensordata[6]);
    // sensor data from gyro (x,y,z)
    fprintf(fid,"%f, %f, %f, ", d->sensordata[7], d->sensordata[8], d->sensordata[9]);
    // sensor data from accelerometer (x,y,z)
    fprintf(fid,"%f, %f, %f, ", d->sensordata[10], d->sensordata[11], d->sensordata[12]);
    // control signals for wheels
    fprintf(fid,"%f, %f, ", d->ctrl[0], d->ctrl[1]);
    // Pitch and PitchVel (x-axis)
    mjtNum quat[4] = {d->qpos[3], d->qpos[4], d->qpos[5], d->qpos[6]}; // using actual angle of pendulum
    mjtNum euler[3];
    QuaternionToEuler(quat, euler);
    fprintf(fid,"%f, %f, ", euler[0], d->qvel[3]);

    //Don't remove the newline
    fprintf(fid,"\n");
}

//**************************
// inject Brownian noise on the control signals
// really I should be injecting noise on the sensor signals, but this is a simple example
// and we want to see the effect of noise on the control signals
void InjectControlNoise(const mjModel* m, mjData* d) {
    // no noise, return
    if (ctrl_noise_std <= 0) {
        return;
    }

    // convert rate and scale to discrete time (Ornstein–Uhlenbeck)
    mjtNum rate = mju_exp(-m->opt.timestep / ctrl_noise_rate);
    mjtNum scale = ctrl_noise_std * mju_sqrt(1-rate*rate);

    for (int i=0; i<m->nu; i++) {
        mjtNum bottom = 0, top = 0, midpoint = 0, halfrange = 1;
        if (m->actuator_ctrllimited[i]) {
            bottom = m->actuator_ctrlrange[2*i];
            top = m->actuator_ctrlrange[2*i+1];
            midpoint =  0.5 * (top + bottom);  // target of exponential decay
            halfrange = 0.5 * (top - bottom);  // scales noise
        }

        // exponential convergence to midpoint at ctrl_noise_rate
        d->ctrl[i] = rate * d->ctrl[i] + (1-rate) * midpoint;

        // add noise
        d->ctrl[i] += scale * halfrange * mju_standardNormal(NULL);

        // clip to range if limited
        if (m->actuator_ctrllimited[i]) {
            d->ctrl[i] = mju_clip(d->ctrl[i], bottom, top);
        }
    }
}

/*
OK. so I'd like to implement a PID controller for the cart-pendulum system.
basically, we are linearizing the system around the equilibrium at the veritical position
then we convert to discrete time.
then we measure the angle of the pendulum from vertical
based on the error of the angle w.r.t. desired angle (need to measure what the desired angle is)
we apply the control of the system...
I don't want the controller to actually execute at the full simulation speed...
so we need to rate-limit to a desired control speed...
in the main function call, we set up a guard check for 100Hz.
This way I know roughly how often the controller will update, and can design my gains using that information
before I actually tune the system

for now, we assume that the wheels are commanded in sync, with the same force applied.

SISO system: angle of the pendulum 
(later on we can perform some smoothing to combine a Gyro and an accelerometer)
(for now we are effectively using a perfect gyro, by getting angle of pendulum)

inputs:
r - reference point (desired output)
y - actual measured angle
kp - proportional gain
kd - derivative gain
ki - integral gain
*/

double myPID(const mjModel* m, mjData* d, double r, double y, double kp, double ki, double kd)
{
    static double e_prior = 0; // history of position error
    static double t_prior = 0; // history of time
    static double integral = 0; // integral of error

    // calculate error
    double e = r - y;

    // calculate derivative of error
    // coarse one-step Euler method for derivative, may need to smooth out 
    // for better performance in the future
    double de = (e - e_prior) / (d->time - t_prior);

    // update history
    e_prior = e;
    t_prior = d->time;

    integral += e * (d->time - t_prior); // approximate integral using rectangle method

    // calculate control signal
    double ctrl = kp * e + kd * de + ki * integral;

    return ctrl;
}

//**************************
void myPIDcontroller(const mjModel* m, mjData* d)
{
    double r = 0.0; // desired angle (reference point)
    // double y = d->sensordata[0]; // sensed angle of the pendulum (e.g., from a gyro or accelerometer)
    // since qpos is now quaternion, we need to convert to angle
    mjtNum euler[3];
    // mjtNum quat[4] = {d->qpos[3], d->qpos[4], d->qpos[5], d->qpos[6]}; // using actual angle of pendulum
    mjtNum quat[4] = {d->sensordata[0], d->sensordata[1], d->sensordata[2], d->sensordata[3]}; // using sensed angle of pendulum
    QuaternionToEuler(quat, euler);
    double xtheta = euler[0]; // use the roll (x-axis) angle as the pendulum angle
    double kp = 600e-3; // proportional gain
    double kd = 20e-3; // derivative gain
    double ki = 60e-3; // integral gain

    static double last_update = 0.0; // last time the control was updated
    
    // guard check to update control at a fixed frequency
    if (d->time - last_update >= 1.0 / ctrl_update_freq)
    {
        // hard-coded PD controller example
        // d->ctrl[0] = -kp*(d->sensordata[0]-0)-kd*d->sensordata[1];
        // d->ctrl[1] = -d->ctrl[0];   // apply the inverse control to both wheels 

        // method-based PID controller
        double ctrl = myPID(m, d, r, xtheta, kp, ki, kd);
        d->ctrl[0] = -ctrl; // apply control to the first actuator (left wheel)
        d->ctrl[1] = ctrl; // apply control to the second actuator (right wheel)

        last_update = d->time;
    }


    // InjectControlNoise(); // inject noise into the control signal

    //write data here (dont change/dete this function call; instead write what you need to save in save_data)
    if ( loop_index%data_frequency==0){
        save_data(m,d);
    }
    loop_index = loop_index + 1;
}

//**************************
// State-space controller example (discrete time, 100Hz):
    // pseudo-code...
    // define A, B, C, D
    // c2d, get G and H matrices
    // define desired pole locations (MATLAB)
    // generate control gains
    // set control: u = -Kx
    // so... it's actually not that hard... we do NOT need to know the A,B,C,D matrices or define G and H for discrete time...
    // in fact, can just do the analysis in MATLAB, then port over the gains here, and confirm whether the system behavior is as expected

    // NOTE: we DO need to figure out what kind of forces are imparted by the wheels turning with Nm torque, since we are not 
    // providing a pure x-axis force on the cart axis.
    // r*f = tau
    // thus, if the "cart" is pushed with a force f, it is equivalent to motor torque over wheel radius: f = tau/r
    // since the MATLAB EOM consider the force applied to the wheel axis, we need to scale our control gains by r to
    // get the equivalent torque input from the motors on the cart axis

// x = [x, xdot, theta, thetadot]'
// u = -Kx
// K matrix (state feedback gain) (from MATLAB place function)
void mySSPcontroller(const mjModel* m, mjData* d)
{

    // since qpos is now quaternion, we need to convert to angle
    mjtNum euler[3];
    mjtNum quat[4] = {d->sensordata[0], d->sensordata[1], d->sensordata[2], d->sensordata[3]}; // using sensed angle of pendulum
    QuaternionToEuler(quat, euler);
    double xtheta = euler[0]; // use the roll (x-axis) angle as the pendulum angle

    static double last_update = 0.0; // last time the control was updated
    
    // guard check to update control at a fixed frequency
    if (d->time - last_update >= 1.0 / ctrl_update_freq)
    {

        // construct state vector
        mjtNum x[4];
        x[0] = d->sensordata[5]; // cart position; "x" is along the y-axis in the mujoco world frame
        x[1] = d->qvel[0]; // cart velocity TODO: estimate cart velocity using Luenberger observer with wheel position and gyro/accelerometer data
        x[2] = xtheta + 0.045; // - 0.036108;           // pendulum angle
        x[3] = d->qvel[3]; // pendulum x-axis angular velocity TODO: use sensordata gyro + accelerometer to get angular velocity
        // matrix multiplication
        double xd[4] = {1, 0, 0, 0}; // desired state vector
        double ctrl = -K_ssp[0] * (x[0] - xd[0])
                    - K_ssp[1] * (x[1] - xd[1])
                    - K_ssp[2] * (x[2] - xd[2])
                    - K_ssp[3] * (x[3] - xd[3]);

        d->ctrl[0] = ctrl; // apply control to the first actuator (left wheel)
        d->ctrl[1] = -ctrl; // apply control to the second actuator (right wheel)

        last_update = d->time;
        save_data(m,d);
    }

    // InjectControlNoise(); // inject noise into the control signal

}

// x = [x, xdot, theta, thetadot]'
// u = -K(x - x_desired)
// K matrix (state feedback gain) (from MATLAB LQR function)

void mySSPdesired_controller(const mjModel* m, mjData* d)
{

    // since qpos is now quaternion, we need to convert to angle
    mjtNum euler[3];
    mjtNum quat[4] = {d->sensordata[0], d->sensordata[1], d->sensordata[2], d->sensordata[3]}; // using sensed angle of pendulum
    QuaternionToEuler(quat, euler);
    double xtheta = euler[0]; // use the roll (x-axis) angle as the pendulum angle

    static double last_update = 0.0; // last time the control was updated
    
    // guard check to update control at a fixed frequency
    if (d->time - last_update >= 1.0 / ctrl_update_freq)
    {
        // cart velocity TODO: estimate cart velocity using Luenberger observer with wheel position and gyro/accelerometer data
        double cart_vel = sqrt((d->qvel[0]) * (d->qvel[0]) + (d->qvel[1]) * (d->qvel[1])); // magnitude of cart velocity in the plane of the ground

        // construct state vector
        mjtNum x[4];
        x[0] = 0; //d->sensordata[5]; // cart position; "x" is along the y-axis in the mujoco world frame
        x[1] = d->qvel[1]; // cart_vel; // TODO still need a way to determine direction of cart velocity
        x[2] = -xtheta - 0.038108;           // pendulum angle, offset so that CoM is over pivot point
        x[3] = -d->qvel[3]; // pendulum x-axis angular velocity TODO: use sensordata gyro + accelerometer to get angular velocity
        // matrix multiplication
        double xd[4] = {0, 0, 0, 0}; // desired state vector
        double ctrl = -K_ssp_desired[0] * (x[0] - xd[0])
                    - K_ssp_desired[1] * (x[1] - xd[1])
                    - K_ssp_desired[2] * (x[2] - xd[2])
                    - K_ssp_desired[3] * (x[3] - xd[3]);

        d->ctrl[0] = -ctrl; // apply control to the first actuator (left wheel)
        d->ctrl[1] = ctrl; // apply control to the second actuator (right wheel)

        last_update = d->time;
        save_data(m,d);
    }

    // InjectControlNoise(); // inject noise into the control signal

}

// SS PI controller - integral action added for position control of the cart
// x = [x, xdot, theta, thetadot, integral of position error]'
// u = -Kx
// K matrix (state feedback gain) (from MATLAB place function)
void mySSPIcontroller(const mjModel* m, mjData* d)
{

    // since qpos is now quaternion, we need to convert to angle
    mjtNum euler[3];
    mjtNum quat[4] = {d->sensordata[0], d->sensordata[1], d->sensordata[2], d->sensordata[3]}; // using sensed angle of pendulum
    QuaternionToEuler(quat, euler);
    double xtheta = euler[0]; // use the roll (x-axis) angle as the pendulum angle

    static double last_update = 0.0; // last time the control was updated
    static double error_integral = 0.0; // integral of the cart position error
    double yd = 0.0; // desired cart position
    static double int_lim = 5.0; // anti-windup limit for integral term

    // guard check to update control at a fixed frequency
    if (d->time - last_update >= 1.0 / ctrl_update_freq)
    {

        // construct state vector
        mjtNum x[4];
        x[0] = d->sensordata[5]; // cart position
        x[1] = d->qvel[1]; // cart velocity TODO: estimate cart velocity using Luenberger observer with wheel position and gyro/accelerometer data
        x[2] = -xtheta - 0.038108;          // pendulum angle + offset to account for CoM not being exactly over the pivot point
            // adding the offset keeps pendulum from running away, but control is very jittery, and seems to have
            // a sharp response during sign changes of cart roll (x-axis) angle
            // perhaps gain is too high? I should inspect what the bandwidth of my controller is and ensure
            // it's <= 25Hz
            // TODO: SS PI control continues to be unstable, need to investigate if possible to control in this way

        x[3] = -d->qvel[3]; // pendulum x-axis angular velocity TODO: use sensordata gyro + accelerometer to get angular velocity

        // limit the integral term to prevent windup
        error_integral = error_integral + yd - x[0] > int_lim ? int_lim : error_integral + yd - x[0];
        error_integral = error_integral < -int_lim ? -int_lim : error_integral; // approximate integral using summation
        
        // matrix multiplication
        // double ctrl = -K_sspi[0]*x[0] - K_sspi[1]*x[1] - K_sspi[2]*x[2] - K_sspi[3]*x[3] - K_sspi[4]*error_integral;
        double ctrl_P = -K_sspi[0]*x[0] - K_sspi[1]*x[1] - K_sspi[2]*x[2] - K_sspi[3]*x[3];
        double ctrl_I = -K_sspi[4]*error_integral;
        double ctrl = ctrl_P + ctrl_I;
        // double ctrl = ctrl_I > int_lim ? int_lim + ctrl_P : ctrl_I + ctrl_P; // saturate integral control to max 5Nm
        // ctrl = ctrl_I < -int_lim ? -int_lim + ctrl_P : ctrl; // saturate integral control to min -5Nm
        d->ctrl[0] = -ctrl; // apply control to the first actuator (left wheel)
        d->ctrl[1] = ctrl; // apply control to the second actuator (right wheel)

        last_update = d->time;
        save_data(m,d);
    }

    // InjectControlNoise(); // inject noise into the control signal

}