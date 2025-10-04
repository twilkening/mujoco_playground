#include "main.h"

/* define globals once here */
double simend = 40; 
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

// gains for state-space controller (using acker/place method)
// mjtNum K_ssp[4] = {-1.2957*r, -2.7642*r, -6.7969*r, -0.5284*r}; // scale by r to convert from force input to torque
mjtNum K_ssp[4] = {-0.0079*r, -1.6837*r, -3.6745*r, -0.1261*r};

// gains for state-space LQR controller with desired state
// mjtNum K_ssp_desired[4] = {-0.0088*r, -3.1188*r, -8.7115*r, -0.7069*r}; // Q = diag([1 1 10 100]); R = 1e4;
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
// mjtNum K_ssp_desired[4] = {-0.1945*r, -2.1261*r, -7.7281*r, -2.0436*r}; // Q = diag([1 1 10 100]); R = 0.1; // joint damping at 0.05 vs. 0.005
// mjtNum K_ssp_desired[4] = {-0.1953*r, -2.1137*r, -7.6722*r, -2.0492*r}; // Q = diag([1 1 10 100]); R = 1e-3; // joint damping at 0.029 vs. 0.05, ctrl/2
// updated MATLAB TWIP parameters 2025-09-27
// mjtNum K_ssp_desired[4] = {-0.1755*r, -0.6057*r, -7.1818*r, -1.8423*r}; // Q = diag([1 1 10 100]); R = 1e-3; 
// mjtNum K_ssp_desired[4] = {-2.9759*r, -2.8853*r, -14.7691*r, -2.0565*r}; // Q = diag([30 1 100 10]); R = 1e-3;
// mjtNum K_ssp_desired[4] = {-2.8984*r, -3.7110*r, -23.7741*r, -2.1660*r}; // Q = diag([30 1 1000 10]); R = 1e-3;
// mjtNum K_ssp_desired[4] = {-2.1149*r, -4.3963*r, -43.9063*r, -2.0583*r}; // Q = diag([30 1 1000 10]); R = 1;
// mjtNum K_ssp_desired[4] = {-0.6704*r, -2.5310*r, -43.1023*r, -1.9079*r}; // Q = diag([30 1 1e5 100]); R = 10;
// mjtNum K_ssp_desired[4] = {-1.1391*r, -2.4679*r, -25.4840*r, -1.3176*r}; // Q = diag([300 1 1e5 100]); R = 100;
mjtNum K_ssp_desired[4] = {-0.0550*r, -1.4363*r, -58.7283*r, -1.6316*r}; // Q = diag([1 300 1e6 100]); R = 100;


// gains for state-space controller with integral action
// mjtNum K_sspi[5] = {-1.9153*r, -3.7396*r, -9.4674*r, -0.7756*r, 0.0088*r}; // Q = diag([1 1 10 100 1]); R = 1e4;
// mjtNum K_sspi[5] = {-3.8108*r, -4.2684*r, -10.0934*r, -0.8324*r, 0.0276*r}; // Q = diag([1 1 10 100 10]); R = 1e4;
// mjtNum K_sspi[5] = {-7.9737*r, -5.2854*r, -11.2643*r, -0.9388*r, 0.0866*r}; // Q = diag([1 1 10 100 100]); R = 1e4;
// mjtNum K_sspi[5] = {-5.9055*r, -4.7983*r, -10.7082*r, -0.8883*r, 0.0550*r}; // Q = diag([1 1 10 100 40]); R = 1e4;
// mjtNum K_sspi[5] = {-28.1303*r, -10.0598*r, -18.4073*r, -1.7473*r, 0.4631*r}; // Q = diag([1 1 10 100 40]); R = 1e2;
// mjtNum K_sspi[5] = {-16.9013*r, -7.6822*r, -15.6549*r, -1.5270*r, 0.2337*r}; // Q = diag([1 1 10 100 10]); R = 1e2;
// mjtNum K_sspi[5] = {-7.4908*r, -5.3718*r, -12.7473*r, -1.3066*r, 0.0747*r}; // Q = diag([1 1 10 100 1]); R = 1e2;
// mjtNum K_sspi[5] = {-3.4935*r, -4.2012*r, -11.1118*r, -1.1909*r, 0.0238*r}; // Q = diag([1 1 10 100 0.1]); R = 1e2;
// mjtNum K_sspi[5] = {-3.4935*r, -4.2012*r, -11.1118*r, -1.1909*r, 0.0010*r}; // Q = diag([1 1 10 100 0.1]); R = 1e2; // modified Ki to reduce oscillations
// mjtNum K_sspi[5] = {-0.4785*r, -3.1247*r, -9.4194*r, -1.0802*r, 0.0008*r}; // Q = diag([1 1 10 100 1e-4]); R = 1e2;
// updated MATLAB TWIP parameters 2025-09-27
// mjtNum K_sspi[5] = {-6.6037*r, -3.2116*r, -11.3553*r, -1.1651*r, 0.0715*r}; // Q = diag([1 1 10 100 1]); R = 1e2; 
// mjtNum K_sspi[5] = {-52.7690*r, -25.7579*r, -67.9563*r, -3.5772*r, 0.5391*r}; // Q = diag([1 300 1e6 100 100]); R = 100;  // better with higher actuator damping
mjtNum K_sspi[5] = {-39.6016*r, -11.5254*r, -22.5286*r, -1.9423*r, 0.6802*r}; // Q = diag([1 300 1e4 100 100]); R = 100;
// mjtNum K_sspi[5] = {-3.0550*r, -1.9080*r, -8.5575*r, -0.6676*r, 0.0270*r}; // Q = diag([1 300 1e4 100 10]); R = 1e4; // still unstable, with damping = 0.0069


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

//**************************
double square_wave(double time, double period, double amplitude)
{
    if (period <= 0.0) {
        // invalid period: return 0 (could also assert or throw)
        return 0.0;
    }

    // fmod may return negative values for negative 'time', so normalize into [0, period)
    double phase_rem = std::fmod(time, period);
    if (phase_rem < 0.0)
        phase_rem += period;
    double phase = phase_rem / period; // now guaranteed in [0, 1)

    return (phase < 0.5) ? amplitude : -amplitude;
}

//**************************
struct RampFilter {
    double target = 0.0;
    double start = 0.0;
    double smoothed = 0.0;
    double t_start = -1e9;

    // call each control update; returns the current smoothed value
    double update(double time, double raw_target, double ramp_time, double eps) {
        if (ramp_time <= 0.0) {
            // immediate switch if no ramping requested
            smoothed = raw_target;
            target = raw_target;
            start = raw_target;
            t_start = time;
            return smoothed;
        }

        // detect meaningful change and start a new ramp
        if (std::fabs(raw_target - target) > eps) {
            target = raw_target;
            start = smoothed;
            t_start = time;
        }

        // compute interpolation factor in [0,1]
        double alpha = (time - t_start) / ramp_time;
        if (alpha <= 0.0) {
            return smoothed; // still at start
        }
        if (alpha >= 1.0) {
            smoothed = target;
            return smoothed;
        }

        smoothed = start + alpha * (target - start);
        return smoothed;
    }
};

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

double myPID(const mjModel* m, mjData* d,
            double r, double y,
            double kp, double ki, double kd,
            double int_limit=0.06,
            bool reset_integrator=false)
{
    static double e_prior = 0; // history of position error
    static double t_prior = 0; // history of time
    static double integral = 0; // integral of error

    // optional integrator reset
    if (reset_integrator) {
        integral = 0.0;
        kd = 0.0; // disable derivative action on reset
    }

    // calculate error
    double e = r - y;

    // calculate derivative of error
    // coarse one-step Euler method for derivative, may need to smooth out 
    // for better performance in the future

    // time step (guard against zero or negative dt)
    double dt = d->time - t_prior;
    if (dt <= 0.0) dt = 0.0;

    // derivative (safe: zero if dt == 0)
    double de = (dt > 0.0) ? (e - e_prior) / dt : 0.0;

    // update integral with anti-windup (clamp to +-int_limit if int_limit>0)
    integral += e * dt;
    if (int_limit > 0.0) {
        if (integral > int_limit) integral = int_limit;
        else if (integral < -int_limit) integral = -int_limit;
    }

    // calculate control signal
    double effort = kp * e + kd * de + ki * integral;

    // update history
    e_prior = e;
    t_prior = d->time;

    return effort;
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
    double integrator_limit = 0.5; // <-- integrator limit (tune as needed)

    static double last_update = 0.0; // last time the control was updated
    
    // guard check to update control at a fixed frequency
    if (d->time - last_update >= 1.0 / ctrl_update_freq)
    {
        // hard-coded PD controller example
        // d->ctrl[0] = -kp*(d->sensordata[0]-0)-kd*d->sensordata[1];
        // d->ctrl[1] = -d->ctrl[0];   // apply the inverse control to both wheels 

        // method-based PID controller
        double ctrl = myPID(m, d, r, xtheta, kp, ki, kd, integrator_limit);
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
        x[1] = d->qvel[1]; // cart velocity TODO: estimate cart velocity using Luenberger observer with wheel position and gyro/accelerometer data
        x[2] = -(xtheta + 0.036108); // - 0.036108;           // pendulum angle
        x[3] = -d->qvel[3]; // pendulum x-axis angular velocity TODO: use sensordata gyro + accelerometer to get angular velocity
        // matrix multiplication
        double xd[4] = {0, 0, 0, 0}; // desired state vector
        double ctrl = -K_ssp[0] * (x[0] - xd[0])
                    - K_ssp[1] * (x[1] - xd[1])
                    - K_ssp[2] * (x[2] - xd[2])
                    - K_ssp[3] * (x[3] - xd[3]);

        d->ctrl[0] = -ctrl; // apply control to the first actuator (left wheel)
        d->ctrl[1] = ctrl; // apply control to the second actuator (right wheel)

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
    static RampFilter offset_ramp;

    // guard check to update control at a fixed frequency
    if (d->time - last_update >= 1.0 / ctrl_update_freq)
    {
        // cart velocity TODO: estimate cart velocity using Luenberger observer with wheel position and gyro/accelerometer data
        // double cart_vel = sqrt((d->qvel[0]) * (d->qvel[0]) + (d->qvel[1]) * (d->qvel[1])); // magnitude of cart velocity in the plane of the ground

        // construct state vector
        mjtNum x[4];
        x[0] = d->sensordata[5]; // cart position; "x" is along the y-axis in the mujoco world frame
        x[1] = d->qvel[1]; // cart_vel; // TODO still need a way to determine direction of cart velocity
        x[2] = -(xtheta + 0.058508);           // pendulum angle, offset so that CoM is over pivot point
        x[3] = -d->qvel[3]; // pendulum x-axis angular velocity TODO: use sensordata gyro + accelerometer to get angular velocity
        // matrix multiplication

        double xd0 = square_wave(d->time, 40.0, 0.2);
        double integrator_limit = 0.02; // anti-windup limit for PID position controller        double raw_offset = myPID(m, d, xd[0], x[0], 0.6, 0.3, 1, integrator_limit); // use PID to determine desired angle offset based on position error
        
        static double prev_xd0 = 0.0;
        const double integral_eps = 1e-6;
        bool reset_integrator = (std::fabs(xd0 - prev_xd0) > integral_eps);
        prev_xd0 = xd0;
        double xd[4] = {xd0, 0, 0, 0}; // desired state vector
        double raw_offset = myPID(m, d, xd[0], x[0], 0.3, 0, 0.8, integrator_limit, reset_integrator);

        double ramp_time = 12.0 / ctrl_update_freq; // 12 controller periods -> 0.1s when ctrl_update_freq==100
        // threshold to detect meaningful changes in offset command, to start a new ramp. otherwise will update offset every control step
        double eps = 0.1; 
        double offset_theta = offset_ramp.update(d->time, raw_offset, ramp_time, eps);
        double ctrl = -K_ssp_desired[0] * (x[0] - xd[0])
                    - K_ssp_desired[1] * (x[1] - xd[1])
                    - K_ssp_desired[2] * (x[2] - xd[2] - offset_theta) // add position control offset to desired angle
                    - K_ssp_desired[3] * (x[3] - xd[3]);
        std::cout << "control: " << ctrl << ", cart pos: " << x[0] << ", desired pos: " << xd[0] << ", offset_theta: " << offset_theta << std::endl;

        d->ctrl[0] = -ctrl/2; // apply control to the first actuator (left wheel, when looking from battery side)
        d->ctrl[1] = ctrl/2; // apply control to the second actuator (right wheel)

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
    double yd = 0; // desired cart position
    static double int_lim = 2.0; // anti-windup limit for integral term

    // guard check to update control at a fixed frequency
    if (d->time - last_update >= 1.0 / ctrl_update_freq)
    {

        // construct state vector
        mjtNum x[4];
        x[0] = d->sensordata[5]; // cart position
        x[1] = d->qvel[1]; // cart velocity TODO: estimate cart velocity using Luenberger observer with wheel position and gyro/accelerometer data
        x[2] = -(xtheta + 0.036108);          // pendulum angle + offset to account for CoM not being exactly over the pivot point
            // adding the offset keeps pendulum from running away, but control is very jittery, and seems to have
            // a sharp response during sign changes of cart roll (x-axis) angle
            // perhaps gain is too high? I should inspect what the bandwidth of my controller is and ensure
            // it's <= 25Hz
            // TODO: SS PI control continues to be unstable, need to investigate if possible to control in this way

        x[3] = -d->qvel[3]; // pendulum x-axis angular velocity TODO: use sensordata gyro + accelerometer to get angular velocity

        // limit the integral term to prevent windup
        error_integral = error_integral + yd - x[0] > int_lim ? int_lim : error_integral + yd - x[0];
        error_integral = error_integral < -int_lim ? -int_lim : error_integral; // approximate integral using summation
        std::cout << "Integral of position error: " << error_integral << ", Actual Position: " << x[0] << std::endl;

        // matrix multiplication
        // double ctrl = -K_sspi[0]*x[0] - K_sspi[1]*x[1] - K_sspi[2]*x[2] - K_sspi[3]*x[3] - K_sspi[4]*error_integral;
        double ctrl_P = -K_sspi[0]*x[0] - K_sspi[1]*x[1] - K_sspi[2]*x[2] - K_sspi[3]*x[3];
        // double ctrl_P = -K_sspi[0]*(x[0] - yd) - K_sspi[1]*x[1] - K_sspi[2]*x[2] - K_sspi[3]*x[3]; // offset position feedback to desired position
        // double ctrl_P = -K_sspi[0]*0 - K_sspi[1]*x[1] - K_sspi[2]*x[2] - K_sspi[3]*x[3]; // zero out position feedback to test integral action
        double ctrl_I = -K_sspi[4]*error_integral;
        double ctrl = ctrl_P + ctrl_I;
        // double ctrl = ctrl_I > int_lim ? int_lim + ctrl_P : ctrl_I + ctrl_P; // saturate integral control to max 5Nm
        // ctrl = ctrl_I < -int_lim ? -int_lim + ctrl_P : ctrl; // saturate integral control to min -5Nm
        d->ctrl[0] = -ctrl/2; // apply control to the first actuator (left wheel)
        d->ctrl[1] = ctrl/2; // apply control to the second actuator (right wheel)

        last_update = d->time;
        save_data(m,d);
    }

    // InjectControlNoise(); // inject noise into the control signal

}