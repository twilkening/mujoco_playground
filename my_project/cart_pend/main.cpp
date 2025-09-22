
#define _USE_MATH_DEFINES
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cstdbool> //for bool
//#include <unistd.h> //for usleep
#include <cmath>
#include <iostream> // For C++ style I/O

#include "mujoco/mujoco.h"
#include <GLFW/glfw3.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif




//simulation end time
double simend = 20;

//related to writing data to a file
FILE *fid;
int loop_index = 0;
const int data_frequency = 50; // #'d loop iteration at which data is written to a file


//Change the path <template_writeData>
//Change the xml file
char path[] = ""; // "../myproject/dbpendulum/";
char xmlfile[] = "cart_pend_meshes.xml";


char datafile[] = "data.csv";


// MuJoCo data structures
mjModel* m = NULL;                  // MuJoCo model
mjData* d = NULL;                   // MuJoCo data
mjvCamera cam;                      // abstract camera
mjvOption opt;                      // visualization options
mjvScene scn;                       // abstract scene
mjrContext con;                     // custom GPU context

// mouse interaction
bool button_left = false;
bool button_middle = false;
bool button_right =  false;
double lastx = 0;
double lasty = 0;

// holders of one step history of time and position to calculate dertivatives
mjtNum position_history = 0;
mjtNum previous_time = 0;

// controller related variables
double ctrl_update_freq = 100.0; // Hz

// keyboard callback
void keyboard(GLFWwindow* window, int key, int scancode, int act, int mods)
{
    // backspace: reset simulation
    if( act==GLFW_PRESS && key==GLFW_KEY_BACKSPACE )
    {
        mj_resetData(m, d);
        mj_forward(m, d);
    }
}

// mouse button callback
void mouse_button(GLFWwindow* window, int button, int act, int mods)
{
    // update button state
    button_left =   (glfwGetMouseButton(window, GLFW_MOUSE_BUTTON_LEFT)==GLFW_PRESS);
    button_middle = (glfwGetMouseButton(window, GLFW_MOUSE_BUTTON_MIDDLE)==GLFW_PRESS);
    button_right =  (glfwGetMouseButton(window, GLFW_MOUSE_BUTTON_RIGHT)==GLFW_PRESS);

    // update mouse position
    glfwGetCursorPos(window, &lastx, &lasty);
}


// mouse move callback
void mouse_move(GLFWwindow* window, double xpos, double ypos)
{
    // no buttons down: nothing to do
    if( !button_left && !button_middle && !button_right )
        return;

    // compute mouse displacement, save
    double dx = xpos - lastx;
    double dy = ypos - lasty;
    lastx = xpos;
    lasty = ypos;

    // get current window size
    int width, height;
    glfwGetWindowSize(window, &width, &height);

    // get shift key state
    bool mod_shift = (glfwGetKey(window, GLFW_KEY_LEFT_SHIFT)==GLFW_PRESS ||
                      glfwGetKey(window, GLFW_KEY_RIGHT_SHIFT)==GLFW_PRESS);

    // determine action based on mouse button
    mjtMouse action;
    if( button_right )
        action = mod_shift ? mjMOUSE_MOVE_H : mjMOUSE_MOVE_V;
    else if( button_left )
        action = mod_shift ? mjMOUSE_ROTATE_H : mjMOUSE_ROTATE_V;
    else
        action = mjMOUSE_ZOOM;

    // move camera
    mjv_moveCamera(m, action, dx/height, dy/height, &scn, &cam);
}


// scroll callback
void scroll(GLFWwindow* window, double xoffset, double yoffset)
{
    // emulate vertical mouse motion = 5% of window height
    mjv_moveCamera(m, mjMOUSE_ZOOM, 0, -0.05*yoffset, &scn, &cam);
}

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
    fprintf(fid,"RightControl, LeftControl,");
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
double ctrl_noise_std = 0.01; // standard deviation of noise to be injected into the control signal
double ctrl_noise_rate = 50; // rate of noise decay, in Hz

void InjectControlNoise() {
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
const double r = 0.04; // radius of the wheels, in meters
// mjtNum K_ssp[4] = {-1.2957*r, -2.7642*r, -6.7969*r, -0.5284*r}; // scale by r to convert from force input to torque
mjtNum K_ssp[4] = {-0.0079*r, -1.6837*r, -3.6745*r, -0.1261*r}; // scale by r to convert from force input to torque
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
// scale by r to convert from force input to torque
// mjtNum K_ssp_desired[4] = {-0.1952*r, -2.1149*r, -7.6775*r, -2.0487*r}; // 
// mjtNum K_ssp_desired[4] = {-0.0273*r, -3.0981*r, -8.7937*r, -0.7602*r}; // 
// mjtNum K_ssp_desired[4] = {-0.0279*r, -3.1258*r, -8.7201*r, -0.7077*r}; // 
// mjtNum K_ssp_desired[4] = {-0.0760*r, -2.9516*r, -9.1199*r, -1.0618*r}; // Q = diag([1 1 10 100]); R = 1e2;
mjtNum K_ssp_desired[4] = {-0.0881*r, -3.1477*r, -8.7474*r, -0.7102*r}; // Q = diag([100 1 10 100]); R = 1e4;
// mjtNum K_ssp_desired[4] = {-1.8645*r, -3.4969*r, -11.9430*r, -2.1176*r}; // Q = diag([100 1 10 100]); R = 1;
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

        // construct state vector
        mjtNum x[4];
        x[0] = d->sensordata[5]; // cart position; "x" is along the y-axis in the mujoco world frame
        x[1] = d->qvel[0]; // cart velocity TODO: estimate cart velocity using Luenberger observer with wheel position and gyro/accelerometer data
        x[2] = xtheta + 0.06;           // pendulum angle, offset so that CoM is over pivot point
        x[3] = d->qvel[3]; // pendulum x-axis angular velocity TODO: use sensordata gyro + accelerometer to get angular velocity
        // matrix multiplication
        double xd[4] = {0, 0, 0, 0}; // desired state vector
        double ctrl = -K_ssp_desired[0] * (x[0] - xd[0])
                    - K_ssp_desired[1] * (x[1] - xd[1])
                    - K_ssp_desired[2] * (x[2] - xd[2])
                    - K_ssp_desired[3] * (x[3] - xd[3]);

        d->ctrl[0] = ctrl; // apply control to the first actuator (left wheel)
        d->ctrl[1] = -ctrl; // apply control to the second actuator (right wheel)

        last_update = d->time;
        save_data(m,d);
    }

    // InjectControlNoise(); // inject noise into the control signal

}

// SS PI controller - integral action added for position control of the cart
// x = [x, xdot, theta, thetadot, integral of position error]'
// u = -Kx
// K matrix (state feedback gain) (from MATLAB place function)
// mjtNum K_sspi[5] = {-10.8317*r, -5.7073*r, -11.0581*r, -0.9188*r, 0.1115*r}; // scale by r to convert from force input to torque
// mjtNum K_sspi[5] = {-19.5492*r, -8.0521*r, -13.6899*r, -1.1253*r, 0.2149*r}; // scale by r to convert from force input to torque
mjtNum K_sspi[5] = {-15.8509*r, -8.4484*r, -20.6371*r, -2.5579*r, 0.1828*r};
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
        x[0] = d->sensordata[4]; // cart position
        x[1] = d->qvel[0]; // cart velocity TODO: estimate cart velocity using Luenberger observer with wheel position and gyro/accelerometer data
        x[2] = xtheta - 0.036108;           // pendulum angle + offset to account for CoM not being exactly over the pivot point
            // adding the offset keeps pendulum from running away, but control is very jittery, and seems to have
            // a sharp response during sign changes of cart roll (x-axis) angle
            // perhaps gain is too high? I should inspect what the bandwidth of my controller is and ensure
            // it's <= 25Hz
            // TODO: SS PI control continues to be unstable, need to investigate if possible to control in this way

        x[3] = d->qvel[3]; // pendulum x-axis angular velocity TODO: use sensordata gyro + accelerometer to get angular velocity
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
        d->ctrl[0] = ctrl; // apply control to the first actuator (left wheel)
        d->ctrl[1] = -ctrl; // apply control to the second actuator (right wheel)

        last_update = d->time;
        save_data(m,d);
    }

    // InjectControlNoise(); // inject noise into the control signal

}

//************************
// main function
int main(int argc, const char** argv)
{

    char xmlpath[100]={};
    char datapath[100]={};

    strcat(xmlpath,path);
    strcat(xmlpath,xmlfile);

    strcat(datapath,path);
    strcat(datapath,datafile);


    // load and compile model
    char error[1000] = "Could not load binary model";

    // check command-line arguments
    if( argc<2 )
        m = mj_loadXML(xmlpath, 0, error, 1000);

    else
        if( strlen(argv[1])>4 && !strcmp(argv[1]+strlen(argv[1])-4, ".mjb") )
            m = mj_loadModel(argv[1], 0);
        else
            m = mj_loadXML(argv[1], 0, error, 1000);
    if( !m )
        mju_error_s("Load model error: %s", error);

    // make data
    d = mj_makeData(m);


    // init GLFW
    if( !glfwInit() )
        mju_error("Could not initialize GLFW");

    // create window, make OpenGL context current, request v-sync
    GLFWwindow* window = glfwCreateWindow(1244, 700, "Demo", NULL, NULL);
    glfwMakeContextCurrent(window);
    glfwSwapInterval(1);

    // initialize visualization data structures
    mjv_defaultCamera(&cam);
    mjv_defaultOption(&opt);
    mjv_defaultScene(&scn);
    mjr_defaultContext(&con);
    mjv_makeScene(m, &scn, 2000);                // space for 2000 objects
    mjr_makeContext(m, &con, mjFONTSCALE_150);   // model-specific context

    // install GLFW mouse and keyboard callbacks
    glfwSetKeyCallback(window, keyboard);
    glfwSetCursorPosCallback(window, mouse_move);
    glfwSetMouseButtonCallback(window, mouse_button);
    glfwSetScrollCallback(window, scroll);

    double arr_view[] = {161.485714, -26.314286, 1.880372, 0.000000, 0.000000, 0.000000};
    cam.azimuth = arr_view[0];
    cam.elevation = arr_view[1];
    cam.distance = arr_view[2];
    cam.lookat[0] = arr_view[3];
    cam.lookat[1] = arr_view[4];
    cam.lookat[2] = arr_view[5];

    // install control callback
    // mjcb_control = myPIDcontroller;
    // mjcb_control = mySSPcontroller;
    mjcb_control = mySSPdesired_controller;
    // mjcb_control = mySSPIcontroller;
    
    fid = fopen(datapath,"w");
    init_save_data();

    // print the number of values in the sensor data; get insight into if gyro and acclerometer return 3 values each
    printf("number of sensor values = %d\n", m->nsensordata);
    // get the # of generalized coordinates, DOF of model, and # of actuators
    printf("# of generalized coordinates: nq = %d,\n# of DOF of model: nv = %d,\n# of actuators: nu = %d\n", m->nq, m->nv, m->nu);

    // set the initial angle of pendulum
    // mjtNum euler_init[3];
    // mjtNum quat_init[4];
    // euler_init[0] = - 15.0 * M_PI / 180.0; // 15 degrees from vertical
    // euler_init[1] = 0.0;
    // euler_init[2] = 0.0;
    // // convert to quaternion
    // quat_init[0] = cos(euler_init[0]/2) * cos(euler_init[1]/2) * cos(euler_init[2]/2) + sin(euler_init[0]/2) * sin(euler_init[1]/2) * sin(euler_init[2]/2);
    // quat_init[1] = sin(euler_init[0]/2) * cos(euler_init[1]/2) * cos(euler_init[2]/2) - cos(euler_init[0]/2) * sin(euler_init[1]/2) * sin(euler_init[2]/2);
    // quat_init[2] = cos(euler_init[0]/2) * sin(euler_init[1]/2) * cos(euler_init[2]/2) + sin(euler_init[0]/2) * cos(euler_init[1]/2) * sin(euler_init[2]/2);
    // quat_init[3] = cos(euler_init[0]/2) * cos(euler_init[1]/2) * sin(euler_init[2]/2) - sin(euler_init[0]/2) * sin(euler_init[1]/2) * cos(euler_init[2]/2);
    // d->qpos[3] = quat_init[0];
    // d->qpos[4] = quat_init[1];
    // d->qpos[5] = quat_init[2];
    // d->qpos[6] = quat_init[3];
    // mj_forward(m, d); // call forward to update the model with the new initial position
    
    // use the first while condition if you want to simulate for a period.
    while( !glfwWindowShouldClose(window))
    {
        // advance interactive simulation for 1/60 sec
        //  Assuming MuJoCo can simulate faster than real-time, which it usually can,
        //  this loop will finish on time for the next frame to be rendered at 60 fps.
        //  Otherwise add a cpu timer and exit this loop when it is time to render.
        mjtNum simstart = d->time;
        while( d->time - simstart < 1.0/60.0 )
        {
            mj_step(m, d);
        }

        if (d->time>=simend)
        {
            fclose(fid);
            break;
        }

       // get framebuffer viewport
        mjrRect viewport = {0, 0, 0, 0};
        glfwGetFramebufferSize(window, &viewport.width, &viewport.height);

        // update scene and render
        // opt.frame = mjFRAME_BODY; // show body frame
        mjv_updateScene(m, d, &opt, NULL, &cam, mjCAT_ALL, &scn);
        mjr_render(viewport, &scn, &con);
        // printf("{%f, %f, %f, %f, %f, %f};\n",cam.azimuth,cam.elevation, cam.distance,cam.lookat[0],cam.lookat[1],cam.lookat[2]);

        // swap OpenGL buffers (blocking call due to v-sync)
        glfwSwapBuffers(window);

        // process pending GUI events, call GLFW callbacks
        glfwPollEvents();

    }

    // free visualization storage
    mjv_freeScene(&scn);
    mjr_freeContext(&con);

    // free MuJoCo model and data, deactivate
    mj_deleteData(d);
    mj_deleteModel(m);

    // terminate GLFW (crashes with Linux NVidia drivers)
    #if defined(__APPLE__) || defined(_WIN32)
        glfwTerminate();
    #endif

    return 1;
}
