

#include<stdbool.h> //for bool
//#include<unistd.h> //for usleep
#include <math.h>

#include <mujoco/mujoco.h>
#include <GLFW/glfw3.h>
#include "stdio.h"
#include "stdlib.h"
#include "string.h"


char filename[] = "pendulum.xml";

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
float_t ctrl_update_freq = 100;
mjtNum last_update = 0.0;
mjtNum ctrl;
double ctrl_noise_std = 0.01; // standard deviation of noise to be injected into the control signal
double ctrl_noise_rate = 10; // rate of noise decay, in Hz

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

// inject Brownian noise on the control signals
// really I should be injecting noise on the sensor signals, but this is a simple example
// and we want to see the effect of noise on the control signals
void InjectNoise() {
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

void set_position_servo(const mjModel*, int actuator_no, double kp) {
    // set position servo control signal
    if (actuator_no < 0 || actuator_no >= m->nu) {
        printf("Invalid actuator number: %d\n", actuator_no);
        return;
    }
    m->actuator_gainprm[10*actuator_no+0] = kp; // set the proportional gain
    m->actuator_biasprm[10*actuator_no+1] = -kp; 
}

void set_velocity_servo(const mjModel*, int actuator_no, double kv) {
    // set position servo control signal
    if (actuator_no < 0 || actuator_no >= m->nu) {
        printf("Invalid actuator number: %d\n", actuator_no);
        return;
    }
    m->actuator_gainprm[10*actuator_no+0] = kv; // set the proportional gain
    m->actuator_biasprm[10*actuator_no+2] = -kv; 
}

void set_torque_control(const mjModel*, int actuator_no, int flag) {
    // set torque control signal
    if (actuator_no < 0 || actuator_no >= m->nu) {
        printf("Invalid actuator number: %d\n", actuator_no);
        return;
    }
    if (flag == 0) {
        m->actuator_gainprm[10*actuator_no+0] = 0;
    } else {
        m->actuator_gainprm[10*actuator_no+0] = 1;
    }
}

void mycontroller(const mjModel* m, mjData* d)
{
    int actuator_no;
    //0 = torque actuator
    actuator_no = 0;
    int flag = 1; // actuator on/off
    set_torque_control(m, actuator_no, flag);
    d->ctrl[0] = -10*(d->qpos[0])-1*d->qvel[0]; // PD controller example
    // d->ctrl[0] = -10*(d->sensordata[0]-0)-1*d->sensordata[1]; // PD controller example
        // here we've zeroed the position and velocity servo gains

    //1 = position servo
    // // print the actuator gainprms
    actuator_no = 1;
    // for (int i = 0; i < 10; i++) {
    //     printf("%f \n", m->actuator_gainprm[10*actuator_no+i]);
    // }
    // printf("********** \n");
    double kp = 1;
    set_position_servo(m, actuator_no, kp); // set the position servo gain to kp
    d->ctrl[1] = 0.5; // set position to 0.5 radians
        // need the gain to be non-zero in the xml definition file

    //2 = velocity servo
    actuator_no = 2;
    double kv = 0;
    set_velocity_servo(m, actuator_no, kv); // set the velocity servo gain to kv
    d->ctrl[2] = 0.2; // set velocity to 0.2 radians per second
        // need the gain to be non-zero in the xml definition file

    // InjectNoise(); // inject noise into the control signal

    // proper PD controller example (easier to read)
    // get current time
    // mjtNum current_time = d->time;

    // // check if enough time has passed since last control update
    // if (current_time - last_update >= 1.0 / ctrl_update_freq) {
    //     // calculate the control signal (e.g., a simple PD controller)
    //     mjtNum position_error = 0 - d->qpos[0]; // target position is 0
    //     mjtNum velocity_error = 0 - d->qvel[0]; // target velocity is 0

    //     // PD control gains
    //     mjtNum kp = 10.0; // proportional gain
    //     mjtNum kd = 1.0;  // derivative gain

    //     // calculate control signal
    //     ctrl = kp * position_error + kd * velocity_error;

    //     // set control signal to actuator
    //     d->ctrl[0] = ctrl;

    //     // update last update time
    //     last_update = current_time;
    // }
}

// main function
int main(int argc, const char** argv)
{

    // load and compile model
    char error[1000] = "Could not load binary model";

    // check command-line arguments
    if( argc<2 )
        m = mj_loadXML(filename, 0, error, 1000);

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

    double arr_view[] = {-90, -15, 4, 0.000000, 0.000000, 0.000000};
    cam.azimuth = arr_view[0];
    cam.elevation = arr_view[1];
    cam.distance = arr_view[2];
    cam.lookat[0] = arr_view[3];
    cam.lookat[1] = arr_view[4];
    cam.lookat[2] = arr_view[5];

    d->qpos[0]=1.57; // pi/2, initial. qpos gives you the initial of the degree of freedom.
    // for a pendulum, this is the angle of the pendulum.

    mjcb_control = mycontroller; // set the controller callback

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

        // get framebuffer viewport
        mjrRect viewport = {0, 0, 0, 0};
        glfwGetFramebufferSize(window, &viewport.width, &viewport.height);

        // update scene and render
        mjv_updateScene(m, d, &opt, NULL, &cam, mjCAT_ALL, &scn);
        mjr_render(viewport, &scn, &con);
        //printf("{%f, %f, %f, %f, %f, %f};\n",cam.azimuth,cam.elevation, cam.distance,cam.lookat[0],cam.lookat[1],cam.lookat[2]);

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
