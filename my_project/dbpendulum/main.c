

#include<stdbool.h> //for bool
//#include<unistd.h> //for usleep
#include <math.h>

#include <mujoco/mujoco.h>
#include <GLFW/glfw3.h>
#include "stdio.h"
#include "stdlib.h"
#include "string.h"

//simulation end time
double simend = 5;

//related to writing data to a file
FILE *fid;
int loop_index = 0;
const int data_frequency = 100; //frequency at which data is written to a file


// char xmlpath[] = "../myproject/template_writeData/pendulum.xml";
// char datapath[] = "../myproject/template_writeData/data.csv";


//Change the path <template_writeData>
//Change the xml file
char path[] = ""; // "../myproject/dbpendulum/";
char xmlfile[] = "dbpendulum.xml";


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
float_t ctrl_update_freq = 100;
mjtNum last_update = 0.0;
mjtNum ctrl;

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


//****************************
//This function is called once and is used to get the headers
void init_save_data()
{
    //write name of the variable here (header)
    fprintf(fid,"t, ");
    fprintf(fid,"PE, KE, TE, ");
    fprintf(fid,"q1, q2, ");
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
    fprintf(fid,"%f, %f, %f, ",
        d->energy[0], d->energy[1], d->energy[0] + d->energy[1]);
    fprintf(fid,"%f, %f, ", d->qpos[0], d->qpos[1]);

    //Don't remove the newline
    fprintf(fid,"\n");
}

//**************************
// inject Brownian noise on the control signals
// really I should be injecting noise on the sensor signals, but this is a simple example
// and we want to see the effect of noise on the control signals
double ctrl_noise_std = 0.01; // standard deviation of noise to be injected into the control signal
double ctrl_noise_rate = 10; // rate of noise decay, in Hz

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

//**************************
void mycontroller(const mjModel* m, mjData* d)
{
    //write control here
    mj_energyPos(m, d);
    mj_energyVel(m, d);
    // printf("time: %f, potential: %f, kinetic: %f, total: %f\n",
    //     d->time, d->energy[0], d->energy[1], 
    //     d->energy[0] + d->energy[1]);

    //check equations that Mujoco has documented for the EOM
    //M*qacc + qfrc_bias = qfrc_applied + ctrl
    //M is the mass matrix
    //qfrc_bias is the bias force = Gravity + Coriolis + Centrifugal forces
    //qfrc_applied is a generalized force applied via the actuator
    //ctrl is the control signal applied to the actuator
    const int nv = 2; // number of actuators
    double dense_M[nv*nv];
    mj_fullM(m,dense_M,d->qM);
    double M[nv][nv];
    M[0][0] = dense_M[0];
    M[0][1] = dense_M[1];
    M[1][0] = dense_M[2];
    M[1][1] = dense_M[3];

    double qddot[nv];
    qddot[0] = d->qacc[0];
    qddot[1] = d->qacc[1];

    double f[nv];
    f[0] = d->qfrc_bias[0];
    f[1] = d->qfrc_bias[1];

    double lhs[nv];
    mju_mulMatVec(lhs,dense_M,qddot,nv,nv); // lhs = M*qddot
    lhs[0] = lhs[0] + f[0]; // lhs = M*qddot + qfrc_bias
    lhs[1] = lhs[1] + f[1];

    // actually applying the control signal now! input force
    // compensates for the Gravity and Coriolis forces
    // d->qfrc_applied[0] = f[0];
    // d->qfrc_applied[1] = f[1];

    double rhs[nv];
    rhs[0] = d->qfrc_applied[0];
    rhs[1] = d->qfrc_applied[1];

    // printf("%f %f \n",lhs[0],rhs[0]);
    // printf("%f %f \n",lhs[1],rhs[1]);
    // printf("******\n");

    // control
    double Kp0 = 100, Kp1 = 100;
    double Kv0 = 10, Kv1 = 10;
    double qref0 = -0.5, qref1 = -1.6; // desired position

    // PD control
    // d->qfrc_applied[0] = -Kp0*(d->qpos[0] - qref0)-Kv0*(d->qvel[0]);
    // d->qfrc_applied[1] = -Kp1*(d->qpos[1] - qref1)-Kv1*(d->qvel[1]);

    // coriolis + gravity + PD control
    // d->qfrc_applied[0] = f[0]-Kp0*(d->qpos[0] - qref0)-Kv0*(d->qvel[0]);
    // d->qfrc_applied[1] = f[1]-Kp1*(d->qpos[1] - qref1)-Kv1*(d->qvel[1]);

    // Feedback linearization control
    // M*(-kp( ...) -kv(...) + f)
    double tau[2];
    tau[0] = -Kp0*(d->qpos[0] - qref0)-Kv0*(d->qvel[0]) + f[0];
    tau[1] = -Kp1*(d->qpos[1] - qref1)-Kv1*(d->qvel[1]) + f[1];
    mju_mulMatVec(tau,dense_M,tau,nv,nv); // tau = M*tau
    tau[0] += f[0];
    tau[1] += f[1];
    d->qfrc_applied[0] = tau[0];
    d->qfrc_applied[1] = tau[1];


    // InjectNoise(); // inject noise into the control signal

    //write data here (dont change/dete this function call; instead write what you need to save in save_data)
    if ( loop_index%data_frequency==0){
        save_data(m,d);
    }
    loop_index = loop_index + 1;
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

    double arr_view[] = {90, -11.588379, 5, 0.000000, 0.000000, 2.000000};
    cam.azimuth = arr_view[0];
    cam.elevation = arr_view[1];
    cam.distance = arr_view[2];
    cam.lookat[0] = arr_view[3];
    cam.lookat[1] = arr_view[4];
    cam.lookat[2] = arr_view[5];

    // install control callback
    mjcb_control = mycontroller;

    fid = fopen(datapath,"w");
    init_save_data();

    d->qpos[0] = 0.5; // initial position
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
