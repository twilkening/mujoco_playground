#ifndef CART_PEND_MAIN_H
#define CART_PEND_MAIN_H

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

/* globals (define them once in control.cpp) */
extern double simend;
extern FILE* fid;
extern int loop_index;
extern const int data_frequency;

extern char path[];
extern char xmlfile[];
extern char datafile[];

/* controllers/history */
extern mjtNum position_history;
extern mjtNum previous_time;
extern double ctrl_update_freq;
extern double ctrl_noise_std;
extern double ctrl_noise_rate;

/* controller gain arrays and constants */
extern const double r;
extern mjtNum K_ssp[4];
extern mjtNum K_ssp_desired[4];
extern mjtNum K_sspi[5];

/* function prototypes */
void QuaternionToEuler(const mjtNum* quat, mjtNum* euler);

void init_save_data(void);
void save_data(const mjModel* m, mjData* d);

void InjectControlNoise(void);

double myPID(const mjModel* m, mjData* d, double r, double y, double kp, double ki, double kd);

void myPIDcontroller(const mjModel* m, mjData* d);
void mySSPcontroller(const mjModel* m, mjData* d);
void mySSPdesired_controller(const mjModel* m, mjData* d);
void mySSPIcontroller(const mjModel* m, mjData* d);

#endif // CART_PEND_MAIN_H
