#include <stdio.h>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <cmath>
#include <algorithm>
//-------------------------------------------------------------------------
//       Program Overview
//-------------------------------------------------------------------------
//       The purpose of this program is to calculate the altitude and
//       velocity of a rocket at a given time t, taking into account
//       air resistance, changing air density, changing gravity,
//       the thrust of the rocket, and the changing mass of the rocket.
//       It would be impossible to determine a closed solution for
//       this problem, so using numerical methods is the only possible
//       way to find the solution. In this program, the Runge-Kutta-Fehlberg
//       method is used for integration (4th and 5th order). 
//-------------------------------------------------------------------------
//       Variables
//-------------------------------------------------------------------------
//       Earthradius [m] Radius of Earth
//       Initgravity [m/s²] Gravity at surface
//       Initdensity [kg/m³] Air density at surface
//       R [Joules/kg*K] Gas constant for air
//       Temp [K] Mean temperature of atmosphere
//       Endtime [s] When to stop plotting
//       Rocketmass [kg] Mass of rocket without fuel
//       Fuelmass [kg] Mass of fuel
//       Impulse [s] Impulse of engine
//                   300 - Kerosene/Oxygen
//                   360 - Hydrogen/Oxygen
//                   490 - Hydrogen/Fluorine
//       Burnrate [kg/s] Amount of fuel burned per second
//       Endthrust [s] Total time of thrust
//       Velocity [m/s] Velocity of material ejected from rocket nozzle
//       Initthrust [N] Thrust of engine (constant)
//       Cd [dimensionless] Coefficient of Drag
//       SurfaceArea [m²] Frontal surface area
//       Initheight [m] Initial altitude of rocket
//       Initvelocity [m/s] Initial velocity of rocket
//       time [s] Current time
//       dt [s] Step size
//       Vold,Vnew [m/s] Velocity
//       Hold, Hnew [m] Altitude
//       Used in the Runge-Kutta-Fehlberg (45) integration method
//       function accel computes the overall acceleration of the rocket at
//       a given time, altitude, and velocity.
//-------------------------------------------------------------------------   
//-------------------------------------------------------------------------
//       Initialization (Constant Physics)
//------------------------------------------------------------------------- 
const double R_univ = 8.3144621;
const double g = 9.80665;
const double T_stp = 273.15;
const double Ts = 216.65;
const double Tamb = 288.15;   
const double Pamb = 101325; 
const double M = 0.0289645;  
const double L = 0.0065;   
const double R_earth = 6378100.0; 
const double rho_air = 1.2754;  
const double gamma_air = 1.4;  
const double pi = 3.14159265359;
//-------------------------------------------------------------------------
//       Initial values (Variables)
//------------------------------------------------------------------------- 
double  Rocketmass = 500.0,
        Fuelmass = 300.0,
        Impulse = 300.0,                        // Specific impulse [s]
        Burnrate = 5.0,                         // Mass flow of propellant [kg/s]
        Initgravity = g,                
        Earthradius = R_earth,    
        Endthrust = Fuelmass/Burnrate,
        Velocity = Initgravity*Impulse,         // Effective exhaust velocity [m/s] 
        Initthrust = Burnrate*Velocity,         // Thrust [N]
        SurfaceArea = 0.1,
        NozzleArea = 0.025,
        NozzlePressure = 0.9*Pamb,
        Initheight = 0.0,
        Initvelocity = 0.0,
        Endtime = 1e100;

//-------------------------------------------------------------------------
//       Functions
//-------------------------------------------------------------------------
double Thrust(double t){
    // This function computes the thrust at a given time. Note that
    // a more complex thrust curve could be put here rather than just
    // a constant value.
    double Thrust = Initthrust;
    if (t >= Endthrust){
        Thrust = 0;
    }

    return Thrust;
}
//-------------------------------------------------------------------------
double Mass(double t){
    // This function computes the mass of the rocket at a given 
    // time t.
    double Mass = Rocketmass;
    if (t < Endthrust){
        Mass = Rocketmass + Fuelmass - Burnrate*t;
    }
        
    return Mass;    
}
//-------------------------------------------------------------------------
double Gravity(double h){
    // This function computes the acceleration of gravity at an
    // altitude h, given that the gravity at the surface is
    // Initgravity.
    double Gravity; 
    Gravity = Initgravity*pow((Earthradius/(Earthradius + h)),2);

    return Gravity;
}
//--------------------------------------------------------------------------
double Density(double h){       
    // Function Density finds the density at a given altitude h.
    // Taken from the Barometric formula. Accurate up to 80,000m above sea level
    // The altitude model used below is based on the standard atmospheric model used in modern meteorology.
    // It takes into accout the different regression rates and properties of the thermoclines.    
    static double pb[7] = {1.2250, 0.36391, 0.08803, 0.01322, 0.00143, 0.00086, 0.000064};
    static double hb[7] = {0.0, 11000, 20000, 32000, 47000, 51000, 71000};
    static double Tb[7] = {288.15, 216.65, 216.65, 228.65, 270.65, 270.65, 214.65};
    static double Lb[7] = {-0.0065, 0.0, 0.001, 0.0028, 0.0, -0.0028, -0.002};
    double Density;
    int i;

    if (h < 11000){
        i = 0;
        Density = pb[i]*pow((Tb[i]/(Tb[i] + Lb[i]*(h - hb[i]))),(1 + (Initgravity*M)/(R_univ*Lb[i])));
    
        return Density;
    }

    else if (h > 11000 and h < 20000){
        i = 1;
        Density = pb[i]*exp((-1*(Initgravity)*M*(h - hb[i]))/(R_univ*Tb[i]));
        
        return Density;
    }
    
    else if (h > 20000 and h < 32000){
        i = 2;
        Density = pb[i]*pow((Tb[i]/(Tb[i] + Lb[i]*(h - hb[i]))),(1 + (Initgravity*M)/(R_univ*Lb[i])));
        
        return Density;
    }

    else if (h > 32000 and h < 47000){
        i = 3;
        Density = pb[i]*pow((Tb[i]/(Tb[i] + Lb[i]*(h - hb[i]))),(1 + (Initgravity*M)/(R_univ*Lb[i])));
    
        return Density;
    }
    else if (h > 47000 and h < 51000){
        i = 4;
        Density = pb[i]*exp((-1 * (Initgravity)*M*(h - hb[i]))/(R_univ*Tb[i]));
        
        return Density;
    }
    
    else if (h > 51000 and h < 71000){
        i = 5;
        Density = pb[i]*pow((Tb[i]/(Tb[i] + Lb[i]*(h - hb[i]))),(1 + (Initgravity*M)/(R_univ*Lb[i])));

        return Density;
    }

    else if (h > 71000 and h < 86000){
        i = 6;
        Density = pb[i]*pow((Tb[i]/(Tb[i] + Lb[i]*(h - hb[i]))),(1 + (Initgravity*M)/(R_univ*Lb[i])));

        return Density;
    }
    
    else{
        return 0.0;
    }
    return Density;
}

//-----------------------------------------------------------------------------------------------------
double Temperature(double h){
    // Calculates air temperature [Kelvin] at altitude [m]
    // from equations at 
    // http://www.grc.nasa.gov/WWW/K-12/airplane/atmosmet.html
    double T;
    if (h <= 11000){
        // Troposphere
        T = Tamb - L*h;
    }
    else if (h <= 25000){
        // Lower Stratosphere
        T = 216.69;
    }
    else{
        T = 141.94 + 0.00299*h;
    }
    return T;
}
//------------------------------------------------------------------------------------------------------
double Pressure(double h){
    // Calculates air pressure [Pa] at altitude [m]"
    // from equations at 
    // http://www.grc.nasa.gov/WWW/K-12/airplane/atmosmet.html
    double p;
    p = Density(h)*(R_univ/M)*Temperature(h);

    return p;
}
//-------------------------------------------------------------------------------------------------------
double Speed_sound(double h){
    // Calculates seeped sound [m/s] at altitude [m]"
    double speed_sound;
    speed_sound = sqrt(gamma_air*(R_univ/M)*Temperature(h));
        
    return speed_sound;
}
//-------------------------------------------------------------------------------------------------------
double Mach_number(double v, double h){
    // Calculates Mach Number
    double Mach;
    Mach = v/Speed_sound(h);
    
    return Mach;
}
//-------------------------------------------------------------------------------------------------------
double to_radians(double degree){
    // Converts degrees to radians
    double radians = degree*pi/180;

    return radians;
}
//--------------------------------------------------------------------------------------------------------
double Cdrag(double v, double h){
    // Calculates drag coefficient [-] at altitude [m] and Mach number[-]"
    // Drag function for V2
    // derived from Sutton, "Rocket Propulsion Elements", 7th ed, p108
    // probably not that relevant to other body types
    // use nose cone formula
    double  theta = to_radians(15);  
    double  c = Speed_sound(h); 
    double  mach = Mach_number(v,c);
    double  cd;
    
    if (mach > 5){
        cd = 2*(pow(sin(theta),2));
    }         
    else if (mach > 1.8 && mach <= 5){
        cd = -0.03125*mach + 0.30625;
    }
    else if (mach > 1.2 && mach <= 1.8){
        cd = -0.25*mach + 0.7;
    }
    else if (mach > 0.8 && mach <= 1.2){
        cd = 0.625*mach - 0.35;
    }
    else if (mach <= 0.8){
        cd = 0.15;
    }   
    return cd;
}
//--------------------------------------------------------------------------
double  Drag(double v,double h){
    // Function Drag computes the air drag on the rocket for a certain
    // velocity and altitude. Note that drag acts in the direction 
    // opposite of the given velocity v. 
    double Drag;
    Drag = Cdrag(v,h)*(0.5)*Density(h)*(pow(v,2))*SurfaceArea;

    return Drag;
}
//-------------------------------------------------------------------------------
double accel(double h, double v, double t){
    // This function computes the total acceleration of the 
    // rocket at a given time, height, and velocity.
    double Acceleration;
    Acceleration = (Thrust(t) - Drag(v,h))/Mass(t) - Gravity(h);

    return Acceleration;
}
//---------------------------------------------------------------------------------
int main(void){
    //-----------------------------------------------------------------------------
    //       Main Program
    //-----------------------------------------------------------------------------
    //       Saves the data in .txt
    //-----------------------------------------------------------------------------
    FILE *file0;
    FILE *file1;
    FILE *file2;
    FILE *file3;
    FILE *file4;

    file0 = fopen("Time.txt", "w");
    file1 = fopen("Velocity.txt", "w");
    file2 = fopen("Altitude.txt", "w");
    file3 = fopen("Relative_error.txt", "w");
    file4 = fopen("Step.txt", "w");
    
    // Coefficients used to compute the independent variable argument of f [time]
    double  a2  =   2.500000000000000e-01;  //  1/4
    double  a3  =   3.750000000000000e-01;  //  3/8
    double  a4  =   9.230769230769231e-01;  //  12/13
    double  a5  =   1.000000000000000e+00;  //  1
    double  a6  =   5.000000000000000e-01;  //  1/2

    // Coefficients used to compute the dependent variable argument of f [velocity and acceleration]
    double  b21 =   2.500000000000000e-01;  //  1/4
    double  b31 =   9.375000000000000e-02;  //  3/32
    double  b32 =   2.812500000000000e-01;  //  9/32
    double  b41 =   8.793809740555303e-01;  //  1932/2197
    double  b42 =  -3.277196176604461e+00;  // -7200/2197
    double  b43 =   3.320892125625853e+00;  //  7296/2197
    double  b51 =   2.032407407407407e+00;  //  439/216
    double  b52 =  -8.000000000000000e+00;  // -8
    double  b53 =   7.173489278752436e+00;  //  3680/513
    double  b54 =  -2.058966861598441e-01;  // -845/4104
    double  b61 =  -2.962962962962963e-01;  // -8/27
    double  b62 =   2.000000000000000e+00;  //  2
    double  b63 =  -1.381676413255361e+00;  // -3544/2565
    double  b64 =   4.529727095516569e-01;  //  1859/4104
    double  b65 =  -2.750000000000000e-01;  // -11/40

    // Coefficients used to compute local truncation error estimate.  These
    // come from subtracting a 4th order RK estimate from a 5th order RK
    // estimate.
    double  r1  =   2.777777777777778e-03;  //  1/360
    double  r3  =  -2.994152046783626e-02;  // -128/4275
    double  r4  =  -2.919989367357789e-02;  // -2197/75240
    double  r5  =   2.000000000000000e-02;  //  1/50
    double  r6  =   3.636363636363636e-02;  //  2/55

    // Coefficients used to compute 4th order RK estimate
    // double   c1  =   1.157407407407407e-01;  //  25/216
    // double   c3  =   5.489278752436647e-01;  //  1408/2565
    // double   c4  =   5.353313840155945e-01;  //  2197/4104
    // double   c5  =  -2.000000000000000e-01;  // -1/5
    
    // Coefficients used to compute 5th order RK estimate
    double  d1  =  1.18518518519000000e-1;  //  16/135
    double  d3  =  5.18986354776000000e-1;  //  6656/12825
    double  d4  =  5.06131490342000000e-1;  //  28561/56430
    double  d5  = -1.80000000000000000e-1;  // -9/50
    double  d6  =  3.63636363636363636e-2;  //  2/55

    
    // Initialize arrays that will be returned
    double  Hold = Initheight;
    double  Vold = Initvelocity;
    double  Time = 0.0;
    double  dtmax = 0.01;
    double  dtmin = 1e-15;
    double  dt = (dtmax - dtmin)/2.0;
    double  tol = dt/pow(10,6);

    while (Time <= Endtime){
        // Compute values needed to compute truncation error estimate and
        // the 4th and 5th order RK estimate.
        //---------------------------------------------------------------------------------------------------------------
        // Runge-Kutta-Fehlberg (4th and 5th order) Integration Method
        //---------------------------------------------------------------------------------------------------------------
        double  k12 = Vold;
        double  k15 = accel(Hold, k12, Time);
        double  k22 = Vold + dt*b21*k15;
        double  k25 = accel(Hold + b21*k12, k22, Time + a2*dt);
        double  k32 = Vold + dt*(b31*k15 + b32*k25);
        double  k35 = accel(Hold + b31*k12 + b32*k22, k32, Time + a3*dt);
        double  k42 = Vold + dt*(b41*k15 + b42*k25 + b43*k35);
        double  k45 = accel(Hold + b41*k12 + b42*k22 + b43*k32, k42, Time + a4*dt);
        double  k52 = Vold + dt*(b51*k15 + b52*k25 + b53*k35 + b54*k45);
        double  k55 = accel(Hold + b51*k12 + b52*k22 + b53*k32 + b54*k42, k52, Time + a5*dt);
        double  k62 = Vold + dt*(b61*k15 + b62*k25 + b63*k35 + b64*k45 + b65*k55);
        double  k65 = accel(Hold + b61*k12 + b62*k22 + b63*k32 + b64*k42 + b65*k52, k62, Time + a6*dt);
        //---------------------------------------------------------------------------------------------------------------        
        //double Hnew = Hold + (c1*k12 + c3*k32 + c4*k42 + c5*k52)  // RK4
        double H_new = Hold + dt*(d1*k12 + d3*k32 + d4*k42 + d5*k52 + d6*k62);  // RK5
        
        //double Vnew = Vold + (c1*k15 + c3*k35 + c4*k45 + c5*k55)  // RK4
        double V_new = Vold + dt*(d1*k15 + d3*k35 + d4*k45 + d5*k55 + d6*k65);  //RK5
        //---------------------------------------------------------------------------------------------------------------
        double  r_H = (r1*k12 + r3*k32 + r4*k42 + r5*k52 + r6*k62)*dt; 
        double  r_V = (r1*k15 + r3*k35 + r4*k45 + r5*k55 + r6*k65)*dt;
        //double    r_H = (H_new - Hnew)/dt;
        //double    r_V = (V_new - Vnew)/dt;
        double  r = sqrt(pow(r_H,2) + pow(r_V,2));
        
        if (r <= tol){
            Time += dt;
            Hold = H_new;  // RK5
            Vold = V_new;  // RK5
        }  

        // Now compute next step size, and make sure that it is not too big or too small.
        dt = dt*fmin(fmax(pow((tol/(2.0*r)),0.25), 0.1), 4.0);
        
        if (dt > dtmax){
            dt = dtmax;
        }
        else if(dt < dtmin){
            printf("Error: stepsize should be smaller than: %lf.\n",dtmin);
            break;
        }
        if (H_new < 0.0){
            printf("The rocket has crashed.\n");
            break;
        }
        // Print's and save's
        std::cout << "\nFlight Data: " << "Time = " << Time << " [secs]" << "  " << "Altitude = " << H_new/pow(10,3) << " [km]" << "  " << "Velocity = " << V_new << " [m/s]" << std::endl;
        std::cout << "Erro rel.= " << r << "     " << "Step = " << dt << std::endl;
        fprintf (file0, " %f\n", Time);
        fprintf (file1, " %f\n", H_new);
        fprintf (file2, " %f\n", V_new);
        fprintf (file3, " %f\n", r);
        fprintf (file4, " %f\n", dt);
    }
}
