# Analytical-Field-Formulations-for-Slotless-Machines

This folder is the result of a pilot research where I tried to quantify
the rsk of demagnetization in slotless PM machines under no-load 
operation. For high-grade PM materials, simply increasing the rotor 
temperature with an Halbach configuration, might lead to a local
irreversible demagnetization; this leads to poor material utilization, and
reduced performance. The goal was to related the "demagnetized volume" 
(which can be done efficiently through the analytical model) to the 
performance reduction. There might be a way to find a relationship between 
the two, but I didn't have time to investigate it any further. What I did
instead, was to develop an optimization algorithm to guide the design
far from designs showing a high demagnetization risk.
This is the publication presenting the work:

[1] M. Leandro and J. K. Nøland, "A Penalty-Based PSO Algorithm for 
Demagnetization Risk-Free Design of Slotless Halbach PM Machines," 2022 
International Conference on Electrical Machines (ICEM), Valencia, Spain, 
2022, pp. 276-281, doi: 10.1109/ICEM51905.2022.9910651.


Short guide through the different folders and scripts in this folder:

BH_curve.m

%--------------------------------------------------------------------------

This is a code I wrote to get the main magnets performance metrics given
the PM BH curve. A lot of them are loaded in the folder PM_data (where the
code takes the data from). The data is written to the PM_data.csv file. 
Here you will find the PM grade name and the temperature of the input data, 
the remanence, the recoil permeability, the knee-poin flux dnesity and the 
demagnetization flux density from the identified airgap flux density


%--------------------------------------------------------------------------

IN THE FOLDER NAMED: Optimization

Demag_optim.m

%--------------------------------------------------------------------------

This is a code you ca use to visualize the demagnetization map over one 
pole of whatever motor you will give as an input. Since I used it to check
the demagnetization map of the designs coming out of the optimization, you 
will find the following three lines:

polepairs = 'p=3'; % available [3 5 7 9 11 13 15 17 19]
mode = 'NO_DEMAG'; % 'NO_DEMAG --> optimization without demagnetization
                   % 'W_DEMAG --> optimization with demagnetization
load(['RESULTS\' polepairs '\' mode '\' polepairs ])

By changing the entry of 'polpairs' and 'mode' you can access different 
designs in the folder 'RESULTS'. But you can also uncomment the lines under 
"test-geometry" and change the geometry to whatever you want to check. You 
can also chnage the PM material data under "PM material properties"

%--------------------------------------------------------------------------

IN THE FOLDER NAMED: NO_DEMAG

MAIN_optim.m

%--------------------------------------------------------------------------

This script initiates the search ranges of the different optimization
variables. The optimization algorithm is contained in the function 
"MAIN_PSO.m". As I said, I haven't spent much time tweaking the algorithm
to make something fancy, feel free to go through the implementation. The 
main flaw of this optimization is the absence of several constraints which
would probably make the algorithm discard several of the generated designs.
But you can also redefine any constraint you'd like as an objective, and 
feed it into the algorithm (probably easier said than done :)).
By default, in "MAIN_PSO" the objective is the sume of the motor active 
weight and the sum of the constraints violation given by the torque and
back-emf defined in "MAIN_optim.m".

%--------------------------------------------------------------------------

IN THE FOLDER NAMED: DEMAG_INCLUDED

MAIN_optim.m

%--------------------------------------------------------------------------

It's exactly the same as the script described previously. The only 
difference is the constraint on the admissible demagnetized volume 
percentage defined along with the same torque and back-emf defined as 
before.

%--------------------------------------------------------------------------

IN THE FOLDER NAMED: DEMAG_INCLUDED

MAIN_optim.m

%--------------------------------------------------------------------------

It's exactly the same as the script described previously. The only 
difference here is that the function "Field_Solution_Demag.m" uses the
GPU computational power to increase the computatioal performance of the
algorithm. If you don't know whether you can use your GPU or not, you can 
run:

canUseGPU()

You can run the script if you get 1 as a logical value.

%--------------------------------------------------------------------------







