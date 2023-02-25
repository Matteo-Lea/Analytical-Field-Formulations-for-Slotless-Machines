# Analytical-Field-Formulations-for-Slotless-Machines
So you have made it this far in the repository, huh?
I appreciate you sticking around :)

Remember that this is public material you can use; if you see anything 
that ought to be improved or anything wrong, you should not hesitate to 
speak up so the model will get better and people needing it will be 
happier :)



Short guide through the different scripts in this folder:

Inrunner.m and Outrunner.m

%--------------------------------------------------------------------------

Either of these two motors are loaded every time you run one of the other 
scripts (just look for the line where they are called in each script and 
you will see  which one is called). 
They do nothing more but loading the motor geometry and some useful motor
parameters needed to solve the field problem; it's like the Parameters
definition under Global Definitions in Comsol if you are familiar with it.
If you feel confident enough for changing the parameters, do so, but make 
sure to understand what they mean first ;)

%--------------------------------------------------------------------------

PM_field_map.m

%--------------------------------------------------------------------------

Probably the most interesting and "colorful" script.
Here the field solution from the permanent magnets is automatically 
computed everywhere on the 2D motor geometry. The output are a lot of
matrices hloding the solution in each and every domain, but the scope of 
the code is to output some nice plots. You can control to plot the field 
map, the torque waveform under sinusoidal supply and the back-efm graph, 
through these three variables:
#####
mapping = "no";
torque = "no";
back_emf = "no";
#####
By default you are not going to see anything.

Easy interaction with the code:
-) Should you see a poor resolution on the field map, you can change the 
radial domains discretization within the section "RADIAL DISCRETIZATION" 
and the angular discretization in "ANGULAR DISCRETIZATION" 
(the computational time will increase)
-) The script is defaulted to output the solution over one pole. If you 
want to see more, you can modify the variable "sec" (yes, sec reminds of 
seconds and it's the built-in function for secant in Matlab, but I was too 
lazy to give it another name)

%--------------------------------------------------------------------------

Armature_field_map.m

%--------------------------------------------------------------------------

Here the field solution from the three-phase statot current is computed 
everywhere. It automatically plots the field map as it's the only feature 
in this script.

For easy interaction with the code, you may want to refer to the section 
above. If you want to change the current value, you can do that within the 
input scripts "Inrunner.m" and "Outrunner.m"

%--------------------------------------------------------------------------

Inductance.m

%--------------------------------------------------------------------------

I mean, the name says it all. This script gives you the motor inductance.
You can find the different terms in the workspace after running it:
-) self-inductance "L_self"
-) mutual inductance "L_mut"
-) synchronous inductance "L_sync"
-) end-winding inductance Ls_ew (which uses an equation I got from my 
electrical machines design notes)
-) terminal inductance assuming a star connection "L_star"

%--------------------------------------------------------------------------

Bemf_torque_constants.m

%--------------------------------------------------------------------------

Not much to add here either. This script gives you the motor back-emf and 
torque constants.
You can find the different terms in the workspace after running it:
-) back-emf/voltage constant "k_v" in [V/rad/s]
-) torque constant "k_t" in [Nm/A_peak]

%--------------------------------------------------------------------------

Speed_test.m

%--------------------------------------------------------------------------

With this script you can run a speed test of all the scripts in this folder
just change the name of the script in the for loop to one of those listed
above. Remember to remove the mapping feature within Armature_field_map.m 
and PM_field_map.m if you enabled them. You can set how many times to run 
the script through the variable "runs"

%--------------------------------------------------------------------------

If you find the topic intriguing and you believe you can contribute to 
improve this work, just do that!

NOTE: 
the major flaw of all these codes is that they fill the workspace with a 
lot of useless stuff. A much cleaner and more efficient way to structure 
everything would be to have all the scripts as callable functions so the 
output is just what you  need. Having just scripts makes it so much easier 
to debug though :)
