# Analytical-Field-Formulations-for-Slotless-Machines

Diving into some more applications of Sub-Domain modeling here we are using
the field solution in the stator yoke to estimate the iron losses.
All the models use predefined lamination loss coefficients based on some 
data I had (unfortunately I cannot share the data), but if you read my 
paper [1] you will get a good guide on how to find the different coefficients
for any laminations loss dataset.

[1] M. Leandro, N. Elloumi, A. Tessarolo and J. K. Nøland, "Analytical Iron
 Loss Evaluation in the Stator Yoke of Slotless Surface-Mounted PM 
Machines," in IEEE Transactions on Industry Applications, vol. 58, no. 4, 
pp. 4602-4613, July-Aug. 2022, doi: 10.1109/TIA.2022.3171528.


A short guide through the different scripts in this folder:

Inrunner.m and Outrunner.m

%--------------------------------------------------------------------------

Either of these two motors are loaded every time you run one of the other 
scripts (just look for the line where they are called in each script and 
you will see  which one is called). 
They do nothing more but load the motor geometry and some practical motor
parameters needed to solve the field problem; it's like the Parameters
definition under Global Definitions in Comsol if you are familiar with it.
If you feel confident enough to change the parameters, do so, but make 
sure to understand what they mean first ;)

%--------------------------------------------------------------------------

Iron_loss_CCM.m

%--------------------------------------------------------------------------

This code uses the iron loss model based on constant coefficients 
describing the lamination's iron loss variation with flux density and 
frequency, according to:
Ke*B^2*f^2 + Kh*B^alpha*f + Ka*B^1.5*f^1.5

'COEFF', 'B_exp', and 'f_exp' hold the loss coefficients, flux density 
exponents and frequency exponents.

'CCM_loss_freq' and 'CCM_loss_time' hold the computed iron losses

%--------------------------------------------------------------------------

Iron_loss_CAL2.m

%--------------------------------------------------------------------------

This code uses the iron loss model based on the two-components equation 
with variable coefficients and constant flux density exponent for the 
hysteresis component:
   Ke(f,B)*B^2*f^2 + Kh(f,B)*B^2*f 

The laminations loss data was divided into two frequency ranges wherein the
dependency on the frequency for the different coefficients is neglected. 
Within the two frequency ranges the coefficients for the polynomial models 
describing the iron loss coefficients are:
-)'Kh_pol1', 'Ke_pol1' --> range 1
-)'Kh_pol2', 'Ke_pol2' --> range 2

'CAL2_loss' holds the iron loss value

%--------------------------------------------------------------------------

Iron_loss_VARCOext.m

%--------------------------------------------------------------------------

This code uses the iron loss model based on the three components equation 
with variable coefficients:
   Ke(B)*B^2*f^2 + Kh(f)*B^alpha(B,f)*f + Ka(B)*B^1.5*f^1.5

The laminations loss data was divided into two frequency ranges wherein the
dependency on the frequency for the different coefficients is neglected. 
Within the two frequency ranges the coefficients for the polynomial models 
describing the iron loss coefficients are:
-)'Ka_pol1', 'Ke_pol1', 'alpha_01f', 'alpha_11f', 'alpha_21f', 'alpha_31f',
  'K_h1' --> range 1
-)'Ka_pol2', 'Ke_pol2', 'alpha_02f', 'alpha_12f', 'alpha_22f', 'alpha_32f',
  'K_h2' --> range 2

'VARCO_loss' holds the iron loss value

%--------------------------------------------------------------------------

Iron_loss_VARCOrot.m

%--------------------------------------------------------------------------

This code uses the iron loss model based on the three components equation 
with variable coefficients with the rotational correction 'rot' which is 
assumed to be dependent on the flux density only:
   (Ke(B)*B^2*f^2 + Kh(f)*B^alpha(B,f)*f + Ka(B)*B^1.5*f^1.5)*rot(B)

The laminations loss data was divided into two frequency ranges wherein the
dependency on the frequency for the different coefficients is neglected. 
Within the two frequency ranges the coefficients for the polynomial models 
describing the iron loss coefficients are:
-)'Ka_pol1', 'Ke_pol1', 'alpha_01f', 'alpha_11f', 'alpha_21f', 'alpha_31f',
  'K_h1', 'rot' --> range 1
-)'Ka_pol2', 'Ke_pol2', 'alpha_02f', 'alpha_12f', 'alpha_22f', 'alpha_32f',
  'K_h2', 'rot' --> range 2

'VARCO_rot' holds the iron loss value

%--------------------------------------------------------------------------

Iron_Loss_speed_sweep.m

%--------------------------------------------------------------------------

This code uses all the different scripts described above to evaluate the 
losses over a speed range defined in the variable 'speed'.
In this way, you will get a comparison of all the different loss models.
Suggestion: In the input file 'Inrunner.m' if you change the stator 
back-iron thickness by changing R_se you might be able to see a bigger 
difference between the different models (since for the predefined motor 
geometry all the iron loss models give a good fit to the experimental data
at the flux density value in the iron yoke)

%--------------------------------------------------------------------------

Speed_test.m

%--------------------------------------------------------------------------

With this script you can run a speed test of all the scripts in this folder
just change the name of the script in the for loop to one of those listed
above. You can set how many times to run the script through the variable 
"runs"

%--------------------------------------------------------------------------


