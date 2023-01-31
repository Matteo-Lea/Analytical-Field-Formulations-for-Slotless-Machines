# Analytical-Field-Formulations-for-Slotless-Machines
Some development and usage of analytical field models for slotless PM machines


%--------------------------------------------------------------------------------------------------
The directory named "Demagnetization Analysis and More" is a collection of codes developed during my 
PhD with the purpose of extending the standard analysis offered by analytical field formulations.
The ambitious initial plan was to account for PM demagnetization in slotless machines thrugh a model
based on the fundamental assumption considering linear material properties.
As no relationship was found between the effect of demagnetization evaluated with non-linear FEA and 
quantities found through the analytical model, an "alternative" optimization procedure is proposed.
The procedure takes a demgnetization related quantity to be considered as a penalty in the optimization 
algorithm.
Are you curious to see what this research led to? Check out this paper:

"M. Leandro and J. K. Nøland, "A Penalty-Based PSO Algorithm for Demagnetization Risk-Free Design
of Slotless Halbach PM Machines," 2022 International Conference on Electrical Machines (ICEM), 2022,
pp. 276-281, doi: 10.1109/ICEM51905.2022.9910651."
%---------------------------------------------------------------------------------------------------
In the directory named "Field Maps and More" you can find the script to plot the permanent magnets 
field map, the script to plot the armature field map and the script to calculate the motor phase
inductance throught he armature field solution. In addition the script plotting the permanet magnets
field map also gives torque waveform and voltage waveform of the loaded motor. These codes are largely
based on the work presented in:

"M. Leandro and J. K. Nøland, "An Approach for Optimal Pre-Conditioning of the Analytical Field 
Solution of Slotless PM Machines," in IEEE Access, vol. 9, pp. 36748-36765, 2021, 
doi: 10.1109/ACCESS.2021.3062769."
%---------------------------------------------------------------------------------------------------



If the topic is intriguing and you believe you can contribute to improve this work, just do that!
