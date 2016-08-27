# hbv_hydromodel
Python implementation of the HBV Hydrological Model
This is a working respository and no code found here is guaranteed to be in a working state.
Check back often for updates.

###Python Dependencies
* The script requires NumPy and Matplotlib (for optional plotting of hydrograph)

###PEST Calibration Files Added
* Read more about how PEST works: [PEST] (http://www.pesthomepage.org/)
* At the above link you can download PEST. Be sure to add the relevant (pest.exe/sceua_p.exe) to you environment path variable.
* This will let you call PEST from the working directory.


#### To perform a calibration run
* Notes: *.pst is the PEST control file. It contains: the intial/min/max parameter ranges, the command line call to the model (in this case the .py file), the template and instruction files. These files tell PEST how to read/write the parameter file as well as how to interpret the model output.
* From the command line:> sceua_p hbv_pestcontrol.pst and follow the prompts (see PEST manual for details)
* PEST will then do the following: 1) run the model with the initial parameter set, 2) read the model output, 3) evaluate the objective function and perturb the parameters according to the algorithm you used until a convergence criteria is met.
