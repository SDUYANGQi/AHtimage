### This is a brief introduction for  all the R codes for simulations in my paper (Additive Hazards Model with Time-Varying Coefficients and Imaging Predictors)


Most important things!!!!!!!!!!!

Please make sure the results are stored in one file and the result analysis reads from the correct result storage file!!!!!!!!!


For Simulation 1:

AHtimage_simu1-1.R  + AHtimage_simu1-1_read.R (result summary code)
1. The setting in this code (for n1c1-) is n = 500, censoring rate approximately 15%, repeated for 1000 times.
2. The result is recorded and written out as txt file for each repeat. 
3. The result can be analyzed by simu1-1_read.R, yielding the CP table and the estimated vs. true functional curves.

AHtimage_simu1-2.R  + AHtimage_simu1-2_read.R (result summary code)
1. The setting in this code (for extra-) is n = 600, censoring rate approximately 60%, repeated for 1000 times.
2. The result is recorded and written out as txt file for each repeat. 
3. The result can be analyzed by simu1-2_read.R, yielding the CP table and the estimated vs. true functional curves.


For Simulation 2:

AHtimage_simu2.R  + AHtimage_simu2_read.R
1. Please decompress the (((slice.zip))) file and make sure all slice images are under: {   (current directory)/slice/   }.
2. There are many reading and write out instructions, please make sure that they are correct in your own computer.
3. gamma1.nii.gz and gamma2.nii.gz are the designed images for the \gamma(v,t), see Simulation 2 for more details.


For questions about these codes, please contact: {qiyang-sdu@outlook.com}.
