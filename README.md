# Geodetic inversion README

## Citation

Yabuki, T. and M. Matsu'ura (1992): A geodetic data inversion using a Bayesian information criterion for spatial distribution of fault slip, *Geophys. J. Int.* **109**, 363-375

## Directories

	|--exe # All compiled scripts are kept under here.
	|--fort # F90 source code
	|--fort_160129 # F90 source code
	|--plate_boundary # Fault geometry data files are kept under here.
	|--result # Setup file (.ctp), sample scripts, and plots of the result.
	   |--data # Output files of the inversion analysis are saved.
	   |--fwd # Input files of the inversion analysis should be kept here.
	|--tmp # Where the temporary files created during inversion analysis are saved.

## How to run geodetic inversion

### Necessary script

*gdinv.exe*: The main script for running a geodetic inversion analysis.  
This script is kept under */exe*, and is the only script necessary for running inversion analysis (other scripts in */exe* are not needed).

### Necessary input files

1. Observation file (*.dat*):
Input file of the coordinates and velocities of the observation points.  
The file should be kept under */result/fwd*. (Sample file: */result/fwd/vel_data_YN.dat*)

2. Fault geometry file (*.cof*):
Input file of the fault geometry.  
The *.cof* file used in the inversion analysis should be kept under */result*.
(Sample file: */result/NE_gsi3.cof*)

3. Setup file (*.ctp*):
Setup file that contains inversion parameters and file names of the input/output.  
(See sample file */result/c0_inversion.ctp* for details)

## Sample script

There is a sample cshell script under */result* that executes inversion analysis using sample observation data and setup file.

	csh inversion.csh

### Inside inversion.csh

	$set ctpname = "c0_inversion.ctp" # Set the ctp file (sample: c0_inversion.ctp).

	$../exe/gdinv.exe ${ctpname} # Geodetic data inversion analysis is executed,
	                             # using the sample setup file (c0_inversion.ctp).

	$sh plot_results.sh # Sample GMT (Generic Mapping Tool) script
	                    # that plots the results of the inversion analysis.

For details of the Generic Mapping Tools, see http://gmt.soest.hawaii.edu/.

### Results of the analysis

The result files of the inversion analysis are saved under */result/data* directory.

	res.scd # output file of the source displacement

	res.cal # output file of the calculated displacements at the input observation points

	res.fmp # output file of the B-splined source displacement

	res_cal.dat # calculated displacements at the input observation points

	res_obs.dat # observed displacements at the input observation points

	res_o_c.dat # residuals of the observed vs calculated displacements

### Results that are plotted

The plots of the results are saved under */result*.

	cnt_res.png # Contour map of the slip deficit rate.

	scd_res.png # Vector map of the source displacement.

	fwd_res.png # Vector map of the observed vs calculated displacements.

	o_c_res.png # Vector map of the residuals of the observed vs calculated displacements.

## Other scripts contained in this package

All compiled scripts are kept under */exe*.

	cgscd.exe # A script for scaling the output data,
	          # used together with the sample GMT script (plot_results.sh).
	          # This script has no effect on the result of the analysis.

	mkscd.exe # A script that creates source displacement files from a given slip distribution data.

	gdjacfw.exe # A script for calculating the Jacobian matrix of forward analysis.

	gdfwd.exe # Main script for running forward analysis.

	gdjac.exe # This script was used to calculate the Jacobian matrix for the inversion analysis,
	          # but has been deprecated since the Jacobian matrix calculation was merged to gdinv.exe.

## Source code

The source codes for the main program of inversion analysis (*gdinv.exe*) are under */fort_160129*. The codes for other scripts (listed above) are kept under */fort*.

## Compiling the scripts

Makefiles of the scripts are under */fort* and */fort_160129*, and the bash shell script for running these makefiles are under */exe*. When compiling all scripts, execute the command below.

	$sh ./exe/compile.sh
