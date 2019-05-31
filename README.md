# EGSim

The Electric power Grid Simulator (EGSim) software toolkit implements algorithms aimed currently at solving static load flow problems for electric power grids. It parses power grid models described in IEEE Common Data Format, and generates solutions for bus voltages and voltage angles, and real and reactive power values through the transmission lines. It is written in C++ and released under a BSD license.

Directory structure:

	- src:          source files
	- dep:          dependencies (netlib libraries)
	- config:       configuration scripts
	- example:      IEEE example models
	- manual.pdf:   theory, examples, installation notes  

To compile, create a build directory, preferably separate from the source directory. Then under the build directory:

	- make a copy of config/config-template.sh 
			cp config/config-template.sh config/config-custom.sh
			then customize for current compilers and paths
	- ../config/config-custom.sh
		- here config-custom.sh is the modified version of the config script
	- make -j 4; make install

Then switch to the example:

	- switch to the example directory
	- from command-line: ./run.sh
		+ this script requires the path to EGSim binary be set through the EGSIMBIN environment variable
		+ this script will run the 118-bus example and 300-bus example, respectively
		+ each computation outputs the following files:
			- admmat.dat  (admittance matrix)
			- admmatR.dat (admittance matrix - real part)
			- admmatI.dat (admittance matrix - imaginary part)
			- linepow.dat (line powers, with flag -p)
			- statsol.dat (power flow solution: voltages magnitude and angles)


