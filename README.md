# R1calc_CPlusPlus
Non linear least squares regression for exponential recovery curve using Gauss Newton Algorithm 

Non linear least squares regression for exponential recovery curve using Gauss Newton Algorithm

Program r1calc This code uses the Gauss Newton non linear regression algorithm to provide the least squares 
fit to data from the NMR inversion recovery Fourier Transform (IRFT) experiment for a single isolated spin

There are 3 optimisation parameters for which initial esimates must be provided.

These are (effectively) peak area at zero time, peak area at long time and relaxation rate R1

Given these estimates the program iterates until the final and optimum least squares solution is found.

Note that if the initial estimates are very poor then the method may not converge

Note that this code is easily modified to fit other functional forms.
Change the function called func to specify the mathematical form of the function to use for fitting. 
Change the function called deriv to specify the gradient of func with respect to the optimisation parameters. 
The optimisation parameters are held in the floating point vector called a

Compile code with command line:
g++ r1calc.cpp -lm -o r1calc

Execute code with command line 
./r1calc 

Enter number of x and y data points when prompted
Enter x (time) and y (peak area) data points when prompted 
Enter initial estimates when prompted 
Code should then iterate reducing least squares residuals 
until convergence criterion < 0.001 
The x y(experimental) and y(calculated) values will then be printed
after the optimisation parameter final values. 

Example input data - 6 data points

0 	-10

1 	-5

2 	0

3 	3

5 	4

10 	-10

Example for initial parameter estimates

x or peak area at zero time = -10

y or peak area at long time = 10

relaxation rate r1 = 0.5

