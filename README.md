# Linear-Regression-Using-Matrices
The main file contains the solution to solve an overdetermined set of equations so as to minimize the sum of square of residual errors. 

### Consider the following overdetermined system of linear equations:
    1. 0.25.x1 + 0.03.x2 = 0.97
    2. -0.3701 +0.17.x2 =0.46
    3. 1.17.x1 + x2 = 2.20
    4. -1.09.x1 -0.17.x2 = -1.19
    5. x1 + 1.19.x2 = 1.73 

The residual is given by,
    
    r= c1. x1 + c2. x2 - c3
 
To find the values of our parameters x1, x2 that minimizes the square of sum of residuals, we perform the following matrix operation:

    X= (A<sup>T</sup>A)<sup>-1</sup>A<sup>T</sup>b , given that (A<sup>T</sup>A) is invertible


