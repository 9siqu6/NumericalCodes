# NumericalCodes
A project developed for building a very basic console application, which will solve various numerical problems such as: Linear equations, Non-linear equations, Differential equations and Matrix inversion
*Description:*

### LINEAR EQUATIONS:
### _Jacobi Iterative Method:_

**Algorithm of the code :**

Steps involved :
1.	Input Setup :
  +	Starting with a matrix mat of size n x (n+1) , where the last column represents the constant vector B . Also have an initial guess vector X, the maximum number of iterations Max_Ite, and a tolerance level tol
2.	Matrix and vector Extraction :
  +	Separating the input matrix into two parts:
  Matrix A containing the coefficients of the equations
  Vector B containing the constants.
3.	Diagonal Dominance Check :
  Checking if the matrix A is diagonally dominant . A matrix is diagonally dominant if for each row, the absolute value of the diagonal element is greater than the sum of the absolute values of the other elements in that row. If the matrix is not diagonally dominant, printing a warning message
4.	 Initialization:
  +	Initializing a previous solution vector preX with the initial guess vector X.
5.	Iteration Process:
  +	For each iteration up to the maximum allowed iterations:
  +	For each row in the matrix:
  +	Setting the initial sum to the corresponding value in vector B.
  +	Subtracting the contributions of all other variables using the previous solution vector preX.
  +	Updating the current variable in vector X by dividing the sum by the diagonal element of the current row.
6.	Convergence Check:
  +	After updating all variables in X, calculate the maximum difference between the current values in X and the previous values in preX.
  +	If the maximum difference is less than the specified tolerance tol, printing a convergence message and return the current solution vector X.
7.	  Updating Previous Values:
  +	Before the next iteration, updating preX to hold the current values of X for comparison in the next iteration.
8.	  Handling Non-Convergence:
  +	If the algorithm reaches the maximum number of iterations without achieving convergence, printing a message indicating non-convergence and return the current approximation of X.

   **Processes Inside:**
1.	bool diagonally_dominant (A) function:

This function checks if a given matrix is diagonally dominant . Diagonal dominance is an important property for convergence in iterative method . A matrix A is diagonally dominant if for each row
     A[i][i] is greater than or equals to the sum of other elements of the same row.
          
Process :
This function check matrix A for being diagonally dominant . If it is diagonally dominant , the function returns true , otherwise false.

2.	Jacobi_iterative(mat,X,Max_ite,tol) function:
  This function performs the jacobi iterative method to approximate a solution to the system of linear equations 
  AX = B

Process:
1.	Extracting Matrix A and vector B from the augmented matrix mat :
  +	A[i][j] = mat[i][j] -> Copies coefficient matrix A.
  +	B[i] = mat[i][n] -> Stores constants vector B.
2.	Diagonal Dominance check :
  +	Calls diagonally_dominant(A) .If A is not diagonally dominant,it warns that the jacobi iteration may not converge reliably.
3.	Initializing previous guess :
  +	preX = X -> sets preX as the previous guess (initially the same as X )
4.	Iteration loop:
  +	For each iteration from 1 to Max_Ite:
  For each row i in A:
  		Calculating new Approximation for X[i]:
  +	sum = B[i] -> initializing sum to the constant term B[i].
  +	For each j≠i , subtracting the product of A[i][j] and preX[j] from sum.
  sum -= A[i][j] * preX[j]
  +	Dividing sum by the diagonal element A[i][i] to update X[i]
  X[i] = sum / A[i][i].
  
  Convergence check:
  + Calculating the difference between X and preX as the maximum absolute difference between corresponding elements.
  + If diff < tol , the algorithm has converged , and it returns X with the number of iterations taken.
  +	Updating preX = X for the next iteration
5.	Non -Convergence handling :
  If the method fails to converge within Max_Ite iterations, it outputs a message and returns X as the best attempt.

Key Points:
  +	Convergence: The algorithm stops either when convergence is reached within the tolerance tol or upon the maximum number of iterations( Max_Ite)
  +	Limitations : The jacobi method may not converge if the matrix is not diagonally dominant .
 



###	_Gauss-Seidel Iterative Method:_
The Gauss-Seidel method is another iterative technique for approximating solutions to a system of linear equations 
AX = B
Like the Jacobi method, it relies on matrix diagonal dominance for convergence. However, Gauss-Seidel can converge faster as it uses updated values immediately within the same iteration.

Algorithms of the code:

Steps Involved:

1.	 Input Setup:
  +	Starting with a matrix mat of size n×(n+1) ,where the last column represents the constant vector B. Also, have an initial guess vector X, the maximum number of iterations Max_Ite and a   tolerance level tol.
2.	 Matrix and Vector Extraction:
  +	Separating the input matrix into two parts:
  +	Matrix A containing the coefficients of the equations.
  +	Vector B containing the constants.
3.	  Diagonal Dominance Check:
  +	Checking if the matrix A is diagonally dominant. This means that for each row, the absolute value of the diagonal element must be greater than the sum of the absolute values of the       other elements in that row. If the matrix is not diagonally dominant, print a warning message indicating that the Gauss-Seidel method may not converge and could produce incorrect results.
4.	  Initialization:
  +	Initializing a previous solution vector preX with the initial guess vector X.
5.	  Iteration Process:
  +	For each iteration up to the maximum allowed iterations:
  +	For each row in the matrix:
  +	Starting by setting the initial sum to the corresponding value in vector B.
  +	Subtracting the contributions of all other variables using the current values in vector X.
  +	Update the current variable in vector X by dividing the adjusted sum by the diagonal element of the current row.
6. Convergence Checking:
+	After updating all variables in X, calculate the maximum difference between the current values in X and the previous values in preX.
  +	If the maximum difference is less than the specified tolerance tol, printing a convergence message and return the current solution vector X.
7.  Updating Previous Values:
  +	Before the next iteration, updating preX to hold the current values of X for comparison in the next iteration.
8.  Handling Non-Convergence:
  +	If the algorithm reaches the maximum number of iterations without achieving convergence, printing a message indicating non-convergence and return the current approximation of X.

**Processes Inside:**

1. GaussSeidel_iterative(mat,X,Max_Ite,tol)
This function performs the Gauss-Seidel iterative method to approximate a solution to the system 
       AX  =  B
-	Parameters:
  +	mat: The augmented matrix containing coefficients and constants of the system.
  +	X: The initial guess vector for the solution.
  +	Max_Ite: The maximum number of iterations allowed for convergence.
  +	tol: The tolerance level for convergence.
-	Process:

  1.	Extracting Matrix A and Vector B from mat:
    +	A[i][j] = mat[i][j] — Copies coefficient matrix A.
    +	B[i] = mat[i][n] — Stores constant terms in vector B.
  2.	Checking Diagonal Dominance:
    +	Calls diagonally_dominant(A). If not diagonally dominant, prints a warning message.
  3.	Initializing preX = X as the previous guess.
  4.	Iteration Loop:
    +	For each iteration from 1 to Max_Ite:
    +	For each row i, calculating sum = B[i].
    +	For each j≠i , subtracting the product of A[i][j] and X[j] from sum.
    +	Updating X[i] = sum / A[i][i].
    +	Calculating diff, the maximum absolute difference between X and preX.
    +	If diff < tol, return X with convergence message.
    +	Update preX = X.
  5.	Non-Convergence Handling:
    +	If no convergence within Max_Ite, print a message and return X.

Key Points:
  +	Convergence: The algorithm stops when it reaches the tolerance tol or the maximum number of iterations (Max_Ite).
  +	Limitations: Gauss-Seidel may not converge if the matrix is not diagonally dominant.



### _Gaussian:_

**Algorithm of the Code:**

This code performs Gaussian Elimination on a given matrix `mat` to convert it into an upper triangular form. The goal of Gaussian Elimination is to simplify the matrix, allowing easier solving for unknowns in a system of linear equations.

Steps involved:
1. Partial Pivoting: Swap rows if necessary to bring the largest element (in absolute terms) in each column to the diagonal position, enhancing numerical stability.
2. Row Scaling: Normalize each row by dividing by its pivot element (the diagonal element).
3. Row Elimination: For each row, eliminate elements below the pivot by subtracting an appropriate multiple of the current row from the rows below.

**Processes Inside:**

1. Initialization: 
   - `n` is set to the matrix's size (number of rows).
   - A vector `piv` is created to store the pivot elements for each row after they are normalized.

2. Loop through Rows (Outer Loop `i`):
   - Partial Pivoting: For each row `i`, a nested loop (`k = i + 1` to `n`) compares each element in the current column with the diagonal element. If any element below the diagonal is larger, it swaps the rows.
   
3. Normalize Pivot Row:
   - The diagonal element (pivot) for row `i` is stored in `pivot`.
   - If the pivot is non-zero, it is saved in `piv`, and the entire row `i` is divided by `pivot` to make the pivot equal to 1.

4. Eliminate Entries Below Pivot (Inner Loop `k`):
   - For each row `k` below the current row `i`, a `factor` is calculated as the value in the `k`th row divided by the pivot element.
   - The `k`th row is then updated by subtracting `factor * (current row i)` from it, effectively making the element below the pivot zero.

5. Return Result:
   - The modified matrix `mat` (now in upper triangular form) is returned.

### _Gauss-Jordan:_

**Algorithm of the Code:**

This code implements the Gauss-Jordan elimination method, an extension of Gaussian Elimination that transforms the matrix into a reduced row echelon form (RREF). In this form, each pivot is 1, and all entries in each pivot's column (both above and below) are zero. This results in a fully simplified matrix, making it straightforward to solve systems of linear equations.

Steps involved:
1. Partial Pivoting: Similar to Gaussian Elimination, rows are swapped to bring the largest element (in absolute terms) in each column to the pivot position.
2. Row Normalization: Normalize each pivot row so that the pivot element becomes 1.
3. Row Elimination Below the Pivot: For each pivot, eliminate all elements below the pivot.
4. Row Elimination Above the Pivot: Starting from the last row upwards, make all elements above each pivot zero.

**Processes Inside:**

1. Initialization:
   - `n` is set to the number of rows in `mat`.
   
2. Loop through Rows (Outer Loop `i`):
   - Partial Pivoting: For each row `i`, a nested loop (`k = i + 1` to `n`) compares each element in the current column with the diagonal element. If a larger element is found below the current pivot, rows are swapped to enhance stability.
   
3. Normalize Pivot Row:
   - The diagonal element (pivot) for row `i` is stored in `pivot`.
   - If the pivot is non-zero, each element in row `i` is divided by `pivot` to make the pivot equal to 1.

4. Eliminate Entries Below Pivot (Forward Elimination):
   - For each row `k` below the current row `i`, calculate the `factor` as the value in the `k`th row divided by the pivot element in row `i`.
   - Update row `k` by subtracting `factor * (current row i)` from it, which makes the element in column `i` zero.

5. Eliminate Entries Above Pivot (Backward Elimination):
   - A backward loop (`i = n-1` down to 0) ensures that each column containing a pivot also has zeroes in all rows above the pivot.
   - For each row `k` above the current row `i`, a `factor` is calculated to zero out the entries above the pivot.

6. Return Result:
   - The fully reduced matrix `mat` in RREF form is returned, with each pivot equal to 1 and zeros both above and below each pivot. 

**Explanation Notes:**
Gauss-Jordan elimination provides a unique solution to a system of linear equations (if it exists) by reducing the matrix to its simplest form, allowing direct identification of variable values.

### _LU Factorization:_

**Algorithm of the Code:**
1. Initialize Variables:
   - Read the number of variables \( n \).
   - Set up the extended matrix A_ext to store both the coefficients and constants.
2. Extract Coefficients and Constants:
   - Separate the coefficients matrix \( A \) and the constants vector \( b \) from A_ext.
3. Initialize \( L \) and \( U \) Matrices:
   - Set up L (lower triangular matrix) and U (upper triangular matrix) with initial values. \( L \) has ones on the diagonal, while \( U \) is initialized to zero.
4. LU Decomposition:
   - Calculate the \( L \) and \( U \) matrices using:
     - *Upper Triangle Calculation*: Calculate the elements in \( U \).
     - *Lower Triangle Calculation*: Calculate the elements in \( L \).
5. Forward Substitution:
   - Solve for \( y \) in the equation \( Ly = b \).
6. Backward Substitution:
   - Solve for \( x \) in the equation \( Ux = y \).
7. Display Solution:
   - Output the solution vector \( x \).

**Processes It Does Inside:**
1. Input Collection:
   - The function prompts the user to enter the number of variables and then takes input for coefficients and constants into an extended matrix A_ext. 
2. Matrix Separation:
   - A_ext is split into two components: A for coefficients and b for constants.
3. Matrix Initialization:
   - L is initialized as a lower triangular matrix with 1s on the diagonal, while U is initialized as an upper triangular matrix filled with zeros.
4. LU Factorization Process:
   - The code iterates through A to compute the values for L and U:
     - Upper Triangular Matrix \( U \): For each row \( i \), it calculates elements of U for columns.
     - Lower Triangular Matrix \( L \): For each row \( i \), it calculates elements of L for rows \( j > i \), using values from U.
5. Forward Substitution for \( Ly = b \):
   - The code solves for y by iterating from the first row to the last, ensuring the values of y satisfy \( Ly = b \).
6. Backward Substitution for \( Ux = y \):
   - The code solves for x by iterating from the last row to the first, ensuring the values of x satisfy \( Ux = y \).
7. Solution Output:
   - The function outputs the computed values of \( x \), representing the solution to the system.

### NON LINEAR:


###	_Bisection Method:_ 

The Bisection Method is an iterative technique for finding a root within a specified interval [a, b]. It relies on the principle that if a continuous function changes sign within an interval, then there must be a root in that interval.

**Algorithms  of the code:**
Steps involved:

1.	 Input Parameters:

The algorithm requires:
  +	A vector of coefficients cof representing a polynomial.
  +	A start value start and an end value end to define the interval for root finding.
  +	A step size that determines the granularity of the intervals checked for roots.
  +	A tolerance level tol to determine the acceptable error margin for root approximation.

2.	  Finding Intervals:

  +	The function intervals is defined to find all intervals in which the polynomial changes sign, indicating the presence of a root.
  +	It initializes two variables, a (starting point) and b (end point for the current interval), and checks pairs of values from start to end using the given step.
  +	For each pair (a, b), it evaluates the polynomial at these points using the function f(cof, x). If the product of f(cof, a) and f(cof, b) is less than or equal to zero, it means there is at least one root in the interval [a, b]. This interval is then added to a list of intervals.

3.	 Bisection Algorithm:

  +	The bisection function takes an interval [a, b] and checks if it is valid for finding a root (i.e., if the function values at a and b have opposite signs).
  +	It calculates the midpoint c of the interval and evaluates the function at this midpoint.
  +	If the absolute value of the function at c is less than the specified tolerance tol, it indicates convergence, and the root c is returned.
  +	If the function changes sign between a and c, the search interval is updated to [a, c]. Otherwise, it updates to [c, b].
  +	This process repeats until the desired tolerance is achieved or the maximum interval width is reached.

4.	  Finding All Roots:

  +	The bisection_method_roots function utilizes the previously defined functions to find all roots of the polynomial within the given interval [a, b].
  +	It first calls the intervals function to obtain all intervals where roots may exist.
  +	For each interval, it calls the bisection function to find the root and adds it to the list of roots if a valid root is found.

5.	 Output:

  +	The final output of the bisection_method_roots function is a vector of roots found within the specified interval.

**Processes Inside :**
1.	Initial function Evaluations:
  +	Calculating fa = f(cof,a) and fb= f(cof,b) to evaluate the function at both ends of the interval
  +	If fa x fb >= 0 printing an error message indicating that the interval is invalid ( no sign change) and returning nan.
2.	Iteration loop:
  +	Initializing c=a to hold the midpoint of the interval
  +	While the interval width (b-a)  >= tol:
  +	Calculating the midpoint c=(b+a)/2 and evaluating fc= f(cof,c).
  +	If |fc|<tol , breaking the loop as the root is sufficiently approximated
  +	If fa x fb < 0 , updating b=c and fb=fc ( the root lies in the left half )
  +	Otherwise updating a=c and fa=fc ( the root lies in the right half)
3.	Returning value:
  +	Returing the estimated root c

Limitations:
          The bisection method requires that the function be continuous and that there is a sign change over the interval [a,b].

###	_False Position method:_

**Algorithms of the code:**

 Steps Involved:

1.	Input Parameters:

•	The function false_position_method takes the following parameters:
  +	A vector cof containing the coefficients of the polynomial.
  +	Two doubles, a and b, representing the endpoints of the interval in which to search for roots.
  +	A tolerance level tol that specifies how close to zero the function value must be for a root to be considered found.
  +	A step size used to identify potential intervals where roots may exist.

2.	 Finding Intervals:

  +	The function begins by finding potential intervals where roots may exist by calling the intervals function, passing the polynomial coefficients and the range defined by a and b, along with the step size.
  +	This returns a vector of pairs representing the intervals.
3.	  Iterating Over Intervals:
  -	For each interval (a, b) in the vector of intervals:
    +	The function evaluates the polynomial at the endpoints a and b, storing the results in f_a and f_b, respectively.
    +	It checks if the function values at the endpoints have opposite signs (f_a * f_b < 0). If they do not, it skips to the next interval since no root exists in this interval.

4.	False Position Method Loop:
  
  +	If a valid interval is found, the method enters a while loop that continues until the width of the interval (b - a) is greater than or equal to the specified tolerance tol.
  +	Inside the loop, it calculates the root using the False Position formula.
  +	It then evaluates the function at this new root (f_root).
  +	If the absolute value of f_root is less than the tolerance (tol), the root is considered found, and it is added to the roots vector, breaking the loop.
  +	If the function value at the root has the same sign as f_a, it updates the lower bound a to root and recalculates f_a. Otherwise, it updates the upper bound b to root and recalculates f_b.

5.	  Returning Roots:

  +	After processing all intervals, the function returns the vector roots, which contains all the roots found using the False Position Method.

**Processes Inside:**
1.	Finding intervals :
  + Calling the intervals function to identify intervals where the function may have roots , passing in the coefficients and the interval limits (a,b) along with the step.
2.	Iterating over intervals:
  +	For each identified interval (a,b):
  +	Calculating fa=f(cof,a) and fb=f(cof,b);
  +	If fa x fb >=0 skiping to the next interval since no root exists in this interval
  3.	Iteration loop:
  +	While the width of the interval (b-a)>=tol
  +	Calculating the root using the false position formula:
Root = b - ( (fb.(b-a))/fb-fa)
  +	Evaluating froot – f(cof, root );
  +	If |froot| < tol considering the root found and breaking the loop
  +	If fa x froot <0
  +	Updating b = root and fb = froot
  +	Otherwise 
  +	Updating a =root and fa = froot
4.	Returning value:
  +	Returning the list of found roots .


Limitations:
The false position method relies on a sign change and can be sensitive to the choice of initial interval, which may affect convergence behavior


### _Newton-Raphson Polynomial Root Finding:_

**Algorithm of the Code:**

This code applies the Newton-Raphson method to find all unique roots of a polynomial function defined by the coefficients in a vector. The function `f(x)` represents a polynomial of degree 4, and `f'(x)` is its derivative, both used to iteratively approximate roots. Starting with a range of initial guesses, the code converges to the unique roots of the polynomial and stores them.

Steps involved:
1. Define Polynomial and Derivative: Functions `f(x)` and `f'(x)` are created to evaluate the polynomial and its derivative for any `x` value.
2. Unique Root Check: Helper function verifies if a computed root is distinct and hasn't been found yet.
3. Newton-Raphson Iteration: For each initial guess, the algorithm refines a root approximation using Newton-Raphson iterations until convergence within a specified tolerance.
4. Collect and Sort Results: Stores unique roots in ascending order and returns them.

**Processes Inside:**

1. Function Definitions:
   - `f(vf, x)` calculates the polynomial value for a given `x` using coefficients `a`, `b`, `c`, `d`, and `e` from vector `vf`.
   - `fp(vf, x)` computes the polynomial's derivative at `x`.

2. Unique Root Verification (`isUniqueRoot`):
   - Takes a list of roots (`vc`) and checks if the calculated root (`root`) is unique by comparing it with each existing root in `vc`.
   - If any stored root differs from the new one by less than `0.0001`, the root is considered non-unique.

3. Newton-Raphson Iteration (`newtonraphson`):
   - Initializes an empty vector `vc` to store unique roots and sets a tolerance `tolr`.
   - For each integer `i` between `-maxAt` and `maxAt`, an initial guess `x0` is calculated as `i - maxAt / 2`.
   - Inner Iteration (for each guess):
     - For each `x0`, computes up to `maxIt` iterations to refine the root approximation.
     - Steps per iteration:
       - Calculates `f(x0)` and `f'(x0)`.
       - If `f'(x0)` is too close to zero, exits to avoid division by zero.
       - Computes the new approximation `x1 = x0 - f(x)/f'(x)`.
       - If `|x1 - x0|` is within `tolr`, it checks if `x1` is unique using `isUniqueRoot`. If unique, adds `x1` to `vc`.
     - Updates `x0 = x1` for further refining if not converged.

4. Result Collection:
   - After collecting all unique roots, sorts `vc` in ascending order and returns it.

**Explanation Notes:**
The Newton-Raphson method iterates over a range of initial guesses to approximate each root of the polynomial, ensuring high accuracy by checking convergence and uniqueness.

### _Secant Method Polynomial Root Finding:_

**Algorithm of the Code:**

This code uses the Secant Method to approximate roots of a polynomial. It iteratively finds the points where the polynomial value approaches zero by using two initial guesses and updating them until convergence to a root. It checks each result for uniqueness before storing it.

Steps involved:
1. Initialize: For each integer `i` in a specified range, two starting points (`x1` and `x2`) are chosen.
2. Secant Iteration: Each pair is refined through Secant iterations to approximate a root.
3. Root Verification: The approximated root is checked for uniqueness before storing.
4. Result Collection: All unique roots are collected, sorted, and returned.

**Processes Inside:**

1. Initialize Variables:
   - Takes the polynomial coefficients in vector `vf`.
   - Defines `vc` to store unique roots and sets tolerance `tol` for convergence.

2. Outer Loop (Starting Guesses):
   - Iterates through initial guesses from `-maxAt` to `maxAt`.
   - For each integer `i`, sets `x1 = i` and `x2 = i + 1` as the starting points.
   - Calculates initial polynomial values `fx1 = f(vf, x1)` and `fx2 = f(vf, x2)`.

3. Secant Iteration (Inner Loop):
   - Up to `maxIt` iterations are performed to find a root approximation.
   - Per Iteration:
     - Checks if `|fx2 - fx1|` is too small (near zero slope); breaks to avoid division by zero.
     - Calculates the next approximation `x3` using the Secant formula: \[ x3 = \frac{x1 \cdot fx2 - x2 \cdot fx1}{fx2 - fx1} \]
     - If `|f(vf, x3)| < tol`, the method considers `x3` a root. It checks if `x3` is unique using `isUniqueRoot` and, if so, adds it to `vc`.
     - Updates `x1 = x2`, `fx1 = fx2`, `x2 = x3`, and `fx2 = f(vf, x3)` for the next iteration if not converged.

4. Result Collection:
   - After identifying all unique roots, the vector `vc` is sorted in ascending order and returned.

**Explanation Notes:**
The Secant method requires two initial guesses but doesn’t need the derivative, making it a good alternative to Newton-Raphson. By iteratively updating the guesses, it approximates the roots, ensuring convergence and uniqueness through checks.

### DIFFERENTIAL EQUATIONS:

### _Runge-Kutta Method for Differential Equation Approximation:_

**Algorithm of the Code:**

This code uses the 4th-order Runge-Kutta method to approximate solutions to a differential equation over a specified range. Given initial conditions and step size `h`, it iteratively calculates the values of `x` and `y` by applying the Runge-Kutta formulas.

Steps involved:
1. Initialization: The initial values of `x` and `y` are set.
2. Runge-Kutta Step Calculation: Four intermediate slopes (`k1`, `k2`, `k3`, and `k4`) are computed using the differential equation.
3. Update Solution: Using the weighted average of the slopes, the method updates `x` and `y` to approximate the solution at each step.
4. Display Results: The results are stored in a vector and printed as `(x, y)` pairs.

**Processes Inside:**

1. Function `rungekutta`:
   - Inputs: `h` (step size), `x` and `y` (current values of the variables), and parameters `a`, `b`, and `t` for the function `f`.
   - Intermediate Slopes:
     - Calculates `k1`, `k2`, `k3`, and `k4` based on `f(x, y, a, b, t)`.
   - Update `y` and `x`:
     - Calculates `y`
     - Updates `x` by adding `h`.
   - Returns: A pair `{x, y}` representing the next approximated values of `x` and `y`.

2. Function `show_rk`:
   - Purpose: Displays Runge-Kutta solution over a range, from 0 to approximately 4*pi.
   - Steps:
     - Initializes `vxy`, a vector of `(x, y)` pairs, and starts with `xy` at `{0.0, 0.0}`.
     - Loop: Continuously calls `rungekutta` to get the next `(x, y)` values and updates `xy`.
     - Output: Prints `x` and `y` values after each step and stores the result in `vxy`.
     - Termination: Ends when `x` exceeds `4π`.

This approach is ideal for simulating solutions to ordinary differential equations where accuracy is crucial. The 4th-order Runge-Kutta method provides high accuracy by averaging multiple slopes, which smoothens the approximation.

### MATRIX INVERSION:

### Matrix Inversion Using Augmented Matrix and Gaussian Elimination

**Algorithm of the Code:**

This code finds the inverse of a given square matrix `mat` using an augmented matrix approach with Gaussian elimination. The algorithm involves creating an augmented matrix by appending the identity matrix to `mat`, then applying row operations to transform the left half into the identity matrix. The result in the right half will then be the inverse of `mat`.

Steps involved:
1. Initialize Identity Matrix: Append the identity matrix to the given matrix, creating an augmented matrix.
2. Forward Elimination: Perform Gaussian elimination to create zeros below the pivot element in each column.
3. Backward Elimination: Create zeros above each pivot, transforming the left half into the identity matrix.
4. Extract Inverse: Copy the right half of the augmented matrix, which now holds the inverse, into a separate matrix.

**Processes Inside:**

1. Identity Matrix Creation:
   - `idt` is initialized as an identity matrix of the same dimension `n` as `mat`.
   - Each row `i` of `mat` is augmented with the corresponding row of `idt`, creating a combined matrix of size \( n \times 2n \).

2. Forward Elimination (Upper Triangular Form):
   - For each row `i`, the code performs a partial pivot by finding the row with the largest absolute value in the current column to avoid numerical instability and swaps it with the current row.
   - The pivot row is normalized by dividing each element by the pivot value, ensuring that the diagonal element becomes 1.
   - For each row below `i`, row operations are performed to create zeros in the current column by subtracting a multiple of the pivot row from each row below it.

3. Backward Elimination (Identity Formation):
   - Starting from the last row, the algorithm creates zeros above each pivot to fully transform the left half of the augmented matrix into an identity matrix.
   - Each operation on row `i` clears the entries above the pivot position in each column.

4. Extracting the Inverse Matrix:
   - Once the left half of `mat` is the identity matrix, the right half represents the inverse of the original matrix.
   - A new matrix `inv` is created to store the values of the inverse, which is then returned.

This function can be used for solving systems of linear equations where matrix inversion is required, as well as for applications in linear algebra where inverse matrices are needed.


**show_matrix():**

+ This function displays a matrix in a readable format.
+ It takes a 2D vector mat as input, which represents the matrix.
+ It iterates through each element, printing each row and column.
+ The output shows "RESULTING MATRIX" followed by the formatted matrix.
+ Useful for visualizing matrices, particularly after transformations.

**returnToMainMenu():**

+ This function pauses the program and waits for user input to return to the main menu.
+ Displays a prompt to press "Enter" to return to the menu.
+ Uses cin.ignore() and cin.get() to wait for input.
+ Clears the console using system("CLS").
+ Provides a smooth transition back to the main menu.

**showMainMenu():**

+ Displays the main menu for selecting a numerical operation.
+ Lists options for different types of numerical equations and matrix inversion.
+ Outputs instructions for the user to enter their choice.
+ Simple, user-friendly interface for navigating the program.
+ Helps guide the user in making selections.

**prompt():**

+ Acts as the main control for the program's operations based on user input.
+ Shows the main menu and retrieves user choices.
+ Based on the choice, it calls functions for solving linear/non-linear equations, differential equations, and matrix inversion.
+ Handles all input and outputs results in a clear format.
+ Uses returnToMainMenu to navigate back to the menu after each operation.
