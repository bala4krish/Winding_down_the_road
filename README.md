# Winding_down_the_road
Numerical Optimization library implementing Gradient Descent and Newton's method using backtracking line search for finding the minimum.

The implementation by default minimizes the famous **Rosenbrock Banana fucntion**, but the user is free to input a function of choice as a string. 

**Code Limitations**

    1.	The code doesn’t check for the convexity of the function. 
    2.	The code can handle a polynomial function with only up to two variables.
    3.	The cost function inputted as a string is required to be in expanded form, without parentheses,

    These limitations will be a reasonable update for the future work.

<p align="center">
  <img width="500" height="380" src="https://user-images.githubusercontent.com/56740627/103165023-71e79b00-47c7-11eb-827b-d862d71622cd.png">
</p>

**Steepest Descent** 

         x(k+1) <- x(k) - a(k) f'(x(k))

**Newton's Method**

         x(k+1) <- x(k) - a(k) f''(x(k))^-1 f'(x(k))

**Backtracking Line Search**
        
         while(f(x(k) - a(k) f'(x)) > (f(x(k)) - c a(k) f'(x(k))' p(k))
                a(k) <- rho * a(k);

         where,
              a(k) = learning rate
              c    = [0,1], usually 0.5
              rho  = [0,1], usually 0.8
              p(k) = f'(x(k)) for Steepest Descent
                   = f''(x(k)) f'(x(k)) for Newton's Method

**Code Structure**

The C++ solution has a class Rosenbrock that has data and member functions to support the implementation of the two algorithms with backtracking line search. The constructor of the class takes the following initialization variables.

**Note** : The code architecture supports optimization of any convex function with two variables provided as a string during constructor initialization.

**Constructor**

    Rosenbrock(float fs, const float rho, std::pair<float, float> fst, std::string fun, std::vector<const char> fin, const float fsd=0)
        : step(fs), rho(rho), start(fst), fn(fun), xyz(fin), slimit(fsd) {}

    Explanatory variables : 
          step   -> learning rate
          rho    -> multiplier in backtracking (provided here for user convenience)
          start  -> initial starting point of minimization
          fn     -> cost function as string
          xyz    -> array of partial derivative variables
          slimit -> stopping limit for f(x) value between two iterations (default 0)
    
    
**Performance Discussion**

    Memory
        Each iteration of Newton’s method requires O(n^2) storage for the n X n Hessian matrix and O(n) gradient storage for n-dimensional gradient vector.
        Steepest descent doesn’t require storage for the Hessian matrix.
    
    Time
        Solving the inverse of Hessian take O(n^3) time. And solving the gradient takes O(n) time. 
        Steepest Descent method doesn’t have to compute inverse. So, it is faster between iterations compared to Newton’s.
    
    Performance
        The steepest descent performs a linear minimization and the Newton’s method performs quadratic approximation. So, Newton’s method is quick to converge.
    
    Saddle trap
        Newton’s method may attract saddle when used without backtracking or other adaptive techniques. Because Newton’s dynamics doesn’t treat saddle appropriately.
        Ratio of saddle and local minima increase with dimensionality. But Steepest Descent repels away saddle. So, they are more robust than Newton’s.

    
    Ill-Conditioning
        Newton’s inverse computation can sometime explode when the Hessian matrix is badly scaled or close to a singularity. This results in a failure to minimize evenwhen a function has global optimum. But Steepest Descent doesn’t get affected by this ill-conditioning. 
        Example – 
                f(x)= (x-y)^2. The hessian matrix for this function is singular.

    Smoothness
        Newton’s method lacks the smoothness that Steepest Descent guarantees. 


**Conclusion**

    Steepest Descent is many cheap steps and Newton’s is few costly steps.
