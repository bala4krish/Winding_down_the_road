# Winding_down_the_road
Numerical Optimization library implementing Gradient Descent and Newton's method using backtracking line search for finding the minimum.

The implementation by default minimizes the famous **Rosenbrock Banana function**, but the user is free to input a function of choice as a string. 

**Code Limitations**

    1.	The code doesn’t check for the convexity of the function. 
    2.	The code can handle a polynomial function with only up to two variables.
    3.	The cost function inputted as a string is required to be in expanded form, without parentheses,

    These limitations will be a reasonable update for the future work.
    
<p float = "left">
    <img width="310" height="210" src="https://user-images.githubusercontent.com/56740627/103203619-2954de00-48aa-11eb-919d-0730e59bde0e.png">
    <img width="310" height="210" src="https://user-images.githubusercontent.com/56740627/103203522-f14d9b00-48a9-11eb-941f-bf543f5f9d44.png">
    <img width="310" height="210" src="https://user-images.githubusercontent.com/56740627/103165023-71e79b00-47c7-11eb-827b-d862d71622cd.png">
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


**Some Results**

    Steepest Descent

            iteration : 0,	        f(x) : 0.168912,	learning rate : 0.000633826

            iteration : 1,	        f(x) : 0.0128133,	learning rate : 0.000792282

            iteration : 2,	        f(x) : 0.0125482,	learning rate : 0.000792282

            iteration : 3,	        f(x) : 0.0125255,	learning rate : 0.0037779

            iteration : 4,	        f(x) : 0.0125076,	learning rate : 0.000990353

            iteration : 5,	        f(x) : 0.0124615,	learning rate : 0.00922338

            iteration : 6,	        f(x) : 0.0124233,	learning rate : 0.000792282

            iteration : 7,	        f(x) : 0.0123059,	learning rate : 0.022518

            iteration : 8,	        f(x) : 0.0122295,	learning rate : 0.000792282
            .
            .
            .
            .
            .
            .
            .
            iteration : 979,	f(x) : 1.56258e-09,	learning rate : 0.000990353

            iteration : 980,	f(x) : 1.55835e-09,	learning rate : 0.00302232

            iteration : 981,	f(x) : 1.55087e-09,	learning rate : 0.00922338

            iteration : 982,	f(x) : 1.55087e-09, 	learning rate : 2.14526e-17


            Optimal points 	        x1 : 1.00004,	 x2 : 1.00008


    Newton’s Method

            iteration : 0,	        f(x) : 0.0383841,	learning rate : 1

            iteration : 1,	        f(x) : 0.0175462,	learning rate : 0.4096

            iteration : 2,	        f(x) : 0.00490727,	learning rate : 1

            iteration : 3,	        f(x) : 0.00083187,	learning rate : 1

            iteration : 4,	        f(x) : 4.2605e-05,	learning rate : 1

            iteration : 5,	        f(x) : 1.93131e-07,	learning rate : 1

            iteration : 6,	        f(x) : 5.05906e-12,	learning rate : 1

            iteration : 7,	        f(x) : 0,	        learning rate : 1

            iteration : 8,	        f(x) : 0,	        learning rate : 1


            Optimal points 	        x1 : 1,	        x2 : 1
