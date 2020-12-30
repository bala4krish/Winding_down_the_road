//
//  main.cpp
//  Numerical_Optimization
//
//  Created by Balaji Krishnamurthy on 12/20/20.
//  Copyright © 2020 Balaji Krishnamurthy. All rights reserved.

//  Rosenbrcok Banana Function
//          f(x)   =  100(x2 - x1^2)^2 + (1 - x1)^2
//          f(x)   =  100x2^2 - 200x2x1^2 + 100x1^4 + 1 - 2x1 + x1^2
//          f'(x)  =  [ -400x1x2 + 400x1^3 - 2 + 2x1,
//                       200x2 - 200x1^2               ]
//          f''(x) =  [-400x2 + 1200x1^2,     -400x1   ],
//                    [-400x1,                  200    ]]

//  Steepest Descent,
//          x(k+1) <-- x(k) - a(k) f'(x(k))
//
//  Newton's Method
//          x(k+1) <-- x(k) - a(k) f''(x(k))^-1 f'(x(k))
//
//  Backtracking Line Search
//          while( f(x(k) - a(k) f'(x)) > ( f(x(k)) - c a(k) f'(x(k))' p(k) )
//                      a(k) = rho * a(k);
//
//          where,
//                a(k) = learning rate
//                c    = [0,1], usually 0.5
//                rho  = [0,1], usually 0.8
//                p(k) = f'(x(k)) for Steepest Descent
//                     = f''(x(k)) f'(x(k)) for Newton's Method


#include "Optimizer.hpp"


int main(int argc, const char * argv[]) {
      
      std::cout << "Please input the cost function in expanded form without parentheses \n\n";
      std::cout << "Example \n\t  x^2 + y^2 \n\nOtherwise, press enter to proceed with the default. The dafault is Rosenbrock Banana function. \n\n";
      std::string func = "100y^2 - 200x^2y + 100x^4 + x^2 - 2x + 1";
//      std::string func = "3x + x^2 - 8y + 2y^2";
      std::string user_fun;
      std::getline(std::cin, user_fun);
      if(user_fun.size())
            func = user_fun;
      const float lr = 1;
      const float rho = 0.8;
      const std::pair<float,float> start_point = {1.2, 1.2};
      const std::vector<const char> xyz = {'x', 'y'};
      const float stop_limit = 0.0;
      int choice = 0;
      std::cout << "\nEnter your choice of optimization \n\t1 - Steepest Descent \n\t2 - Newton's Method \n\t3 - Both\n\n";
      std::cin >> choice;
      Rosenbrock objRs(lr, rho, start_point, func, xyz, stop_limit);
      switch (choice) {
            case 1:
                  objRs.steepestDescent();
                  break;
            case 2:
                  objRs.NewtonMethod();
                  break;
            default:
                  objRs.steepestDescent();
                  objRs.NewtonMethod();
                  break;
      }
      
      return 0;
}
