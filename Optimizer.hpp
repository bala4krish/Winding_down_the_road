//
//  Optimizer.hpp
//  Numerical_Optimization
//
//  Created by Balaji Krishnamurthy on 12/20/20.
//  Copyright © 2020 Balaji Krishnamurthy. All rights reserved.
//

#ifndef Optimizer_hpp
#define Optimizer_hpp
#include <iostream>
#include <limits>
#include <vector>
#include <unordered_map>

#define DOUBLE_MAX std::numeric_limits<double>::max()
#define VECTOROFPAIR std::vector<std::pair<std::string,float>>

class Rosenbrock
{
private:
      float step;
      const float rho;
      std::pair<float, float> start;
      std::string fn;
      const std::vector<const char> xyz;
      const float slimit;
      std::vector<std::pair<float,float>> SD_ans;
      std::vector<std::pair<float,float>> NM_ans;
      VECTOROFPAIR fvec;
      std::vector<VECTOROFPAIR> gradvec;
      std::vector<std::vector<VECTOROFPAIR>> hessmat;
      
public:
      Rosenbrock(float fs, float rho, std::pair<float, float> fst, std::string fun, std::vector<const char> fin, float fsd=0)
      : step(fs), rho(rho), start(fst), fn(fun), xyz(fin), slimit(fsd) {}
      std::unordered_map<std::string, int> Derivate(const char& ) const;
      double computePoly(const std::vector<std::pair<std::string,float>>& , const std::pair<float, float>& ) const;
      std::vector<float> Gradient(const std::pair<float, float>& );
      void computeFx();
      void SD_BackTrack(const std::pair<float, float>& );
      std::pair<float,float> steepestDescent();
      std::vector<std::vector<double>> Hessian(const std::pair<float, float>& );
      void Inverse(std::vector<std::vector<double>>& ) const;
      void NM_BackTrack(const std::pair<float, float>& );
      std::pair<float,float> NewtonMethod();
      ~Rosenbrock()
      {
            SD_ans.clear();
            NM_ans.clear();
            fvec.clear();
            gradvec.clear();
            hessmat.clear();
      }
};

#endif /* Optimizer_hpp */
