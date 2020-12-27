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
//          x(k+1) <-- x(k) - a(k) f''(x(k)) f'(x(k))
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

#include <iostream>
#include <limits>
#include <vector>
#include <unordered_map>

#define DOUBLE_MAX std::numeric_limits<double>::max()

class Rosenbrock
{
private:
      float step;
      const float rho;
      std::pair<float, float> start;
      std::string fn;
      std::vector<const char> xyz;
      const float slimit;
      std::vector<std::pair<float,float>> SD_ans;
      std::vector<std::pair<float,float>> NM_ans;
      std::vector<std::pair<std::string,float>> fvec;
      std::vector<std::vector<std::pair<std::string,float>>> gradvec;
      std::vector<std::vector<std::vector<std::pair<std::string,float>>>> hessmat;
      
public:
      Rosenbrock(float fs, const float rho, std::pair<float, float> fst, std::string fun, std::vector<const char> fin, const float fsd=0)
      : step(fs), rho(rho), start(fst), fn(fun), xyz(fin), slimit(fsd) {}
      std::unordered_map<std::string, int> Derivate(const char& );
      double computePoly(std::vector<std::pair<std::string,float>>& , std::pair<float, float>& );
      std::vector<float> Gradient(std::pair<float, float>& );
      void computeFx();
      void SD_BackTrack(std::pair<float, float>& );
      std::pair<float,float> steepestDescent();
      std::vector<std::vector<double>> Hessian(std::pair<float, float>& );
      void Inverse(std::vector<std::vector<double>>& );
      void NM_BackTrack(std::pair<float, float>& );
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

std::unordered_map<std::string, int> Rosenbrock::Derivate(const char& dx)
{
      std::unordered_map<std::string, int> fmap;
      std::string temp;
      bool nflag = false;
      bool dflag = false;
      int coeff = 1;
      int mul = 1;
      auto sz = fn.length();
      for(int i = 0; i <= sz; i++)
      {
            if(i==sz || fn[i] == '-' || fn[i] == '+')
            {
                  if(nflag)
                  {
                        coeff *= -1;
                        nflag = !nflag;
                  }
                  if(dflag)
                  {
                        if(temp.back() == dx && (fn[i-1] == dx || fn[i-2] == dx))
                        {
                              temp.pop_back();
                              if(temp.empty())
                                    temp = "1";
                        }
                        fmap[temp] += coeff;
                        dflag = !dflag;
                  }
                  coeff = 1;
                  temp.clear();
                  mul = 1;
                  if(i != sz && fn[i] == '-')
                        nflag = !nflag;
            }
            else if(fn[i] == dx)
            {
                  dflag = !dflag;
                  temp += fn[i];
            }
            else if(std::isdigit(fn[i]))
            {
                  int t = fn[i] - '0';
                  if(dflag)
                  {
                        coeff *= t;
                        if(t > 2)
                              temp += (t - 1 + '0');
                        else
                              temp.pop_back();
                  }
                  else if(temp.size())
                        temp += fn[i];
                  else
                  {
                        coeff *= mul;
                        if(coeff == 1)
                        {
                              coeff = 0;
                              mul = 10;
                        }
                        coeff += t;
                  }
            }
            else if(!std::isspace(fn[i]))
            {
                  if(dflag)
                  {
                        if(fn[i] != '^' && fn[i-1] == dx)
                              temp.pop_back();
                  }
                  temp += fn[i];
            }
      }
      return fmap;
}

double Rosenbrock::computePoly(std::vector<std::pair<std::string,float>>& vect, std::pair<float, float>& pt)
{
      double grad_sum = 0;
      
      for(auto& x : vect)
      {
            double ans = 1;
            float last = 0;
            for(auto& p : x.first)
            {
                  if(p == xyz.front())
                  {
                        ans *= pt.first;
                        last = pt.first;
                  }
                  else if(p == xyz.back())
                  {
                        ans *= pt.second;
                        last = pt.second;
                  }
                  else if (std::isdigit(p))
                  {
                        int u = p - '0';
                        while(--u)
                              ans *= last;
                  }
            }
            ans *= x.second;
            grad_sum += ans;
      }
      
      return grad_sum;
}

std::vector<float> Rosenbrock::Gradient(std::pair<float, float>& point)
{
      std::vector<float> ftemp;
      
      if(gradvec.empty())
      {
            for(const char& c : xyz)
            {
                  auto fd = Derivate(c);
                  std::vector<std::pair<std::string, float>> vec;
                  for(auto& x : fd)
                        vec.emplace_back(x.first, x.second);
                  gradvec.emplace_back(vec);
            }
      }
      
      for(auto& vect : gradvec)
      {
            ftemp.emplace_back(computePoly(vect, point));
      }
      
      return ftemp;
}

void Rosenbrock::computeFx()
{
      std::string temp;
      bool nflag = false;
      int coeff = 1;
      int mul = 1;
      auto sz = fn.size();
      
      for(int i=0; i<=sz; i++)
      {
            if(i==sz || fn[i] == '-' || fn[i] == '+')
            {
                  if(nflag)
                  {
                        coeff *= -1;
                        nflag = !nflag;
                  }
                  if(temp.empty())
                        temp = "1";
                  fvec.emplace_back(temp, coeff);
                  temp.clear();
                  coeff = 1;
                  mul = 1;
                  if(i != sz && fn[i] == '-')
                        nflag = !nflag;
            }
            else if(std::isdigit(fn[i]))
            {
                  if(temp.size())
                        temp += fn[i];
                  else
                  {
                        int t = fn[i] - '0';
                        coeff *= mul;
                        if(coeff == 1)
                        {
                              coeff = 0;
                              mul = 10;
                        }
                        coeff += t;
                  }
            }
            else if(!std::isspace(fn[i]))
            {
                  temp += fn[i];
            }
      }
}

void Rosenbrock::SD_BackTrack(std::pair<float, float>& btpt)
{
      if(fvec.empty())
      {
            computeFx();
      }
      
      float c = 0.5;
      std::vector<float> btgrad;
      std::pair<float, float> btfx;
      double btfdx;
      float mul = 1;
      double bttemp = 0;
      double left = 1;
      double right = 0;
     
      while(left > right)
      {
            step *= mul;
            btgrad = Gradient(btpt);
            btfx = {btpt.first - step * btgrad[0], btpt.second - step * btgrad[1]};
            left = computePoly(fvec, btfx);
            btfdx = computePoly(fvec, btpt);
            for(auto& bt : btgrad)
                  bttemp += bt * bt;
            bttemp =  c * step * bttemp;
            right = btfdx - bttemp;
            mul = rho;
      }
}

std::pair<float,float> Rosenbrock::steepestDescent()
{
      double sd_curr_fx = DOUBLE_MAX/100;
      double sd_prev_fx = DOUBLE_MAX;
      double sd_grad_sq = 1;
      auto sd_start = start;
      std::pair<float, float> sd_nstate = start;
      std::pair<float, float> sd_residue;
      std::vector<float> sd_grad;
      int sd_iter = 0;
      SD_ans.emplace_back(sd_nstate);
      
      while((sd_prev_fx - sd_curr_fx) > slimit)
      {
            start = sd_nstate;
            sd_prev_fx = sd_curr_fx;
            sd_grad_sq = 0;
            sd_grad = Gradient(start);
            // update learning rate until backtrack inequality is true;
            SD_BackTrack(start);
            sd_residue = {-1* sd_grad.front() * step, -1 * sd_grad.back() * step};
            sd_nstate = {start.first + sd_residue.first, start.second + sd_residue.second};
            sd_curr_fx = computePoly(fvec, sd_nstate);
            std::cout << "iteration : " << sd_iter << ",\t" << "f(x) : " << sd_curr_fx << ",\t" << "learning rate : " << step << "\n\n";
            sd_iter++;
            step = 1;
            SD_ans.emplace_back(sd_nstate);
      }
      std::cout << "\nOptimal points \t x1 : " <<  sd_nstate.first << ",\t x2 : " << sd_nstate.second << "\n\n\n\n" << std::endl;
      start = sd_start;
      return sd_nstate;
}

std::vector<std::vector<double>> Rosenbrock::Hessian(std::pair<float, float>& Hpt)
{
      std::vector<std::vector<double>> Hess;
      
      if(hessmat.empty())
      {
            std::string Hfn = fn;
            for(auto& vect : gradvec)
            {
                  fn.clear();
                  auto sz = vect.size();
                  for(int i = 0; i < sz; i++)
                  {
                        if(vect[i].second < 0)
                        {
                              if(fn.size())
                                    fn.pop_back();
                              fn += '-';
                        }
                        fn += ' ';
                        fn += std::to_string(std::abs(int(vect[i].second)));
                        if(vect[i].first != "1")
                           fn += vect[i].first;
                        if(i < sz - 1)
                              fn += " +";
                  }
                  std::vector<std::vector<std::pair<std::string, float>>> Htemp;
                  for(const char& c : xyz)
                  {
                        auto Hd = Derivate(c);
                        std::vector<std::pair<std::string, float>> Hvec;
                        for(auto& x : Hd)
                              Hvec.emplace_back(x.first, x.second);
                        Htemp.emplace_back(Hvec);
                  }
                  hessmat.emplace_back(Htemp);
            }
            fn = Hfn;
      }
      
      for(auto& hvect : hessmat)
      {
            std::vector<double> htemp;
            for(auto& hv : hvect)
            {
                  double grad_sum = computePoly(hv, Hpt);
                  htemp.emplace_back(grad_sum);
            }
            Hess.emplace_back(htemp);
      }
      
      return Hess;
}

void Rosenbrock::Inverse(std::vector<std::vector<double>>& Inv)
{
      double det = Inv.front().front() * Inv.back().back() - Inv.front().back() * Inv.back().front();
      std::swap(Inv.front().front(), Inv.back().back());
      Inv.front().front() /= det;
      Inv.back().back() /= det;
      Inv.front().back() *= (-1/det);
      Inv.back().front() *= (-1/det);
}

void Rosenbrock::NM_BackTrack(std::pair<float, float>& nmbtpt)
{
      if(fvec.empty())
      {
            computeFx();
      }
      
      float c = 0.5;
      std::vector<float> nm_grad;
      std::vector<std::vector<double>> nm_hess;
      std::pair<float, float> left_res;
      std::pair<float, float> left_arg;
      std::vector<float> nm_grad_tr;
      double right_res;
      double right_fx;
      float mul = 1;
      double left = 1;
      double right = 0;
      
      while(left > right)
      {
            step *= mul;
            nm_grad = Gradient(start);
            nm_hess = Hessian(start);
            Inverse(nm_hess);
            left_res = {(nm_hess.front().front() * nm_grad.front() + nm_hess.front().back() * nm_grad.back()),
                         (nm_hess.back().front() * nm_grad.front() + nm_hess.back().back() * nm_grad.back())};
            left_arg = {start.first - step * left_res.first, start.second - step * left_res.second};
            left = computePoly(fvec, left_arg);
            
            nm_grad_tr = nm_grad;
            right_res = nm_grad_tr.front() * left_res.first + nm_grad_tr.back() * left_res.second;
            right_fx = computePoly(fvec, nmbtpt);
            right = right_fx - c * step * right_res;
            
            mul = rho;
      }
}

std::pair<float,float> Rosenbrock::NewtonMethod()
{
      double nm_curr_fx = DOUBLE_MAX/100;
      double nm_prev_fx = DOUBLE_MAX;
      double nm_grad_sq = 1;
      std::vector<std::vector<double>> nm_hess;
      std::pair<float, float> nm_nstate = start;
      std::pair<float, float> nm_residue;
      std::vector<float> nm_gradient;
      int nm_iter = 0;
      
      NM_ans.emplace_back(nm_nstate);
      
      while(std::abs(nm_prev_fx - nm_curr_fx) > slimit)
      {
            start = nm_nstate;
            nm_prev_fx = nm_curr_fx;
            nm_grad_sq = 0;
            nm_gradient = Gradient(start);
            nm_hess = Hessian(start);
            Inverse(nm_hess);
            // update step unti backtrack inequality is true;
            NM_BackTrack(start);
            nm_residue = {step * (nm_hess.front().front() * nm_gradient.front() + nm_hess.front().back() * nm_gradient.back()),
                          step * (nm_hess.back().front() * nm_gradient.front() + nm_hess.back().back() * nm_gradient.back())};
            nm_nstate = {start.first - nm_residue.first, start.second - nm_residue.second};
            nm_curr_fx = computePoly(fvec, nm_nstate);
            std::cout << "iteration : " << nm_iter << ",\t" << "f(x) : " << nm_curr_fx << ",\t" << "learning rate : " << step << "\n\n";
            nm_iter++;
            step = 1;
            NM_ans.emplace_back(nm_nstate);
      }
      std::cout << "\nOptimal points \t x1 : " <<  nm_nstate.first << ",\t x2 : " << nm_nstate.second << "\n" << std::endl;
      
      return nm_nstate;
}

int main(int argc, const char * argv[]) {
      std::cout << "Please input the cost function in expanded form without parentheses \n\n";
      std::cout << "Example \n\t  x^2 + y^2 \n\nOtherwise, press enter to proceed with the default. The dafault is Rosenbrock Banana function. \n\n";
      std::string user_fun;
      std::string func = "100y^2 - 200x^2y + 100x^4 + x^2 - 2x + 1";
//            std::string func = "3x + x^2 - 8y + 2y^2";
      std::getline(std::cin, user_fun);
      if(user_fun.size())
            func = user_fun;
      float lr = 1;
      const float rho = 0.8;
      std::pair<float,float> start_point = {1.2, 1.2};
      std::vector<const char> xyz = {'x', 'y'};
      float stop_limit = 0.0;
      int choice = 0;
      std::cout << "\nEnter your choice for the optimization method \n\t1 - Steepest Descent \n\t2 - Newton's Method \n\t3 - Both\n\n";
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





