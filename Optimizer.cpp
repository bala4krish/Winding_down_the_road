//
//  Optimizer.cpp
//  Numerical_Optimization
//
//  Created by Balaji Krishnamurthy on 12/20/20.
//  Copyright © 2020 Balaji Krishnamurthy. All rights reserved.
//

#include "Optimizer.hpp"


// Function to compute derivative, used for Gradient and Hessian
std::unordered_map<std::string, int> Rosenbrock::Derivate(const char& dx)
{
      std::unordered_map<std::string, int> fmap;
      std::string temp;
      bool nflag = false;
      bool dflag = false;
      int coeff = 1;
      int mul = 1;
      const auto sz = this->fn.length();
      
      for(int i = 0; i <= sz; i++)
      {
            if(i == sz || this->fn[i] == '-' || this->fn[i] == '+')
            {
                  if(nflag)
                  {
                        coeff *= -1;
                        nflag = !nflag;
                  }
                  if(dflag)
                  {
                        if(temp.back() == dx && (this->fn[i-1] == dx || this->fn[i-2] == dx))
                        {
                              temp.pop_back();
                              if(temp.empty())
                              {
                                    temp = "1";
                              }
                        }
                        fmap[temp] += coeff;
                        dflag = !dflag;
                  }
                  coeff = 1;
                  temp.clear();
                  mul = 1;
                  if(i != sz && this->fn[i] == '-')
                  {
                        nflag = !nflag;
                  }
            }
            else if(this->fn[i] == dx)
            {
                  dflag = !dflag;
                  temp += this->fn[i];
            }
            else if(std::isdigit(this->fn[i]))
            {
                  int t = this->fn[i] - '0';
                  if(dflag)
                  {
                        coeff *= t;
                        if(t > 2)
                        {
                              temp += (t - 1 + '0');
                        }
                        else
                        {
                              temp.pop_back();
                        }
                  }
                  else if(temp.size())
                  {
                        temp += this->fn[i];
                  }
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
            else if(!std::isspace(this->fn[i]))
            {
                  if(dflag)
                  {
                        if(this->fn[i] != '^' && this->fn[i-1] == dx)
                        {
                              temp.pop_back();
                        }
                  }
                  temp += this->fn[i];
            }
      }
      
      return fmap;
}

// This function computes value of a polynomial at given points
double Rosenbrock::computePoly(const std::vector<std::pair<std::string,float>>& vect, const std::pair<float, float>& pt)
{
      double grad_sum = 0;
      
      for(const auto& x : vect)
      {
            double ans = 1;
            float last = 0;
            for(const auto& p : x.first)
            {
                  if(p == this->xyz.front())
                  {
                        ans *= pt.first;
                        last = pt.first;
                  }
                  else if(p == this->xyz.back())
                  {
                        ans *= pt.second;
                        last = pt.second;
                  }
                  else if (std::isdigit(p))
                  {
                        int u = p - '0';
                        while(--u)
                        {
                              ans *= last;
                        }
                  }
            }
            ans *= x.second;
            grad_sum += ans;
      }
      
      return grad_sum;
}

// Function to compute Gradient
std::vector<float> Rosenbrock::Gradient(const std::pair<float, float>& point)
{
      std::vector<float> ftemp;
      
      if(this->gradvec.empty())
      {
            for(const char& c : this->xyz)
            {
                  const auto fd = this->Derivate(c);
                  std::vector<std::pair<std::string, float>> vec;
                  for(const auto& x : fd)
                  {
                        vec.emplace_back(x.first, x.second);
                  }
                  this->gradvec.emplace_back(vec);
            }
      }
      
      for(const auto& vect : this->gradvec)
      {
            ftemp.emplace_back(this->computePoly(vect, point));
      }
      
      return ftemp;
}

// Function to compute f(x)
void Rosenbrock::computeFx()
{
      std::string temp;
      bool nflag = false;
      int coeff = 1;
      int mul = 1;
      const auto sz = this->fn.size();
      
      for(int i = 0; i <= sz; i++)
      {
            if(i == sz || this->fn[i] == '-' || this->fn[i] == '+')
            {
                  if(nflag)
                  {
                        coeff *= -1;
                        nflag = !nflag;
                  }
                  if(temp.empty())
                  {
                        temp = "1";
                  }
                  this->fvec.emplace_back(temp, coeff);
                  temp.clear();
                  coeff = 1;
                  mul = 1;
                  if(i != sz && this->fn[i] == '-')
                  {
                        nflag = !nflag;
                  }
            }
            else if(std::isdigit(this->fn[i]))
            {
                  if(temp.size())
                  {
                        temp += this->fn[i];
                  }
                  else
                  {
                        int t = this->fn[i] - '0';
                        coeff *= mul;
                        if(coeff == 1)
                        {
                              coeff = 0;
                              mul = 10;
                        }
                        coeff += t;
                  }
            }
            else if(!std::isspace(this->fn[i]))
            {
                  temp += this->fn[i];
            }
      }
}

// Function to bactrack for Steepest Descent
void Rosenbrock::SD_BackTrack(const std::pair<float, float>& btpt)
{
      if(this->fvec.empty())
      {
            this->computeFx();
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
            this->step *= mul;
            btgrad = this->Gradient(btpt);
            btfx = {btpt.first - this->step * btgrad[0], btpt.second - this->step * btgrad[1]};
            left = this->computePoly(this->fvec, btfx);
            btfdx = this->computePoly(this->fvec, btpt);
            for(const auto& bt : btgrad)
            {
                  bttemp += bt * bt;
            }
            bttemp =  c * this->step * bttemp;
            right = btfdx - bttemp;
            mul = this->rho;
      }
}

// Function for the Steepest Descent algorithm
std::pair<float,float> Rosenbrock::steepestDescent()
{
      double sd_curr_fx = DOUBLE_MAX/100;
      double sd_prev_fx = DOUBLE_MAX;
      const auto sd_start = this->start;
      std::pair<float, float> sd_nstate = this->start;
      std::pair<float, float> sd_residue;
      std::vector<float> sd_grad;
      int sd_iter = 0;
      this->SD_ans.emplace_back(sd_nstate);
      
      while((sd_prev_fx - sd_curr_fx) > slimit)
      {
            this->start = sd_nstate;
            sd_prev_fx = sd_curr_fx;
            sd_grad = this->Gradient(this->start);
            // update learning rate until backtrack inequality is true;
            this->SD_BackTrack(this->start);
            sd_residue = {-1* sd_grad.front() * this->step, -1 * sd_grad.back() * this->step};
            sd_nstate = {this->start.first + sd_residue.first, this->start.second + sd_residue.second};
            sd_curr_fx = this->computePoly(this->fvec, sd_nstate);
            std::cout << "iteration : " << sd_iter << ",\t" << "f(x) : " << sd_curr_fx << ",\t" << "learning rate : " << this->step << "\n\n";
            sd_iter++;
            this->step = 1;
            this->SD_ans.emplace_back(sd_nstate);
      }
      std::cout << "\nOptimal points \t x1 : " <<  sd_nstate.first << ",\t x2 : " << sd_nstate.second << "\n\n\n\n" << std::endl;
      this->start = sd_start;
      
      return sd_nstate;
}

// Function to compute Hessian matrix
std::vector<std::vector<double>> Rosenbrock::Hessian(const std::pair<float, float>& Hpt)
{
      std::vector<std::vector<double>> Hess;
      
      if(this->hessmat.empty())
      {
            const std::string Hfn = this->fn;
            for(const auto& vect : this->gradvec)
            {
                  this->fn.clear();
                  const auto sz = vect.size();
                  for(int i = 0; i < sz; i++)
                  {
                        if(vect[i].second < 0)
                        {
                              if(this->fn.size())
                              {
                                    this->fn.pop_back();
                              }
                              this->fn += '-';
                        }
                        this->fn += ' ';
                        this->fn += std::to_string(std::abs(int(vect[i].second)));
                        if(vect[i].first != "1")
                        {
                           this->fn += vect[i].first;
                        }
                        if(i < sz - 1)
                        {
                              this->fn += " +";
                        }
                  }
                  
                  std::vector<std::vector<std::pair<std::string, float>>> Htemp;
                  
                  for(const char& c : this->xyz)
                  {
                        const auto Hd = this->Derivate(c);
                        std::vector<std::pair<std::string, float>> Hvec;
                        for(const auto& x : Hd)
                        {
                              Hvec.emplace_back(x.first, x.second);
                        }
                        Htemp.emplace_back(Hvec);
                  }
                  this->hessmat.emplace_back(Htemp);
            }
            this->fn = Hfn;
      }
      
      for(const auto& hvect : this->hessmat)
      {
            std::vector<double> htemp;
            for(const auto& hv : hvect)
            {
                  double grad_sum = this->computePoly(hv, Hpt);
                  htemp.emplace_back(grad_sum);
            }
            Hess.emplace_back(htemp);
      }
      
      return Hess;
}

// Function to compute Inverse of the Hessian matrix
void Rosenbrock::Inverse(std::vector<std::vector<double>>& Inv)
{
      const double det = Inv.front().front() * Inv.back().back() - Inv.front().back() * Inv.back().front();
      std::swap(Inv.front().front(), Inv.back().back());
      Inv.front().front() /= det;
      Inv.back().back() /= det;
      Inv.front().back() *= (-1/det);
      Inv.back().front() *= (-1/det);
}

// Function to bactrack for Newton's Method
void Rosenbrock::NM_BackTrack(const std::pair<float, float>& nmbtpt)
{
      if(this->fvec.empty())
      {
            this->computeFx();
      }
      
      const float c = 0.5;
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
            this->step *= mul;
            nm_grad = this->Gradient(this->start);
            nm_hess = this->Hessian(this->start);
            this->Inverse(nm_hess);
            left_res = {(nm_hess.front().front() * nm_grad.front() + nm_hess.front().back() * nm_grad.back()),
                         (nm_hess.back().front() * nm_grad.front() + nm_hess.back().back() * nm_grad.back())};
            left_arg = {this->start.first - this->step * left_res.first, this->start.second - this->step * left_res.second};
            left = this->computePoly(this->fvec, left_arg);
            
            nm_grad_tr = nm_grad;
            right_res = nm_grad_tr.front() * left_res.first + nm_grad_tr.back() * left_res.second;
            right_fx = this->computePoly(this->fvec, nmbtpt);
            right = right_fx - c * this->step * right_res;
            
            mul = this->rho;
      }
}

// Function to implement Newton's algorithm
std::pair<float,float> Rosenbrock::NewtonMethod()
{
      double nm_curr_fx = DOUBLE_MAX/100;
      double nm_prev_fx = DOUBLE_MAX;
      std::vector<std::vector<double>> nm_hess;
      std::pair<float, float> nm_nstate = this->start;
      std::pair<float, float> nm_residue;
      std::vector<float> nm_gradient;
      int nm_iter = 0;
      
      this->NM_ans.emplace_back(nm_nstate);
      
      while(std::abs(nm_prev_fx - nm_curr_fx) > this->slimit)
      {
            this->start = nm_nstate;
            nm_prev_fx = nm_curr_fx;
            nm_gradient = this->Gradient(this->start);
            nm_hess = this->Hessian(this->start);
            this->Inverse(nm_hess);
            // update step unti backtrack inequality is true;
            this->NM_BackTrack(start);
            nm_residue = {this->step * (nm_hess.front().front() * nm_gradient.front() + nm_hess.front().back() * nm_gradient.back()),
                          this->step * (nm_hess.back().front() * nm_gradient.front() + nm_hess.back().back() * nm_gradient.back())};
            nm_nstate = {this->start.first - nm_residue.first, this->start.second - nm_residue.second};
            nm_curr_fx = this->computePoly(this->fvec, nm_nstate);
            std::cout << "iteration : " << nm_iter << ",\t" << "f(x) : " << nm_curr_fx << ",\t" << "learning rate : " << this->step << "\n\n";
            nm_iter++;
            this->step = 1;
            this->NM_ans.emplace_back(nm_nstate);
      }
      
      std::cout << "\nOptimal points \t x1 : " <<  nm_nstate.first << ",\t x2 : " << nm_nstate.second << "\n" << std::endl;
      
      return nm_nstate;
}

/* Optimizer_hpp */
