#include "Interpolation/chebyshev_grid.hh"
#include <stdexcept>

namespace Interpolation
{
namespace Chebyshev
{
    StandardGrid::StandardGrid(size_t p)
    {
        _p=p;
        for (size_t j=0; j<=_p; ++j) {
            _tj.push_back(cos(j*M_PI/ static_cast<double>(p)));
        }

        // beta_j = (-1)^j * (1 if j!=0, p otherwise)
        for (size_t j=0; j<=p; j++) {
            double sign = j % 2 == 0 ? +1 : -1;
            // Is j even? If yes +1, if not -1.
            double scaling = 1.;
            if (j == 0 || j == p) scaling = 0.5;

            _betaj.push_back(sign*scaling);

        }

        _Dij.resize(p+1, vector_d(p+1, 0.));
        // D is a vector of dimension (p+1) of vectors of dimension (p+1) of zeros.
        _Dij[0][0] = (2. * p * p+1)/6.;
        _Dij[p][p] = -_Dij[0][0];

        for(size_t j = 1; j < p; j++){
            _Dij[j][j] = -0.5 * _tj[j] / (1. - pow(_tj[j], 2));
            // NB: Do not use ^2!! It means another thing. Use pow() or just multiply by itself.
        }

        for (size_t i=0; i<= p; i++) {
            for (size_t j=0; j<= p; j++) {
                if ( j== 1) continue;
                _Dij[i][j] = - (_betaj[i] / _betaj[j]) / (_tj[i] - _tj[j]);
            }
        }

    }

    /*
    double StandardGrid::poli_weight(double t, size_t j) const
    {
        double den = 0.;
        for(int i=0; i< this->_betaj.size(); i++){
            den += this->_betaj[i] / (t - this->_tj[i]);
        }
        
        double b_i = den * this->_betaj[j] / (t - this->_tj[j]);
        return b_i;
    }
    */

/*
  double StandardGrid::poli_weight(double t, size_t j, double den) const
    {
        double b_i = den * this->_betaj[j] / (t - this->_tj[j]);
        return b_i;
    }
*/
    double StandardGrid::interpolate(double t, const vector_d &fj, size_t start, size_t end) const
    // the "&" means "pass by reference", i.e. do not copy the vector, just pass a reference
    {
        if(t<-1 | t > 1){
            throw std::domain_error("StandardGrid::interpolate: t must be in [-1, 1]");
        }
        if(end - start != _p){
            throw std::domain_error("StandardGrid::interpolate: end - start must be equal to p");
        }
        /*
        double res = 0.;
        for(size_t i = 0, i<=_p; i++){
            res += poli_weight(t, i) * fj[i + start];
        }
        return res;
        */

        double den = 0.;
        for(size_t j=0; j<= _p; j++){
            den += _betaj[j] / (t - _tj[j]);
        }
        double res = 0.;
        for(size_t i = 0;i<=_p; i++){
            res += poli_weight(t, i, den) * fj[i + start];
        }
        return res;

    }

    double StandardGrid::poli_weight(double t, size_t j, double den) const
    {
        if (std::abs(t - _tj[j]) < 1e-15) 
            // If t is very close to _tj[j], return the limit value of the weight, which is 1.
            return 1.;
            double res = 0.;
            res = _betaj[j] / (t - _tj[j]) /den;
            return res;
    }

    double StandardGrid::poli_weight(double t, size_t j) const
    {
       if (std::abs(t - _tj[j]) < 1e-15) 
            // If t is very close to _tj[j], return the limit value of the weight, which is 1.
            return 1.;
        double den = 0.;
        for(size_t j=0; j<= _p; j++){
            if (std::abs(t - _tj[j]) < 1e-15) 
                // If t is very close to _tj[j], return the limit value of the weight, which is 1.
                return 0.;
            den += _betaj[j] / (t - _tj[j]);
        }
        double res = 0.;
        res = _betaj[j] / (t - _tj[j]) /den;
        return res;
    }

    vector_d StandardGrid::discretize(const std::function<double(double)> &f) const
    {
        vector_d fj(_p+1, 0.);
        for (size_t i=0; i<= _p;i++){
            fj[i]=f(_tj[i]);
        }
        return fj;
    }

    double StandardGrid::poli_weight_der(double t, size_t j) const{

        double l=0;
        for(size_t i=0; i <= _p; i++){
            l += poli_weight(t, i)*_Dij[i][j];
        }
        return l;
    }

    double StandardGrid::interpolate_der(double t, const vector_d &fj, size_t start, size_t end) const
    {
        if(t<-1 | t > 1){
            throw std::domain_error("StandardGrid::interpolate: t must be in [-1, 1]");
        }
        if(end - start != _p){
            throw std::domain_error("StandardGrid::interpolate: end - start must be equal to p");
        }

        double p = 0.;
        for (size_t i= start, j=0; i<= end; i++, j++){
            p+= fj[i] * poli_weight_der(t, j);

        }
        return p;

    }

    double StandardGrid::poli_weight_der(double t, size_t j, double den) const
    {
        double l=0;
        for(size_t i=0; i <= _p; i++){
            l += poli_weight(t, i, den)*_Dij[i][j];
        }
        return l;
    }

    void StandardGrid::apply_D(vector_d &fj, size_t start, size_t end) const
    // Guarda la soluzione da Rodini
    {
       if (end - start != _p) {
      throw std::invalid_argument("[StandardGrid::apply_D]: cannot apply "
                                  "derivative matrix to partial vector.");
        }

        vector_d f_tilde = fj;
        for(size_t j= start; j<= end; j++){
            for (size_t i= start; i<= end; i++){
                fj[i] += f_tilde[i] *_Dij[j][i]  ;
            } 
            }
   } 

   vector_d StandardGrid::discretize(const std::function<double(double)> &fnc) const
{
   vector_d fj(_p + 1, 0.);
   vector_d result(_p + 1, 0.);
   for (size_t i = 0; i <= _p; i++) {
      fj[i] = fnc(_tj[i]);
      result[i] = fnc(_tj[i]);
   }
   return fj;
   return result;
}

} // namespace Cebyshev
} // namespace Interpolation