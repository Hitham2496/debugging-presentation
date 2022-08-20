#include <iostream>
#include <iomanip>
#include <vector>
#include <math.h>
#include <memory>

//! This class contains a bare-bones implementation of four-momentum
//! for use in basic high energy physics (HEP) calculations.
class FourMomentum {
public:

  FourMomentum(
    double E,
    double pX,
    double pY,
    double pZ
  ):
    p({E, pX, pY, pZ})
  {}

  FourMomentum(
    const std::vector<double>& momentum
  ):
    p(momentum)
  {}

  double E() const{
    return p[0];
  }
  double px() const{
    return p[1];
  }
  double py() const{
    return p[2];
  }
  double pz() const{
    return p[3];
  }
  double pperp2() const{
    return p[1]*p[1] + p[2]*p[2];
  }
  double pperp() const{
    return sqrt(pperp2());
  }
  double rap() const{
    return 0.5*log((p[0]+p[3])/(p[0]-p[3]));
  }
  double phi() const{
    return atan2(p[2],p[1]);
  }
  double m2() const{
    return p[0]*p[0] - p[1]*p[1] - p[2]*p[2] - p[3]*p[3];
  }
  double m() const{
    return sqrt(m2());
  }
  std::vector<double> P() const{
    return p;
  }

  //! Overloaded += operator:
  //!   adds a four momentum to the current (p = p + second).
  FourMomentum& operator+=(const FourMomentum& second){
    for (int j=0; j<p.size(); ++j){
      	  this->p[j] += second.P()[j];
    }
    return *this;
  }

  //! Overloaded -= operator:
  //!   subtracts a four momentum from the current (p = p - second).
  FourMomentum& operator-=(const FourMomentum& second){
    for (int j=0; j<p.size(); ++j){
      	  this->p[j] -= second.P()[j];
    }
    return *this;
  }

  ~FourMomentum() {};

private:
	
  std::vector<double> p;

};


//! Overloaded insertion operator to allow printing of FourMomentum
//! objects - this will be useful during the exercise.
std::ostream& operator<<(std::ostream & os, FourMomentum const & f){
  os << "Four Momentum with components\n";
  os << std::setprecision(9) << "E = " << f.P()[0]
                             << " px = " << f.P()[1]
      		             << " py = " << f.P()[2]
      		             << " pz = " << f.P()[3]
			     << " || (mass)^2 = " << f.m2() << std::endl;
  return os;
}

//! Returns the product of two four vectors under the Minkowski metric.
double dot(const FourMomentum& pa, const FourMomentum& pb){
  return pa.E()*pb.E() - pa.px()*pb.px()-pa.py()*pb.py()-pa.pz()*pb.pz();
}


//! Transformation of input momenta before performing our calculation.
FourMomentum pre_calc_transform(const FourMomentum& a,
  const FourMomentum& b, const FourMomentum& c){

  FourMomentum x = b;
  FourMomentum y = x;

  // x = b + c
  x+=c;

  // y = b - c
  y-=c;

  return FourMomentum(a.E()+x.E(),
		      a.px()+x.pperp()*cos(x.phi())+y.pperp()*sin(y.phi()),
		      a.py()+x.pperp()*sin(x.phi())+y.pperp()*cos(y.phi()),
		      a.pz()+y.pz());

}

//! Performs our calculation, the product of the logarithm of the
//! transformed mass and the transformed transver momentum (normalised
//! by scale q).
double do_calculation(const FourMomentum& transformed, double q2){

  double lSoft = log(transformed.pperp2()/q2);
  double lHE = log(transformed.m2()/q2);

  return lSoft*lHE;

}


//! Performs our HEP calculation:
//!  1. transform momenta a,b,c into t according to pre_calc_transform
//!  2. find some numerical factor l according to do_calculation
//!  3. print the final result = l*dot(a,b+c)
int main() {

  FourMomentum a(200.,0.,0.,200.);
  FourMomentum b(90.,30.,30.,2000.);
  FourMomentum c(45.,15.,20.,1000.);

  double q2 = 100.;

  std::cout << "Performing a horrible calculation with momenta:\n"
	    << a << "\n" << b << "\n" << c << std::endl;

  FourMomentum t = pre_calc_transform(a, b, c);

  double log_product = do_calculation(t, q2);

  b+=c;

  double result = log_product*dot(a,b);

  // Final result is product of logs multiplied by a.(b+c)
  std::cout << "Answer is : " << result << std::endl;
		   
  return 0;
}
