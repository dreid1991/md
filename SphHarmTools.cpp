#include "SphHarmTools.h"
#include <iostream>
double SphHarmTools::PLM(int l, int m, double x)
{
  double fact,pll,pmm,pmmp1,somx2;
  int i,ll;
  pll = 0.0;
  if (m < 0 || m > l || fabs(x) > 1.0) {
	  cerr << "impossible combination of l and m" << "\n";
	}
  pmm=1.0;
  if (m > 0)
  {
    somx2=sqrt((1.0-x)*(1.0+x));
    fact=1.0;
    for (i=1;i<=m;i++)
    {
      pmm *= -fact*somx2;
      fact += 2.0;
    }
  }
  if (l == m)
    return pmm;
  else
  {
    pmmp1=x*(2*m+1)*pmm;
    if (l == (m+1))
      return pmmp1;
    else
    {
      for (ll=m+2;ll<=l;ll++)
        {
          pll=(x*(2*ll-1)*pmmp1-(ll+m-1)*pmm)/(ll-m);
          pmm=pmmp1;
          pmmp1=pll;
        }
        return pll;
     }
  }
}

int SphHarmTools::fact(int x) {
	if (x==0) {
		return 1;
	} else {
		int prod = 1;
		for (int i=2; i<=x; i++) {
			prod *= i;
		}
		return prod;
	}
}

complex<double> SphHarmTools::evalSphHarm(int l, int m, Vector &v, int i) {
	double r = v.len();
	double PI = 3.141592654;
	double theta = acos(v.z / r);
	double phi = atan2(v.y, v.x);

	if (phi < 0) {
		phi += 2 * PI;
	}
	int mAbs = fabs(m);
	double legendreVal = PLM(l, fabs(m), cos(theta));
	double coef = pow(-1., m) * sqrt((2 * l + 1) * fact(l - mAbs) / (4 * PI * fact(l + mAbs))) * legendreVal;
	return complex<double>(coef * cos(m * phi), coef * sin(m * phi));
}


vector<complex<double> > SphHarmTools::calcQlms(Atom *a, int l) {
	vector<complex<double> > qlms;
	complex<double> blank;
	int numNeighbors = a->neighbors.size();
	for (int m=-l; m<=l; m++) {
		complex<double> sphHarmSum;//make sure is zero by default

		for (int i=0; i<numNeighbors; i++) {
			Vector vToN = a->vToNeighbor(i);
			sphHarmSum += evalSphHarm(l, m, vToN, i);
		}

		if (numNeighbors) {
			qlms.push_back(sphHarmSum / (double) numNeighbors);
		} else {

			qlms.push_back(blank);
		}
	}
	return qlms;
}

double SphHarmTools::calcQl(Atom *a, int l) {
	double PI = 3.141592654;
	double sumQlmSqrs = 0;
	vector<complex<double> > qlms = calcQlms(a, l);
	for (unsigned int i=0; i<qlms.size(); i++) {
		sumQlmSqrs += qlms[i].real() * qlms[i].real() + qlms[i].imag() * qlms[i].imag();
	}
	return sqrt(4 * PI / (2 * l + 1) * sumQlmSqrs);
}

complex<double> SphHarmTools::calcQlmBar(Atom *a, int l, int m) {
	complex<double> sumQlms;
	for (unsigned int i=0; i<a->neighbors.size(); i++) {
		sumQlms += a->neighbors[i]->qms[l][m + l];
	}
	sumQlms += a->qms[l][m + l];
	return sumQlms / (double) (a->neighbors.size() + 1);
}

void SphHarmTools::calcQlAvg(vector<Atom *> &atoms, int l) {
	double PI = 3.141592654;
	for (unsigned int i=0; i<atoms.size(); i++)	{
		atoms[i]->qms[l] = calcQlms(atoms[i], l);	
	}
	vector<vector<complex<double> > > qlmBars;
	for (unsigned int i=0; i<atoms.size(); i++) {
		vector<complex<double> > qlmBarNext;
		for (int m=-l; m<=l; m++) {
			qlmBarNext.push_back(calcQlmBar(atoms[i], l, m));
		}
		qlmBars.push_back(qlmBarNext);
	}
	for (unsigned int i=0; i<atoms.size(); i++) {
		double sumQlmBarSqrs;
		for (int m=-l; m<=l; m++) {
			sumQlmBarSqrs += pow(qlmBars[i][m + l].real(), 2) + pow(qlmBars[i][m + l].imag(), 2);
		}
		atoms[i]->qs[l] = sqrt(4 * PI / (2 * l + 1) * sumQlmBarSqrs);
	}
}
