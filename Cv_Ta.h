#ifndef _CV_TA_H
#define _CV_TA_H

#include <tulip/TulipPlugin.h>
#include <utility>

class Cv_Ta:public tlp::Algorithm {
private:
	unsigned int iter_max;
	double lambda;
	tlp::DoubleProperty *f0, *fn, *fnp1, *beta;
	std::pair< double, double > in_out_means;

public:
  Cv_Ta(const tlp::AlgorithmContext& context);
  ~Cv_Ta() {};
  bool run();
  bool check(std::string &) { return true; }

private:
	void fnToSelection();
	void exportSelection(const int i);
	void computeMeanValues();	

};
#endif
