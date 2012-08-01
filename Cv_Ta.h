#ifndef _CV_TA_H
#define _CV_TA_H

#include <tulip/TulipPlugin.h>
#include <utility>

class Cv_Ta:public tlp::Algorithm {
private:
	unsigned int iter_max;
	double lambda1, lambda2;
	tlp::DoubleProperty *fn, *fnp1, *beta, *w;
	tlp::DoubleVectorProperty *f0;
	tlp::BooleanProperty *roi;
	std::pair< std::vector< double >, std::vector< double > > in_out_means;
	unsigned int f0_size;

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
