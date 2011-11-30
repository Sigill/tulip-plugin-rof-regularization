#ifndef _CV_TA_H
#define _CV_TA_H

#include <tulip/TulipPlugin.h>
#include <utility>

class Cv_Ta:public tlp::Algorithm {

public:
  Cv_Ta(const tlp::AlgorithmContext& context);
  ~Cv_Ta();
  bool run();
  bool check(std::string &);

private:
	void fnToSelection();
	void exportSelection(const int i);
	void computeMeanValues(std::pair<double, double>* c);	

};
#endif
