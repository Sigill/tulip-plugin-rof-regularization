#ifndef _CV_TA_H
#define _CV_TA_H

#include <tulip/TulipPlugin.h>
#include <utility>
#include <stdexcept>

class Cv_Ta:public tlp::Algorithm {
private:
	unsigned int iter_max;
	double lambda1, lambda2;
	tlp::DoubleProperty *fn, *fnp1, *beta, *w;
	tlp::DoubleVectorProperty *f0;
	tlp::BooleanProperty *seed, *roi;

	unsigned int export_interval;
	std::string export_directory;

	std::pair< std::vector< double >, std::vector< double > > in_out_means;
	unsigned int f0_size;

public:
	Cv_Ta(const tlp::AlgorithmContext& context);
	~Cv_Ta() {}
	bool run();
	bool check(std::string &);

private:
	void fnToSelection();
	void exportSelection(const int i);
	void computeMeanValues();
};

class export_exception : public std::runtime_error {
public:
	export_exception( const std::string &err ) : std::runtime_error( err ) {}
};

#endif
