#ifndef _REG_CV_H
#define _REG_CV_H

#include <tulip/TulipPluginHeaders.h>
#include <tulip/StringCollection.h>
#include <stdexcept>
#include <vector>
#include <utility>

class Reg_CV:public tlp::Algorithm {
private:
	unsigned int iter_max;
	double lambda1, lambda2;
	tlp::PropertyInterface *f0;
	tlp::DoubleProperty *seed, *result, *fn, *fnp1, *beta, *w, *f0_scalar;
	tlp::DoubleVectorProperty *f0_vector;
	tlp::BooleanProperty *segmentation_result;
	std::string fn_name, fnp1_name, beta_name;

	unsigned int export_interval;
	std::string export_directory;

	std::pair< std::vector< double >, std::vector< double > > in_out_means;
	unsigned int f0_size;

public:
	PLUGININFORMATIONS("ChanVese Regularization", "Cyrille FAUCHEUX", "2013-08-11", "", "1.0", "Regularization")

	Reg_CV(const tlp::PluginContext *context);
	~Reg_CV() {}
	bool check(std::string &);
	bool run();

private:
	std::pair<double, double> computeFnMinMax();
	void fnToSelection();
	void exportIntermediateResult(const int i);
	void computeMeanValues();
};

PLUGIN(Reg_CV)

#define CHECK_PROP_PROVIDED(PROP, STOR) \
	do { \
		if(!dataSet->get(PROP, STOR)) \
			throw std::runtime_error(std::string("No \"") + PROP + "\" provided"); \
	} while(0)

class export_exception : public std::runtime_error {
public:
	export_exception( const std::string &err ) : std::runtime_error( err ) {}
};

#endif
