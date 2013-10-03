#ifndef _REG_CV_H
#define _REG_CV_H

#include <tulip/TulipPluginHeaders.h>
#include <tulip/StringCollection.h>
#include <stdexcept>
#include <vector>
#include <utility>

class Reg_ROF:public tlp::Algorithm {
private:
	unsigned int iter_max;
	double lambda;
	tlp::DoubleProperty *seed, *result, *fn, *fnp1, *beta, *w, *f0;
	tlp::BooleanProperty *segmentation_result;
	std::string fn_name, fnp1_name, beta_name;

	unsigned int export_interval;
	std::string export_directory;

public:
	PLUGININFORMATIONS("Rudin-Osher-Fatemi Regularization", "Cyrille FAUCHEUX", "2013-10-02", "", "1.0", "Regularization")

	Reg_ROF(const tlp::PluginContext *context);
	~Reg_ROF() {}
	bool check(std::string &);
	bool run();

private:
	std::pair<double, double> computeFnMinMax();
	void fnToSelection();
	void exportIntermediateResult(const int i);
};

PLUGIN(Reg_ROF)

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
