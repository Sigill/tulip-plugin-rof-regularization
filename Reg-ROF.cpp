#include "Reg-ROF.h"
#include <cmath>
#include <sstream>
#include <iostream>
#include <iomanip>
#include <algorithm>
#include <cfloat>

#include <QFileInfo>
#include <QDir>

using namespace std;
using namespace tlp;

template <class T>
std::string to_string(T t, std::ios_base & (*f)(std::ios_base&))
{
	std::ostringstream oss;
	oss << f << t;
	return oss.str();
}

std::string random_string(const size_t len) {
	static const char alphanum[] =
		"0123456789"
		"ABCDEFGHIJKLMNOPQRSTUVWXYZ"
		"abcdefghijklmnopqrstuvwxyz";

	std::string result(len, '-');

	for (size_t i = 0; i < len; ++i) {
		result[i] = alphanum[rand() % 61];
	}

	return result;
}

namespace {
	const char * paramHelp[] = {
		// 0 Seed
		HTML_HELP_OPEN() \
			HTML_HELP_DEF( "type", "DoubleProperty" ) \
			HTML_HELP_BODY() \
			"The double property to regularize." \
			HTML_HELP_CLOSE(),

		// 1 Result
		HTML_HELP_OPEN() \
			HTML_HELP_DEF( "type", "DoubleProperty" ) \
			HTML_HELP_BODY() \
			"The double property where the resulting double property will be stored." \
			HTML_HELP_CLOSE(),

		// 2 Data
		HTML_HELP_OPEN() \
			HTML_HELP_DEF( "type", "DoubleProperty" ) \
			HTML_HELP_BODY() \
			"Specify which property holds the data associated to each nodes." \
			HTML_HELP_CLOSE(),

		// 3 Similarity measure
		HTML_HELP_OPEN() \
			HTML_HELP_DEF( "type", "DoubleProperty" ) \
			HTML_HELP_BODY() \
			"Specify which property holds the weight associated to each edge." \
			HTML_HELP_CLOSE(),

		// 4 Number of iterations
		HTML_HELP_OPEN() \
			HTML_HELP_DEF( "type", "Unsigned int" ) \
			HTML_HELP_BODY() \
			"Defines the number of iterations the algorithm will perform." \
			HTML_HELP_CLOSE(),

		// 5 Lambda
		HTML_HELP_OPEN() \
			HTML_HELP_DEF( "type", "Double" ) \
			HTML_HELP_BODY() \
			"Lambda parameter." \
			HTML_HELP_CLOSE(),

		// 6 Export interval
		HTML_HELP_OPEN() \
			HTML_HELP_DEF( "type", "Unsigned int" ) \
			HTML_HELP_BODY() \
			"Specify at which interval the processed graph must be exported (0 to disable)." \
			HTML_HELP_CLOSE(),

		// 7 Export directory
		HTML_HELP_OPEN() \
			HTML_HELP_DEF( "type", "Directory pathname" ) \
			HTML_HELP_BODY() \
			"This parameter is used to specify where the processed graph must be exported." \
			HTML_HELP_CLOSE(),

		// 8 Segmentation result
		HTML_HELP_OPEN() \
			HTML_HELP_DEF( "type", "BooleanProperty" ) \
			HTML_HELP_BODY() \
			"Indicates in which BooleanProperty the segmentation result must be saved." \
			HTML_HELP_CLOSE()
	};
}

//======================================================
Reg_ROF::Reg_ROF(const tlp::PluginContext *context):
	Algorithm(context),
	seed(NULL), result(NULL), w(NULL), iter_max(0), lambda(0.0), export_interval(0), export_directory(),
	fn(NULL), fnp1(NULL), beta(NULL), f0(NULL), segmentation_result(NULL),
	fn_name(), fnp1_name(), beta_name()
{
	addInParameter< DoubleProperty >      ("seed",                  paramHelp[0], "viewMetric");
	addInParameter< DoubleProperty >      ("result",                paramHelp[1], "viewMetric");
	addInParameter< BooleanProperty >     ("segmentation result",   paramHelp[8], "viewSelection");
	addInParameter< DoubleProperty >      ("data",                  paramHelp[2], "viewMetric");
	addInParameter< DoubleProperty >      ("similarity measure",    paramHelp[3], "viewMetric");
	addInParameter< unsigned int >        ("number of iterations",  paramHelp[4], "100");
	addInParameter< double >              ("lambda",                paramHelp[5], "1");
	addInParameter< unsigned int >        ("export interval",       paramHelp[6], "0");
	addInParameter< string >              ("dir::export directory", paramHelp[7], "");
}

bool Reg_ROF::check(std::string &err)
{
	try {
		if(dataSet == NULL)
			throw std::runtime_error("No dataset provided.");

		CHECK_PROP_PROVIDED("seed", this->seed);

		CHECK_PROP_PROVIDED("result", this->result);

		CHECK_PROP_PROVIDED("segmentation result", this->segmentation_result);

		CHECK_PROP_PROVIDED("data", this->f0);

		CHECK_PROP_PROVIDED("similarity measure", this->w);

		CHECK_PROP_PROVIDED("number of iterations", this->iter_max);

		CHECK_PROP_PROVIDED("lambda", this->lambda);

		CHECK_PROP_PROVIDED("export interval", this->export_interval);

		if(this->iter_max < 0) {
			std::ostringstream m;
			m << "Invalid number of iterations: " << this->iter_max;
			throw std::runtime_error(m.str());
		}

		// No need for an export directory if export is disabled
		if(this->export_interval > 0) {
			CHECK_PROP_PROVIDED("dir::export directory", this->export_directory);
		}

		/*
		 * Find an unused property name for fn, fnp1 & beta
		 */
		do {
			this->fn_name = random_string(6);
		} while(graph->existLocalProperty(this->fn_name));
		do {
			this->fnp1_name = random_string(6);
		} while(graph->existLocalProperty(this->fnp1_name));
		do {
			this->beta_name = random_string(6);
		} while(graph->existLocalProperty(this->beta_name));

		this->fn   = graph->getLocalProperty< DoubleProperty >(this->fn_name);
		this->fnp1 = graph->getLocalProperty< DoubleProperty >(this->fnp1_name);
		this->beta = graph->getLocalProperty< DoubleProperty >(this->beta_name);

		std::cerr << "Processing graph " << graph->getName() << std::endl
		          << "Number of iterations: " << this->iter_max << std::endl
		          << "Lambda: " << this->lambda << std::endl
		          << "Export interval: " << this->export_interval << std::endl
		          << "Export directory: " << this->export_directory << std::endl;
	} catch (std::runtime_error &ex) {
		err.assign(ex.what());
		return false;
	}

	return true;
}

//======================================================
bool Reg_ROF::run() {
	this->fn->setAllNodeValue(0.0);
	this->fnp1->setAllNodeValue(0.0);

	this->fn->copy(this->seed);

	{
		DoubleProperty *tmp;
		Iterator<node> *itNodesU;
		Iterator<edge> *itEdges;
		node u, v;
		edge e;
		double num, denum, b, u0;
		bool continueProcess = true;

		if(pluginProgress)
			pluginProgress->setComment("Processing");

		for(unsigned int i = 0; i < iter_max; ++i) {
			itEdges = graph->getEdges();
			while(itEdges->hasNext()) {
				e = itEdges->next();
				u = graph->source(e);
				v = graph->target(e);
				this->beta->setEdgeValue(e, this->w->getEdgeValue(e) / (fabs(fn->getNodeValue(u) - fn->getNodeValue(v)) + 1));
			}
			delete itEdges;

			itNodesU = graph->getNodes();
			while(itNodesU->hasNext()) {
				u = itNodesU->next();

				itEdges = graph->getInOutEdges(u);
				num = 0; denum = 0;
				while(itEdges->hasNext()) {
					e = itEdges->next();
					v = graph->opposite(e, u);

					b = beta->getEdgeValue(e);

					num += b * fn->getNodeValue(v);
					denum += b;
				}
				delete itEdges;

				u0 = 1 - 2 * this->f0->getNodeValue(u);

				this->fnp1->setNodeValue(u,
						std::max( std::min( (num - lambda * u0) / denum, 1.0 ), 0.0 )
						);
			}
			delete itNodesU;

			tmp = fn;
			fn = fnp1;
			fnp1 = tmp;

			if(pluginProgress) {
				if(pluginProgress->state() != TLP_CONTINUE)
					continueProcess = false;

				if((i + 1) % 10 == 0) {
					pluginProgress->progress(i+1, iter_max);
					std::cerr << "Iteration " << (i+1) << "/" << iter_max << std::endl;
				}

				if((this->export_interval > 0) && ((i + 1) % this->export_interval) == 0) {
					std::cerr << "Iteration " << (i+1) << std::endl;
					try {
						this->result->copy(this->fn);
						fnToSelection();
						exportIntermediateResult(i+1);
					} catch (export_exception &ex) {
						std::cerr << "Export failed: " << ex.what() << std::endl;
					}
				}
			}

			if(!continueProcess)
				break;
		}

		this->result->copy(this->fn);

		fnToSelection();

		graph->delLocalProperty(this->beta_name);
		graph->delLocalProperty(this->fn_name);
		graph->delLocalProperty(this->fnp1_name);
	}

	if(pluginProgress && iter_max > 0) {
		pluginProgress->progress(iter_max, iter_max);
		pluginProgress->setComment("Computing selection");
	}

	if( this->export_interval > 0 )
	{
		try {
			exportIntermediateResult(iter_max);
		} catch (export_exception &ex) {
			std::cerr << "Export failed: " << ex.what() << std::endl;
		}
	}

	return true;
}
//=======================================================================

std::pair<double, double> Reg_ROF::computeFnMinMax() {
	double fn_min = DBL_MAX, fn_max = -DBL_MAX;

	Iterator<node> *itN = graph->getNodes();

	while (itN->hasNext()) {
		node itn=itN->next();

		double tmp = this->fn->getNodeValue(itn);

		if (tmp > fn_max) fn_max = tmp;

		if (tmp < fn_min) fn_min = tmp;
	}
	delete itN;

	return std::pair<double, double>(fn_min, fn_max);
}

void Reg_ROF::fnToSelection() 
{
	Iterator<node> *itNodesU;
	node u;
	const std::pair<double, double> fn_min_max = computeFnMinMax();
	const double threshold = (fn_min_max.first + fn_min_max.second) / 2.0;

	std::cout << "Threshold: " << threshold << std::endl;

	itNodesU = graph->getNodes();
	while(itNodesU->hasNext()) {
		u = itNodesU->next();
		if(fn->getNodeValue(u) > threshold) {
			this->segmentation_result->setNodeValue(u, 1);
		} else {
			this->segmentation_result->setNodeValue(u, 0);
		}
	}
	delete itNodesU;
}

void Reg_ROF::exportIntermediateResult(const int i) {
	std::ostringstream directory_name;
	directory_name << this->export_directory << "/" << std::setfill('0') << std::setw(6) << i << ".tlp";

	if(!tlp::saveGraph(graph, directory_name.str()))
	{
		throw export_exception(directory_name.str() + " cannot be written");
	}
}
