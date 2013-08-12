#include "Reg-CV.h"
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
			HTML_HELP_DEF( "type", "PropertyInterface" ) \
			HTML_HELP_BODY() \
			"Specify which property holds the data associated to each nodes. Can be a Double[Vector]Property." \
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

		// 5 Lambda 1
		HTML_HELP_OPEN() \
			HTML_HELP_DEF( "type", "Double" ) \
			HTML_HELP_BODY() \
			"Lambda 1 parameter." \
			HTML_HELP_CLOSE(),

		// 6 Lambda 2
		HTML_HELP_OPEN() \
			HTML_HELP_DEF( "type", "Double" ) \
			HTML_HELP_BODY() \
			"Lambda 2 parameter." \
			HTML_HELP_CLOSE(),

		// 7 Export interval
		HTML_HELP_OPEN() \
			HTML_HELP_DEF( "type", "Unsigned int" ) \
			HTML_HELP_BODY() \
			"Specify at which interval the processed graph must be exported (0 to disable)." \
			HTML_HELP_CLOSE(),

		// 8 Export directory
		HTML_HELP_OPEN() \
			HTML_HELP_DEF( "type", "Directory pathname" ) \
			HTML_HELP_BODY() \
			"This parameter is used to specify where the processed graph must be exported." \
			HTML_HELP_CLOSE(),

		// 9 Segmentation result
		HTML_HELP_OPEN() \
			HTML_HELP_DEF( "type", "BooleanProperty" ) \
			HTML_HELP_BODY() \
			"Indicates in which BooleanProperty the segmentation result must be saved." \
			HTML_HELP_CLOSE()
	};
}

//======================================================
Reg_CV::Reg_CV(const tlp::PluginContext *context):
	Algorithm(context),
	seed(NULL), result(NULL), f0(NULL), w(NULL), iter_max(0), lambda1(0.0), lambda2(0.0), export_interval(0), export_directory(),
	fn(NULL), fnp1(NULL), beta(NULL), f0_size(0), f0_scalar(NULL), f0_vector(NULL), segmentation_result(NULL),
	fn_name(), fnp1_name(), beta_name()
{
	addInParameter< DoubleProperty >      ("seed",                  paramHelp[0], "viewMetric");
	addInParameter< DoubleProperty >      ("result",                paramHelp[1], "viewMetric");
	addInParameter< BooleanProperty >     ("segmentation result",   paramHelp[9], "viewSelection");
	addInParameter< PropertyInterface* >  ("data",                  paramHelp[2], "");
	addInParameter< DoubleProperty >      ("similarity measure",    paramHelp[3], "");
	addInParameter< unsigned int >        ("number of iterations",  paramHelp[4], "100");
	addInParameter< double >              ("lambda1",               paramHelp[5], "1");
	addInParameter< double >              ("lambda2",               paramHelp[6], "1");
	addInParameter< unsigned int >        ("export interval",       paramHelp[7], "0");
	addInParameter< string >              ("dir::export directory", paramHelp[8], "");
}

bool Reg_CV::check(std::string &err)
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

		CHECK_PROP_PROVIDED("lambda1", this->lambda1);

		CHECK_PROP_PROVIDED("lambda2", this->lambda2);

		CHECK_PROP_PROVIDED("export interval", this->export_interval);

		if(this->iter_max < 0) {
			std::ostringstream m;
			m << "Invalid number of iterations: " << this->iter_max;
			throw std::runtime_error(m.str());
		}

		if (dynamic_cast< DoubleProperty* >(this->f0)) {
			this->f0_size = 1;
			this->f0_scalar = dynamic_cast< DoubleProperty* >(this->f0);
		} else if (dynamic_cast< DoubleVectorProperty* >(this->f0)) {
			this->f0_vector = dynamic_cast< DoubleVectorProperty* >(this->f0);
			this->f0_size = this->f0_vector->getNodeValue(graph->getOneNode()).size();
		} else {
			throw std::runtime_error("\"data\" must be a Double[Vector]Property.");
		}

		// No need for an export directory if export is disabled
		if(this->export_interval > 0) {
			CHECK_PROP_PROVIDED("dir::export directory", this->export_directory);

			// Checking if we can write in the export directory
			QString qstring_export_directory = QString::fromStdString(this->export_directory);
			QFileInfo info_export_directory(qstring_export_directory);
			QDir qdir_export_directory(qstring_export_directory);

			if(info_export_directory.exists()) {
				if(info_export_directory.isDir()) {
					if(qdir_export_directory.entryInfoList(QDir::NoDotAndDotDot | QDir::AllEntries).count() != 0) {
						throw std::runtime_error("Export directory (" + qdir_export_directory.absolutePath().toStdString() + ") is not empty.");
					}
				} else {
					throw std::runtime_error("Export directory (" + qdir_export_directory.absolutePath().toStdString() + ") already exists but is not a directory.");
				}
			} else {
				if(!qdir_export_directory.mkpath(".")) {
					throw std::runtime_error("Export directory (" + qdir_export_directory.absolutePath().toStdString() + ") cannot be created.");
				}
			}
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

		std::cout << "Processing graph " << graph->getName() << std::endl
		          << "Number of iterations: " << this->iter_max << std::endl
		          << "Length of the data: " << this->f0_size << std::endl
		          << "Lambda1: " << this->lambda1 << std::endl
		          << "Lambda2: " << this->lambda2 << std::endl
		          << "Export interval: " << this->export_interval << std::endl
		          << "Export directory: " << this->export_directory << std::endl;
	} catch (std::runtime_error &ex) {
		err.assign(ex.what());
		return false;
	}

	return true;
}

//======================================================
bool Reg_CV::run() {
	this->fn->setAllNodeValue(0.0);
	this->fnp1->setAllNodeValue(0.0);

	this->fn->copy(this->seed);

	this->in_out_means.first = std::vector< double >(this->f0_size, 0.0);
	this->in_out_means.second = std::vector< double >(this->f0_size, 0.0);

	{
		DoubleProperty *tmp;
		Iterator<node> *itNodesU;
		Iterator<edge> *itEdges;
		node u, v;
		edge e;
		double num, denum, b, cv_criteria, cv_criteria_cumulated;
		std::vector< double > u0(this->f0_size);
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

			computeMeanValues();

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

				if(this->f0_scalar) {
					u0[0] = this->f0_scalar->getNodeValue(u);
				} else {
					u0 = this->f0_vector->getNodeValue(u);
				}
				cv_criteria_cumulated = 0;
				for(unsigned int j = 0; j < f0_size; ++j) {
					cv_criteria = in_out_means.first[j] - u0[j];
					cv_criteria_cumulated += lambda1 * cv_criteria * cv_criteria;
					cv_criteria = in_out_means.second[j] - u0[j];
					cv_criteria_cumulated -= lambda2 * cv_criteria * cv_criteria;
				}

				this->fnp1->setNodeValue(u,
						std::max(
							std::min(
								(num - cv_criteria_cumulated / f0_size) / denum,
								1.0
							   ),
							0.0
						   )
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

std::pair<double, double> Reg_CV::computeFnMinMax() {
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

void Reg_CV::fnToSelection() 
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

void Reg_CV::exportIntermediateResult(const int i) {
	std::ostringstream directory_name;
	directory_name << this->export_directory << "/" << std::setfill('0') << std::setw(6) << i << ".tlp";

	if(!tlp::saveGraph(graph, directory_name.str()))
	{
		throw export_exception(directory_name.str() + " cannot be written");
	}
}

void Reg_CV::computeMeanValues()
{
	int n1 = 0, n2 = 0;
	node n;
	unsigned int i;
	std::vector< double > f0_value(1);

	std::fill(this->in_out_means.first.begin(), this->in_out_means.first.end(), 0.0);
	std::fill(this->in_out_means.second.begin(), this->in_out_means.second.end(), 0.0);

	Iterator<node> *itNodes = graph->getNodes();
	while(itNodes->hasNext()) {
		n = itNodes->next();

		if(this->f0_scalar) {
			f0_value[0] = this->f0_scalar->getNodeValue(n);
		} else {
			f0_value = this->f0_vector->getNodeValue(n);
		}
		if(this->fn->getNodeValue(n) >= 0.5) {
			for(i = 0; i < f0_size; ++i) {
				this->in_out_means.first[i] += f0_value[i];
			}
			++n1;
		} else {
			for(i = 0; i < f0_size; ++i) {
				this->in_out_means.second[i] += f0_value[i];
			}
			++n2;
		}
	}
	delete itNodes;

	for(i = 0; i < f0_size; ++i) {
		if(n1 == 0)
			this->in_out_means.first[i] = 1.0;
		else
			this->in_out_means.first[i] /= n1;

		if(n2 == 0)
			this->in_out_means.second[i] = 0.0;
		else
			this->in_out_means.second[i] /= n2;
	}
}
