#include "Cv_Ta.h"
#include <sstream>
#include <iostream>
#include <iomanip>

#include <QFileInfo>
#include <QDir>

ALGORITHMPLUGIN(Cv_Ta, "Cv_Ta", "Cyrille FAUCHEUX","20/11/2011","Alpha","1.1");

using namespace std;
using namespace tlp;

inline double max(double v1, double v2) {
	return (v1 > v2 ? v1 : v2);
}

inline double min(double v1, double v2) {
	return (v1 > v2 ? v2 : v1);
}

template <class T>
std::string to_string(T t, std::ios_base & (*f)(std::ios_base&))
{
	std::ostringstream oss;
	oss << f << t;
	return oss.str();
}

namespace {
	const char * paramHelp[] = {
		HTML_HELP_OPEN() \
			HTML_HELP_DEF( "type", "BooleanProperty" ) \
			HTML_HELP_DEF( "default", "mask" ) \
			HTML_HELP_BODY() \
			"The boolean property that defines the initial selection." \
			HTML_HELP_CLOSE(),

		HTML_HELP_OPEN() \
			HTML_HELP_DEF( "type", "Unsigned int" ) \
			HTML_HELP_DEF( "default", "1000" ) \
			HTML_HELP_BODY() \
			"Defines the number of iterations the algorithm should do." \
			HTML_HELP_CLOSE(),

		HTML_HELP_OPEN() \
			HTML_HELP_DEF( "type", "Double" ) \
			HTML_HELP_DEF( "default", "0.25" ) \
			HTML_HELP_BODY() \
			"Lambda 1 parameter from the original formula." \
			HTML_HELP_CLOSE(),

		HTML_HELP_OPEN() \
			HTML_HELP_DEF( "type", "DoubleVectorProperty" ) \
			HTML_HELP_DEF( "default", "data" ) \
			HTML_HELP_BODY() \
			"Specify which property holds the data associated to each pixels/nodes." \
			HTML_HELP_CLOSE(),

		HTML_HELP_OPEN() \
			HTML_HELP_DEF( "type", "DoubleProperty" ) \
			HTML_HELP_DEF( "default", "weight" ) \
			HTML_HELP_BODY() \
			"Specify which property holds the weight associated to each edge." \
			HTML_HELP_CLOSE(),

		HTML_HELP_OPEN() \
			HTML_HELP_DEF( "type", "BooleanProperty" ) \
			HTML_HELP_DEF( "default", "roi" ) \
			HTML_HELP_BODY() \
			"The boolean property that defines the region of interest." \
			HTML_HELP_CLOSE(),

		HTML_HELP_OPEN() \
			HTML_HELP_DEF( "type", "Double" ) \
			HTML_HELP_DEF( "default", "0.25" ) \
			HTML_HELP_BODY() \
			"Lambda 2 parameter from the original formula." \
			HTML_HELP_CLOSE(),

		HTML_HELP_OPEN() \
			HTML_HELP_DEF( "type", "Unsigned int" ) \
			HTML_HELP_DEF( "default", "50" ) \
			HTML_HELP_BODY() \
			"Specify at which interval the processed image must be exported. (0 to disable intermediate export)" \
			HTML_HELP_CLOSE(),

		HTML_HELP_OPEN() \
			HTML_HELP_DEF( "type", "Directory pathname" ) \
			HTML_HELP_DEF( "default", "/tmp" ) \
			HTML_HELP_BODY() \
			"This parameter is used to specify where the processed images must be exported." \
			HTML_HELP_CLOSE(),
	};
}

//======================================================
Cv_Ta::Cv_Ta(const tlp::AlgorithmContext &context):Algorithm(context), fn(NULL), fnp1(NULL), beta(NULL), w(NULL), f0(NULL), mask(NULL), roi(NULL) {
	addDependency<Algorithm>("Export image 3D mask", "1.0");
	addParameter<DoubleVectorProperty>("Data", paramHelp[3]);
	addParameter<BooleanProperty>("Mask", paramHelp[0]);
	addParameter<unsigned int>("Number of iterations", paramHelp[1], "1000");
	addParameter<double>("Lambda1", paramHelp[2], "0.25");
	addParameter<double>("Lambda2", paramHelp[6], "0.25");
	addParameter<DoubleProperty>("Weight", paramHelp[4]);
	addParameter<BooleanProperty>("Roi", paramHelp[5]);
	addParameter<unsigned int>("Export interval", paramHelp[7], "50");
  addParameter<string>("dir::Export directory", paramHelp[8]);
}

template<typename T>
void PrintVector(const std::vector<T>& t){
  for(typename std::vector<T>::size_type i=0; i<t.size(); ++i){
    std::cout << t[i] << " - ";
  }
  std::cout << std::endl << std::endl;
}

#define CHECK_PROP_PROVIDED(PROP, STOR) \
  if(!dataSet->get(#PROP, STOR)) \
    throw std::runtime_error(std::string("No \"") + #PROP + "\" provided")

bool Cv_Ta::check(std::string &err)
{
  try {
    if(dataSet == NULL)
      throw std::runtime_error("No dataset provided.");

    CHECK_PROP_PROVIDED(Data, this->f0);

    CHECK_PROP_PROVIDED(Mask, this->mask);

    CHECK_PROP_PROVIDED(Number of iterations, this->iter_max);

    CHECK_PROP_PROVIDED(Lambda1, this->lambda1);

    CHECK_PROP_PROVIDED(Lambda2, this->lambda2);

    CHECK_PROP_PROVIDED(Weight, this->w);

    CHECK_PROP_PROVIDED(Roi, this->roi);

    CHECK_PROP_PROVIDED(Export interval, this->export_interval);

    CHECK_PROP_PROVIDED(dir::Export directory, this->export_directory);

    if(this->iter_max <= 0) {
      std::ostringstream err;
      err << "Invalid number of iterations: " << this->iter_max;
      throw std::runtime_error(err.str());
    }

    { // Checking if we can write in the export directory
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

    std::cout << "Number of iterations: " << this->iter_max << std::endl;
    std::cout << "Lambda1: " << this->lambda1 << std::endl;
    std::cout << "Lambda2: " << this->lambda2 << std::endl;
    std::cout << "Export interval: " << this->export_interval << std::endl;
    std::cout << "Export directory: " << this->export_directory << std::endl;
  } catch (std::runtime_error &ex) {
    err.assign(ex.what());
    return false;
  }

  return true;
}

//======================================================
bool Cv_Ta::run() {
	this->fn = graph->getLocalProperty<DoubleProperty>("fn");
	this->fnp1 = graph->getLocalProperty<DoubleProperty>("fnp1");
	this->fn->setAllNodeValue(0.0);

	{ // Initializing Fn from the mask
		Iterator<node> *itNodes = graph->getNodes();
		node n;

		if(itNodes->hasNext()) {
			n = itNodes->next();
			this->f0_size = this->f0->getNodeValue(n).size();
		}
		delete itNodes;


		this->in_out_means.first = std::vector< double >(this->f0_size, 0.0);
		this->in_out_means.second = std::vector< double >(this->f0_size, 0.0);

		itNodes = graph->getNodes();
		while(itNodes->hasNext()) {
			n = itNodes->next();
			if(this->roi->getNodeValue(n))
				this->fn->setNodeValue(n, this->mask->getNodeValue(n) ? 1.0 : 0.0);
		}
		delete itNodes;

		computeMeanValues();
	}

  std::cout << "Fn initialized" << std::endl;

	{
		DoubleProperty *tmp;
		Iterator<node> *itNodesU;
		Iterator<edge> *itEdges;
		node u, v;
		edge e;
		double num, denum, b, cv_criteria, cv_criteria_cumulated;
		std::vector< double > u0;
		bool continueProcess = true;

		if(pluginProgress)
			pluginProgress->setComment("Processing");

		beta = graph->getLocalProperty<DoubleProperty>("beta");
		for(unsigned int i = 0; i < iter_max; ++i) {
			itEdges = graph->getEdges();
			while(itEdges->hasNext()) {
				e = itEdges->next();
				u = graph->source(e);
				v = graph->target(e);
				if(this->roi->getNodeValue(u) && this->roi->getNodeValue(v))
					beta->setEdgeValue(e, this->w->getEdgeValue(e) / (fabs(fn->getNodeValue(u) - fn->getNodeValue(v)) + 1));
			}
			delete itEdges;

			computeMeanValues();

			itNodesU = graph->getNodes();
			while(itNodesU->hasNext()) {
				u = itNodesU->next();

				if(!this->roi->getNodeValue(u))
					continue;

				itEdges = graph->getInOutEdges(u);
				num = 0; denum = 0;
				while(itEdges->hasNext()) {
					e = itEdges->next();
					b = beta->getEdgeValue(e);

					v = graph->opposite(e, u);

					if(this->roi->getNodeValue(v)) {
						num += b * fn->getNodeValue(v);
						denum += b;
					}
				}
				delete itEdges;

				u0 = f0->getNodeValue(u);
				cv_criteria_cumulated = 0;
				for(unsigned int j = 0; j < f0_size; ++j) {
					cv_criteria = in_out_means.first[j] - u0[j];
					cv_criteria_cumulated += lambda1 * cv_criteria * cv_criteria;
					cv_criteria = in_out_means.second[j] - u0[j];
					cv_criteria_cumulated -= lambda2 * cv_criteria * cv_criteria;
				}

				fnp1->setNodeValue(u,
						max(
							min(
								(num - cv_criteria_cumulated / f0_size) / denum,
								1
								),
							0
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

        if(i % 10 == 0)
          pluginProgress->progress(i, iter_max);

				if(/*pluginProgress->isPreviewMode() &&*/ ((i+1) % this->export_interval == 0)) {
					std::cerr << "Iteration " << (i+1) << std::endl;
					fnToSelection();
          try {
  					exportSelection(i+1);
          } catch (export_exception &ex) {
            std::cerr << "Export failed: " << ex.what() << std::endl;
          }
				}
			}

      if(!continueProcess)
        break;
		}
	}

	if(pluginProgress) {
		pluginProgress->progress(iter_max, iter_max);
		pluginProgress->setComment("Computing selection");
	}

	fnToSelection();
  try {
    exportSelection(iter_max);
  } catch (export_exception &ex) {
    std::cerr << "Export failed: " << ex.what() << std::endl;
  }

	graph->delLocalProperty("fn");
	graph->delLocalProperty("fnp1");
	graph->delLocalProperty("beta");

	std::cout << "FINISHED" << std::endl;

	return true;
}
//=======================================================================

void Cv_Ta::fnToSelection() 
{
	Iterator<node> *itNodesU;
	node u;
	BooleanProperty *selection = graph->getLocalProperty<BooleanProperty>("viewSelection");
	DoubleProperty *fn = graph->getLocalProperty<DoubleProperty>("fn");
	double threshold = (fn->getNodeMax(graph)-fn->getNodeMin(graph)) / 2.0;

	std::cout << "Threshold: " << threshold << std::endl;

	itNodesU = graph->getNodes();
	while(itNodesU->hasNext()) {
		u = itNodesU->next();
		if(fn->getNodeValue(u) > threshold) {
			selection->setNodeValue(u, 1);
		} else {
			selection->setNodeValue(u, 0);
		}
	}
	delete itNodesU;
}

void Cv_Ta::exportSelection(const int i) {
  int depth = -1;

  if(!graph->getAttribute<int>("depth", depth)) {
    throw export_exception("The graph does not seems to be an image");
  }

  std::ostringstream directory_name;
  directory_name << std::setfill('0') << std::setw(6) << i;

  std::string directory = this->export_directory;
  std::string pattern = "%1.bmp";

  if(depth == 1) {
    pattern = directory_name.str() + "_" +  pattern;
  } else {
    directory += "/" + directory_name.str();

    QDir qdir_directory(directory.c_str());

    if(!qdir_directory.mkpath(".")) {
      throw export_exception("The directory " + QDir::toNativeSeparators(QString(directory.c_str())).toStdString() + " cannot be created");
    }
  }
  
	BooleanProperty *selection = graph->getLocalProperty<BooleanProperty>("viewSelection");
	dataSet->set<PropertyInterface*>("Mask property", selection);
	dataSet->set<string>("dir::Output folder", directory);
	dataSet->set<string>("Pattern", pattern);

	string message;

	if(!graph->applyAlgorithm("Export image 3D mask", message, dataSet, NULL))
    throw export_exception("Unable to export the image (" + message + ")");
}

void Cv_Ta::computeMeanValues()
{
	this->in_out_means.first = std::vector< double >(this->f0_size, 0.0);
	this->in_out_means.second = std::vector< double >(this->f0_size, 0.0);
	int n1 = 0, n2 = 0;
	node n;
	unsigned int i;
	std::vector< double > f0_value;

	Iterator<node> *itNodes = graph->getNodes();
	while(itNodes->hasNext()) {
		n = itNodes->next();
		
		if(!this->roi->getNodeValue(n))
			continue;

		f0_value = this->f0->getNodeValue(n);
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
		this->in_out_means.first[i] /= n1 + 1;
		this->in_out_means.second[i] /= n2 + 1;
	}
}


