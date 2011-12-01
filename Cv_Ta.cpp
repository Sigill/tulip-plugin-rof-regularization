#include "Cv_Ta.h"
#include <sstream>
#include <iostream>

ALGORITHMPLUGIN(Cv_Ta, "Cv_Ta", "Cyrille FAUCHEUX","20/11/2011","Alpha","1.0");

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
			HTML_HELP_DEF( "default", "4.0" ) \
			HTML_HELP_BODY() \
			"Lambda parameter from the original formula." \
			HTML_HELP_CLOSE(),

		HTML_HELP_OPEN() \
			HTML_HELP_DEF( "type", "DoubleProperty" ) \
			HTML_HELP_BODY() \
			"Specify which property holds the data associated to each pixels/nodes." \
			HTML_HELP_CLOSE(),
	};
}

//======================================================
Cv_Ta::Cv_Ta(const tlp::AlgorithmContext &context):Algorithm(context) {
	addDependency<Algorithm>("Export image 3D mask", "1.0");
	addParameter<DoubleProperty>("Data", paramHelp[3]);
	addParameter<BooleanProperty>("Mask", paramHelp[0]);
	addParameter<unsigned int>("Number of iterations", paramHelp[1], "1000");
	addParameter<double>("Lambda", paramHelp[2], "4.0");
}

//======================================================
bool Cv_Ta::run() {
	this->iter_max = 1000;
	this->lambda = 4.0;
	this->fn = graph->getLocalProperty<DoubleProperty>("fn");
	this->fnp1 = graph->getLocalProperty<DoubleProperty>("fnp1");

	BooleanProperty *mask = NULL;

	if (dataSet != 0) {
		dataSet->get<DoubleProperty*>("Data", this->f0);
		dataSet->get<BooleanProperty*>("Mask", mask);
		dataSet->get("Number of iterations", iter_max);
		dataSet->get("Lambda", lambda);
	}

	{
		Iterator<node> *itNodes = graph->getNodes();
		node n;
		while(itNodes->hasNext()) {
			n = itNodes->next();
			fn->setNodeValue(n, !mask->getNodeValue(n) ? 1 : 0);
		}
		delete itNodes;

		computeMeanValues();
	}

	{
		DoubleProperty *tmp;
		Iterator<node> *itNodesU;
		Iterator<edge> *itEdges;
		node u, v;
		edge e;
		double num, denum, b, v1, v2, u0;
		bool continueProcess = true;

		if(pluginProgress)
			pluginProgress->setComment("Processing");

		beta = graph->getLocalProperty<DoubleProperty>("beta");
		for(unsigned int i = 0; i < iter_max && continueProcess; ++i) {
			itEdges = graph->getEdges();
			while(itEdges->hasNext()) {
				e = itEdges->next();
				beta->setEdgeValue(e, 1/(fabs(fn->getNodeValue(graph->source(e)) - fn->getNodeValue(graph->target(e))) + 1));
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
					b = beta->getEdgeValue(e);

					v = graph->opposite(e, u);

					num += b * fn->getNodeValue(v);
					denum += b;
				}
				delete itEdges;

				u0 = f0->getNodeValue(u);
				v2 = in_out_means.first - u0;
				v1 = in_out_means.second - u0;

				fnp1->setNodeValue(u,
						max(
							min(
								(num - lambda * (v1 * v1 - v2 * v2)) / denum,
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

				pluginProgress->progress(i, iter_max);
				if(/*pluginProgress->isPreviewMode() &&*/ (i % 20 == 0)) {
					fnToSelection();
					exportSelection(i);
					//std::cout << "Threshold: " << (fn->getNodeMax(graph) / 2.0) << std::endl;
				}
			}
		}
	}

	if(pluginProgress) {
		pluginProgress->progress(iter_max, iter_max);
		pluginProgress->setComment("Computing selection");
	}

	fnToSelection();
	exportSelection(iter_max);

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
	double threshold = fn->getNodeMax(graph) / 2.0;

	std::cout << "Threshold: " << threshold << std::endl;

	itNodesU = graph->getNodes();
	while(itNodesU->hasNext()) {
		u = itNodesU->next();
		if(fn->getNodeValue(u) < threshold) {
			selection->setNodeValue(u, 1);
		} else {
			selection->setNodeValue(u, 0);
		}
	}
	delete itNodesU;
}

void Cv_Ta::exportSelection(const int i) {
	BooleanProperty *selection = graph->getLocalProperty<BooleanProperty>("viewSelection");
	string message;
	dataSet->set<PropertyInterface*>("Mask property", selection);
	dataSet->set<string>("dir::Output folder", "/tmp");
	dataSet->set<string>("Pattern", to_string<int>(i, std::dec) + "_%1.bmp");
	if(!applyAlgorithm(graph, message, dataSet, "Export image 3D mask", pluginProgress))
	{
		std::cerr << "Error: Unable to load the mask" << std::endl;
		if(pluginProgress) pluginProgress->setError("Error: Unable to load the mask");
	}
}

void Cv_Ta::computeMeanValues()
{
	double c1 = 0, c2 = 0;
	int n1 = 0, n2 = 0;
	node n;

	Iterator<node> *itNodes = graph->getNodes();
	while(itNodes->hasNext()) {
		n = itNodes->next();
		if(this->fn->getNodeValue(n) <= 0.5) {
			c1 += this->f0->getNodeValue(n);
			++n1;
		} else {
			c2 += this->f0->getNodeValue(n);
			++n2;
		}
	}
	delete itNodes;

	this->in_out_means.first = c1 / (n1 + 1);
	this->in_out_means.second = c2 / (n2 + 1);
}


