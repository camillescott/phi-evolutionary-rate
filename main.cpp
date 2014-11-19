#include <stdio.h>
#include <stdlib.h>

#include <vector>
#include <map>
#include <math.h>
#include <time.h>
#include <iostream>
#include <thread>
#include <functional> // for ref
#include "globalConst.h"
#include "tHMM.h"
#include "tAgent.h"
#include "tGame.h"
#include "params/params.h"

#ifdef _WIN32
	#include <process.h>
#else
	#include <unistd.h>
#endif

#define randDouble ((double)rand()/(double)RAND_MAX)

using namespace std;
using namespace Params;

double perSiteMutationRate=0.005;
int update=0;
size_t repeats=1;
int maxAgent=100;
int totalGenerations=200;
const int mp=2;

void computeLOD(FILE *f,FILE *g, tAgent *agent,tGame *game);
const char* cstr(string s) { return s.c_str(); }

void threadedExecuteGame(int chunkBegin, int chunkEnd, const vector<tAgent*>& agent, tGame& game) {
	for (int i=chunkEnd-chunkBegin-1; i>=0; --i) {
		game.executeGame(agent[chunkBegin+i], 2, nullptr, true, -1, -1);
		agent[chunkBegin+i]->fitnesses.push_back(agent[chunkBegin+i]->fitness);
	}
}

int main(int argc, char *argv[]) {
	string experimentID;
	int replicateID=0;
	vector<tAgent*>agent;
	vector<tAgent*>nextGen;
	int who=0;
	size_t i, j;
	double maxFitness;
	tAgent *masterAgent=nullptr;
	tGame *game=nullptr;
	FILE *LODFile=nullptr;
	FILE *genomeFile=nullptr;
	bool showhelp;
	float evolvePhiLimit;
	int evolvePhiGenLimit;
	float evolveRLimit;
	int evolveRGenLimit;
	float evolveTopologyLimit;
	int evolveTopologyGenLimit;
	bool stopOnLimit;
	string filenameLOD, filenameGenome, filenameStartWith;
	int nthreads=1;
	vector<thread> threads;

	addp(TYPE::BOOL, &showhelp);
	addp(TYPE::STRING, &filenameLOD, "--LOD", "filename to save Line of Descent.");
	addp(TYPE::STRING, &experimentID, "--experiment", "unique identifier for this experiment, shared by all replicates.");
	addp(TYPE::INT, &replicateID, "--replicate", "unique number to identify this replicate in this experiment.");
	addp(TYPE::STRING, &filenameGenome, "--genome", "filename to save genomes of the LODFile.");
	addp(TYPE::FLOAT, &evolvePhiLimit, "-1.0", false, "--evolvePhiLimit", "phi threshold for brain selection during evolution before switching to task fitness. Define to enable.");
	addp(TYPE::INT, &totalGenerations, "200", false, "--generations", "number of generations to simulate (updates).");
	addp(TYPE::STRING, &filenameStartWith, "none", false, "--startwith", "specify a genome file used to seed the population.");
	addp(TYPE::INT, &evolvePhiGenLimit, "-1", false, "--evolvePhiGenLimit", "instead of evolving to a value of phi, number of generations.");
	addp(TYPE::FLOAT, &evolveRLimit, "-1.0", false, "--evolveRLimit", "R threshold for brain selection during evolution before switching to task fitness. Define to enable.");
	addp(TYPE::INT, &evolveRGenLimit, "-1", false, "--evolveRGenLimit", "instead of evolving to a vlue of R, number of generations.");
	addp(TYPE::FLOAT, &evolveTopologyLimit, "-1.0", false, "--evolveTopologyLimit", "Topology (avg path length) threshold for brain selection during evolution before switching to task fitness. Define to enable.");
	addp(TYPE::INT, &evolveTopologyGenLimit, "-1", false, "--evolveTopologyGenLimit", "instead of evolving to a value of Topology, number of generations.");
	addp(TYPE::INT, &nthreads, "1", false, "--nthreads", "number of worker threads");
	addp(TYPE::BOOL, &stopOnLimit, "false", false, "--stopOnLimit", "if a limit is specified, then the simulation will stop at the limit.");
	argparse(argv);
	if (showhelp) {
		cout << argdetails() << endl;
		exit(0);
	}

    srand(getpid());
    masterAgent=new tAgent();

	LODFile=fopen(cstr(filenameLOD),"w+t");
	genomeFile=fopen(cstr(filenameGenome),"w+t");	
	srand(getpid());
	agent.resize(maxAgent);
	masterAgent=new tAgent;
	vector<vector<int> > data;
	game=new tGame;
	masterAgent->setupRandomAgent(5000);
	
	masterAgent->setupPhenotype();

	if (filenameStartWith != "none") {
		masterAgent->loadAgent(filenameStartWith.c_str());
		for(i=0;i<agent.size();i++){
			agent[i]=new tAgent;
			agent[i]->inherit(masterAgent,0.0,0); // small mutation rate to preserve genome
		}
	} else {
		for(i=0;i<agent.size();i++){
			agent[i]=new tAgent;
			agent[i]->inherit(masterAgent,0.5,0); // large mutation rate to spread the genotypic population
		}
	}
	nextGen.resize(agent.size());
	masterAgent->nrPointingAtMe--;
	cout<<"setup complete"<<endl;
	printf("%s	%s	%s	%s	%s	%s	%s\n", "update","(double)maxFitness","maxPhi","r", "maxTopology", "agent[who]->correct","agent[who]->incorrect");

	while(update<totalGenerations){
		for(i=0;i<agent.size();i++){
			agent[i]->fitness=0.0;
			agent[i]->phitness=0.0;
			agent[i]->fitnesses.clear();
		}

		threads.clear();
		int chunksize=agent.size()/nthreads;
        for (int i=0; i<nthreads; ++i) {
		    threads.push_back(thread(threadedExecuteGame, i*chunksize, (i*chunksize)+chunksize,
                ref(agent), ref(*game)));
        }
		//threadedExecuteGame(chunksize, agent.size(), ref(agent), ref(*game));
		for (thread& t : threads) t.join(); // sync threads

		maxFitness=0.0;
		double maxPhi=0.0;
		double maxR=0.0;
		double maxTopology=0;
		bool evolvePhiPhase=false;
		bool evolvingR = false;
		bool evolvingTopology = false;
		
		for(i=0;i<agent.size();i++){
			agent[i]->fitness=agent[i]->fitnesses[0];
			if(agent[i]->fitness>maxFitness) {
				maxFitness=agent[i]->fitness;
				who=i;
			}
			if(pow(1.1,mp*agent[i]->phi)>maxPhi)
				maxPhi=pow(1.1,mp*agent[i]->phi);
			if (pow(1.1,mp*agent[i]->R)>maxR)
				maxR=pow(1.1,mp*agent[i]->R);
			if (pow(1.1,mp*agent[i]->Topology)>maxTopology)
				maxTopology=pow(1.1,mp*agent[i]->Topology);
		}
		if (evolvePhiLimit > 0.0f) evolvePhiPhase=true;
		if ((evolvePhiGenLimit != -1) && (update < evolvePhiGenLimit)) evolvePhiPhase = true;
		if (evolveRLimit > 0.0f) evolvingR = true;
		if ((evolveRGenLimit > 0) && (update < evolveRGenLimit)) evolvingR = true;
		if (evolveTopologyLimit > 0) evolvingTopology = true;
		if ((evolveTopologyGenLimit > 0) && (update < evolveTopologyGenLimit)) evolvingTopology = true;
		  printf("%i	%f	%f	%f	%f	%i	%i\n", update, (double)maxFitness, (log10(maxPhi)/log10(1.1))/mp, (log10(maxR)/log10(1.1))/mp, (log10(maxTopology)/log10(1.1))/mp, agent[who]->correct, agent[who]->incorrect);

		  threads.clear();
		if (evolvePhiPhase || evolvingR || evolvingTopology) {
			int j=0;
			for(i=0;i<agent.size();i++) {
				tAgent *d;
				d=new tAgent;
				float selectedValue=0.0f;
				float maxValue=0.0f;
				if (evolvingR) maxValue = maxR;
				else if (evolvingTopology) maxValue = maxTopology;
				else maxValue = maxPhi;
				if(maxValue<=0.0){
					j=rand()%(int)agent.size();
				} else {
					do{
						j=rand()%(int)agent.size();
						if (evolvingR) {
							selectedValue = agent[j]->R;
							maxValue = maxR;
						} else if (evolvePhiPhase) {
							selectedValue = agent[j]->phi;
							maxValue = maxPhi;
						} else {
							selectedValue = (float)agent[j]->Topology;
							maxValue = (float)maxTopology;
						}
					} while((j==(i))||(randDouble>(pow(1.1,mp*selectedValue)/maxValue)));
				}
				d->inherit(agent[j],perSiteMutationRate,update);
				nextGen[i]=d;
			}
			if ((log10(maxPhi)/log10(1.1))/mp > evolvePhiLimit) {
				if (evolvePhiLimit > 0) {
					evolvePhiLimit = -1.0f;
					if (stopOnLimit) {
						break;
					}
				}
			}
			if ((log10(maxR)/log10(1.1))/mp > evolveRLimit) {
				if (evolveRLimit > 0) {
					evolveRLimit = -1.0f;
					if (stopOnLimit) {
						break;
					}
				}
			}
			if ((log10(maxTopology)/log10(1.1))/mp > evolveTopologyLimit) {
				if (evolveTopologyLimit > 0) {
					evolveTopologyLimit = -1.0f;
					if (stopOnLimit) {
						break;
					}
				}
			}
		} else {
			for(i=0;i<agent.size();i++)
			{
				tAgent *d;
				d=new tAgent;
				if(maxFitness<=0.0){
					j=rand()%(int)agent.size();
				} else
				do{ j=rand()%(int)agent.size(); } while((j==i)||(randDouble>(agent[j]->fitness/maxFitness)));
				d->inherit(agent[j],perSiteMutationRate,update);
				nextGen[i]=d;
			}
		}
		for(i=0;i<agent.size();i++){
			agent[i]->retire();
			agent[i]->nrPointingAtMe--;
			if(agent[i]->nrPointingAtMe==0)
				delete agent[i];
			agent[i]=nextGen[i];
		}
		agent=nextGen;
		update++;
	}
	
	agent[0]->ancestor->saveLOD(LODFile,genomeFile, experimentID, replicateID, -1); // -1 to tell saveLOD to make header for csv
	if (stopOnLimit) {
		float maxPhi=0.0;
		tAgent* bestAgent=nullptr;
		for (tAgent* a : agent) {
			if (a->phi > maxPhi) {
				maxPhi = a->phi;
				bestAgent = a;
			}
		}
		if (bestAgent) {
			bestAgent->saveGenome(genomeFile);
		}
	} else {
		agent[0]->ancestor->ancestor->saveGenome(genomeFile);
	}
//	agent[0]->ancestor->saveToDot(argv[3]);
	agent.clear();
	nextGen.clear();
	delete masterAgent;
	delete game;
	return 0;
}

void computeLOD(FILE *f,FILE *g, tAgent *agent,tGame *game){
	/*vector<vector<int> > table;
	double R,oldR;
	if(agent->ancestor!=NULL)
		computeLOD(f,g,agent->ancestor,game);
	agent->setupPhenotype();
	table=game->executeGame(agent, 2, NULL,false,-1,-1);
	R=game->computeR(table,0);
	oldR=game->computeOldR(table);
	fprintf(f,"%i	%i	%i	%f	%f",agent->ID,agent->correct,agent->incorrect,agent->extra);
	fprintf(f,"\n");
	*/
}

