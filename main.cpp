#include <stdio.h>
#include <stdlib.h>

#include <vector>
#include <map>
#include <math.h>
#include <time.h>
#include <iostream>
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
int totalGenerations=400;

void computeLOD(FILE *f,FILE *g, tAgent *agent,tGame *game);
const char* cstr(string s) { return s.c_str(); }

int main(int argc, char *argv[]) {
	string experimentID;
	int replicateID=0;
	vector<tAgent*>agent;
	vector<tAgent*>nextGen;
	tAgent *masterAgent;
	int who=0;
	size_t i, j;
	tGame *game;
	double maxFitness;
	FILE *LODFile;
	FILE *genomeFile;
	bool showhelp;
	float preselectPhiLimit;
	float evolvePhiLimit;
	string filenameLOD, filenameGenome;

	addp(TYPE::BOOL, &showhelp);
	addp(TYPE::STRING, &filenameLOD, "--LOD", "filename to save Line of Descent.");
	addp(TYPE::STRING, &experimentID, "--experiment", "unique identifier for this experiment, shared by all replicates.");
	addp(TYPE::INT, &replicateID, "--replicate", "unique number to identify this replicate in this experiment.");
	addp(TYPE::STRING, &filenameGenome, "--genome", "filename to save genomes of the LODFile.");
	addp(TYPE::FLOAT, &preselectPhiLimit, "-1.0", false, "--preselectPhi", "phi threshold for brain selection before evolution. Define to enable.");
	addp(TYPE::FLOAT, &evolvePhiLimit, "-1.0", false, "--evolvePhi", "phi threshold for brain selection during evolution before switching to task fitness. Define to enable.");
	argparse(argv);
	if (showhelp) {
		cout << argdetails() << endl;
		exit(0);
	}

    srand(getpid());
    masterAgent=new tAgent;
	/*LODFile=fopen(cstr(filenameLOD),"r+t");
	if(LODFile!=NULL){
		fclose(LODFile);
		exit(0);
	}*/

	LODFile=fopen(cstr(filenameLOD),"w+t");
	genomeFile=fopen(cstr(filenameGenome),"w+t");	
	srand(getpid());
	agent.resize(maxAgent);
    masterAgent=new tAgent;
    vector<vector<int> > data;
	game=new tGame;
	masterAgent=new tAgent;
	masterAgent->setupRandomAgent(5000);
	
	masterAgent->setupPhenotype();
	//masterAgent->loadAgent("startGenome.txt");

	for(i=0;i<agent.size();i++){
		agent[i]=new tAgent;
		agent[i]->inherit(masterAgent,0.5,0);
	}
	nextGen.resize(agent.size());
	masterAgent->nrPointingAtMe--;
	cout<<"setup complete"<<endl;
	printf("%s	%s	%s	%s	%s	%s\n", "update","(double)maxFitness","maxPhi","agent[who]->phi","agent[who]->correct","agent[who]->incorrect");
	while(update<totalGenerations){
		//*
		for(i=0;i<agent.size();i++){
			agent[i]->fitness=1.0;
			agent[i]->phitness=0.0;
			agent[i]->fitnesses.clear();
		}
		for(i=0;i<agent.size();i++){
			for(j=0;j<repeats;j++){
				game->executeGame(agent[i], 2, NULL, true, -1, -1);
				//agent[i]->fitness*=agent[i]->fitness;
				//agent[i]->fitness+=1.0;
				agent[i]->fitnesses.push_back(agent[i]->fitness);
				/*agent[i]->fitness*=(1.0+agent[i]->phi);*/
				//agent[i]->phitness+=(agent[i]->phi / (double)repeats);
				//agent[i]->fitnesses.push_back(agent[i]->phi);
			}
		}
		maxFitness=0.0;
		double maxPhi=0.0;
		bool evolvePhiPhase=false;
		if (evolvePhiLimit > 0.0f) evolvePhiPhase=true;
		
		for(i=0;i<agent.size();i++){
			agent[i]->fitness=agent[i]->fitnesses[0];
			//agent[i]->fitness=agent[i]->phitness;
			if(agent[i]->fitness>maxFitness)
				maxFitness=agent[i]->fitness;
			if(pow(agent[i]->phi,2)>maxPhi)
				maxPhi=pow(agent[i]->phi,2);
		}
        //if((update&31)==0)
		  printf("%i	%f	%f	%f	%i	%i\n", update, (double)maxFitness, sqrt(maxPhi), agent[who]->phi, agent[who]->correct, agent[who]->incorrect);

		if (evolvePhiPhase) {
			for(i=0;i<agent.size();i++)
			{
				tAgent *d;
				d=new tAgent;
				if(maxFitness<=0.0){
					j=rand()%(int)agent.size();
				} else
				do{ j=rand()%(int)agent.size(); } while((j==i)||(randDouble>(pow(agent[j]->phi,2)/maxPhi)));
				d->inherit(agent[j],perSiteMutationRate,update);
				nextGen[i]=d;
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
		if (maxPhi > evolvePhiLimit) {
			evolvePhiPhase = false;
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
	agent[0]->ancestor->ancestor->saveGenome(genomeFile);
//	agent[0]->ancestor->saveToDot(argv[3]);
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

