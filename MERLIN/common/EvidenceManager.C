#include <fstream>
#include <iostream>
#include <cstring>
#include <math.h>
#include "Error.H"
#include "Variable.H"
#include "VariableManager.H"
#include "Evidence.H"
#include "EvidenceManager.H"


EvidenceManager::EvidenceManager()
{
	foldCnt=1;
	preRandomizeSplit=false;
	randseed=0;
	dataMat = NULL;
}

EvidenceManager::~EvidenceManager()
{
	if (dataMat!=NULL)
	{
		delete dataMat;
	}
}

int
EvidenceManager::setVariableManager(VariableManager* aPtr)
{
	vMgr=aPtr;
	return 0;
}

Error::ErrorCode
EvidenceManager::loadEvidenceFromFile_Continuous(const char* inFName)
{
	ifstream inFile(inFName);
	char* buffer=NULL;
	string buffstr;
	int bufflen=0;
	int lineNo=0;

	// skip the first line.
	if(inFile.good())
	{
		getline(inFile,buffstr);
	}

	while(inFile.good())
	{
		getline(inFile,buffstr);

		if(buffstr.length()<=0)
		{
			continue;
		}
		if(bufflen<=buffstr.length())
		{
			if(buffer!=NULL)
			{
				delete[] buffer;
			}
			bufflen=buffstr.length()+1;
			buffer=new char[bufflen];
		}
		strcpy(buffer,buffstr.c_str());

		//All the evidences for each variable are stored in a map, indexed by the varId
		EMAP* evidMap=new EMAP;
		char* tok=strtok(buffer,"\t");
		//The toks take the form of varid and value

		int vId = 0;
		while(tok!=NULL)
		{
			Evidence* evid = new Evidence;
			evid->assocVariable(vId);
			//double varVal=log(atof(tok));
			double varVal=atof(tok);
			if(isinf(varVal) || isnan(varVal))
			{
				//cout <<"Found nan! " << tok << endl;
				cerr << "Please remove zero's (" << tok << ") from the expression data. " << endl;
				exit(-1);	
			}
			evid->setEvidVal(varVal);
			//(*evidMap)[vId]=evid;
			evidMap->push_back(evid);
			tok=strtok(NULL,"\t");
			vId++;
		}
		evidenceSet.push_back(evidMap);
		lineNo++;
	}

	inFile.close();
	//updateDataMat();

	cout <<"Read " << evidenceSet.size() << " different datapoints " << endl;

	return Error::SUCCESS;
}

//We create a matrix of randomized evidence, where each evidence has some value
//for the random variables. We populate the matrix one random variable at a time.
//We first generate a vector of permuted indices, in the range of 0 to the total
//number of evidences. Then we populate the part of the matrix associated with
//this variable by querying values from the original matrix in the order specified
//by the permuted indices

int
EvidenceManager::randomizeEvidence(gsl_rng* r)
{
	//First create all the evidence sets
	for(int i=0;i<evidenceSet.size();i++)
	{
		EMAP* evidMap=new EMAP;
		randEvidenceSet.push_back(evidMap);
	}
	//Populate variable wise
	VSET& variableSet=vMgr->getVariableSet();
	int* randInds=new int[trainIndex.size()];
	for(VSET_ITER vIter=variableSet.begin();vIter!=variableSet.end();vIter++)
	{
		//generate a random vector of indices ranging from 0 to evidenceSet.size()-1
		populateRandIntegers(r,randInds,trainIndex,trainIndex.size());	
		int j=0;
		for(int i=0;i<evidenceSet.size();i++)
		{
			EMAP* evidMap=NULL;
			if(trainIndex.find(i)!=trainIndex.end())
			{	
				evidMap=evidenceSet[randInds[j]];
				j++;
			}
			else
			{
				evidMap=evidenceSet[i];
			}
			EMAP* randEvidMap=randEvidenceSet[i];
			Evidence* evid=(*evidMap)[vIter->first];
			(*randEvidMap)[vIter->first]=evid;
		}
		string& geneName=(string&)vIter->second->getName();
		if((strcmp(geneName.c_str(),"FBgn0002631")==0) 
		|| (strcmp(geneName.c_str(),"FBgn0000411")==0) 
		|| (strcmp(geneName.c_str(),"FBgn0004915")==0) 
		|| (strcmp(geneName.c_str(), "FBgn0002573")==0)
		|| (strcmp(geneName.c_str(),"FBgn0005596")==0)  
		|| (strcmp(geneName.c_str(),"FBgn0035769")==0) 
		|| (strcmp(geneName.c_str(),"FBgn0011655")==0)
		|| (strcmp(geneName.c_str(),"FBgn0000576")==0))
		{
			cout <<geneName<<"IDs";
			for(int i=0;i<trainIndex.size();i++)
			{
				cout <<"\t" <<randInds[i];
			}
			cout << endl;
			cout <<geneName;
			for(INTINTMAP_ITER tIter=trainIndex.begin();tIter!=trainIndex.end();tIter++)
			{
				EMAP* emap=randEvidenceSet[tIter->first];
				cout <<"\t" << (*emap)[vIter->first]->getEvidVal();
			}
			cout << endl;
			cout <<geneName;
			for(INTINTMAP_ITER tIter=trainIndex.begin();tIter!=trainIndex.end();tIter++)
			{
				EMAP* emap=evidenceSet[tIter->first];
				cout <<"\t" << (*emap)[vIter->first]->getEvidVal();
			}
			cout << endl;
		}
	}
	return 0;
}

int 
EvidenceManager::getNumberOfEvidences()
{
	return evidenceSet.size();
}

EMAP* 
EvidenceManager::getEvidenceAt(int evId)
{
	if((evId>=evidenceSet.size()) && (evId<0))
	{
		return NULL;
	}
	return evidenceSet[evId];
}

EMAP* 
EvidenceManager::getRandomEvidenceAt(int evId)
{
	if((evId>=randEvidenceSet.size()) && (evId<0))
	{
		return NULL;
	}
	return randEvidenceSet[evId];
}

int
EvidenceManager::addEvidence(EMAP* evidSet)
{
	evidenceSet.push_back(evidSet);
	return 0;
}

int 
EvidenceManager::setFoldCnt(int f)
{
	foldCnt=f;
	return 0;
}

int
EvidenceManager::generateValidationSet(const char* vFName, int vSetSize,gsl_rng* r)
{
	ifstream inFile(vFName);
	if(inFile.good())
	{
		char buffer[256];
		while(inFile.good())
		{
			inFile.getline(buffer,255);
			if(strlen(buffer)<=0)
			{
				continue;
			}
			int dId=atoi(buffer);
			validationIndex[dId]=0;
		}
		inFile.close();
	}
	else
	{
		populateRandIntegers(r,validationIndex,evidenceSet.size(),vSetSize);
		ofstream oFile(vFName);
		for(INTINTMAP_ITER vIter=validationIndex.begin();vIter!=validationIndex.end();vIter++)
		{
			oFile << vIter->first << endl;
		}
		oFile.close();
	}
	return 0;
}


int 
EvidenceManager::setPreRandomizeSplit()
{
	preRandomizeSplit=true;
	return 0;
}

int 
EvidenceManager::splitData(int s)
{
	int testSetSize=(evidenceSet.size()-validationIndex.size())/foldCnt;
	int testStartIndex=s*testSetSize;
	int testEndIndex=(s+1)*testSetSize;
	if(s==foldCnt-1)
	{
		testEndIndex=evidenceSet.size()-validationIndex.size();
	}
	if(foldCnt==1)
	{
		testStartIndex=-1;
		testEndIndex=-1;
	}
	trainIndex.clear();
	testIndex.clear();
	int m=0;
	int* randInds=NULL;
	if(preRandomizeSplit)
	{
		randInds=new int[evidenceSet.size()];
		//generate a random vector of indices ranging from 0 to evidenceSet.size()-1
		gsl_rng* r=gsl_rng_alloc(gsl_rng_default);
		randseed=getpid();
		gsl_rng_set(r,randseed);
		populateRandIntegers(r,randInds,evidenceSet.size());	
		gsl_rng_free(r);
		cout <<"Random seed " << randseed << endl;
	}
	for(int i=0;i<evidenceSet.size();i++)
	{
		int eInd=i;
		if(randInds!=NULL)
		{
			eInd=randInds[i];
		}
		if(validationIndex.find(eInd)!=validationIndex.end())
		{
			continue;
		}
		if((m>=testStartIndex) && (m<testEndIndex))
		{
			testIndex[eInd]=0;
		}
		else
		{
			trainIndex[eInd]=0;
		}
		m++;
	}
	if(preRandomizeSplit)
	{
		delete[] randInds;
	}
	return 0;
}

INTINTMAP& 
EvidenceManager::getTrainingSet()
{
	return trainIndex;
}

INTINTMAP& 
EvidenceManager::getTestSet()
{
	return testIndex;
}

INTINTMAP&
EvidenceManager::getValidationSet()
{	
	return validationIndex;
}


int 
EvidenceManager::standardizeData()
{
	VSET& varSet=vMgr->getVariableSet();
	for(VSET_ITER vIter=varSet.begin();vIter!=varSet.end();vIter++)
	{
		double mean=0;
		for(int i=0;i<evidenceSet.size();i++)
		{
			EMAP* evidMap=evidenceSet[i];
			mean=mean+(*evidMap)[vIter->first]->getEvidVal();
		}
		mean=mean/evidenceSet.size();
		double std=0;
		for(int i=0;i<evidenceSet.size();i++)
		{
			EMAP* evidMap=evidenceSet[i];
			double diff=mean-(*evidMap)[vIter->first]->getEvidVal();
			std=std+(diff*diff);
		}
		std=sqrt(std/(evidenceSet.size()-1));
		//Standardize
		for(int i=0;i<evidenceSet.size();i++)
		{
			EMAP* evidMap=evidenceSet[i];
			Evidence* evid=(*evidMap)[vIter->first];
			double tval=evid->getEvidVal();
			double sval=(tval-mean)/std;
			evid->setEvidVal(sval);
		}
	}
	return 0;
}

int
EvidenceManager::partitionData(int numberOfComponents,map<int,EvidenceManager*>& evMgrSet,int& rseed,map<int,INTINTMAP*>& datasetInd)
{
	int datasubsetSize=evidenceSet.size()/numberOfComponents;
	int* randInds=new int[evidenceSet.size()];
	//generate a random vector of indices ranging from 0 to evidenceSet.size()-1
	gsl_rng* r=gsl_rng_alloc(gsl_rng_default);
	//randseed=getpid();
	randseed=18721;
	rseed=randseed;
	gsl_rng_set(r,randseed);
	populateRandIntegers(r,randInds,evidenceSet.size());	
	for(int i=0;i<evidenceSet.size();i++)
	{
		cout << randInds[i] << endl;
	}
	gsl_rng_free(r);
	int ind=0;
	for(int n=1;n<=numberOfComponents;n++)
	{
		int startIndex=(n-1)*(datasubsetSize);
		int endIndex=n*datasubsetSize;
		if(n==numberOfComponents)
		{
			endIndex=evidenceSet.size();
		}
		EvidenceManager* localManager=new EvidenceManager;
		INTINTMAP* origIDs=new INTINTMAP;
		int currind=(int)pow(2.0,ind);
		datasetInd[currind]=origIDs;
		evMgrSet[currind]=localManager;
		int dId=0;
		for(int e=startIndex;e<endIndex;e++)
		{
			//int rId=randInds[e];
			int rId=e;
			EMAP* evidSet=evidenceSet[rId];
			localManager->addEvidence(evidSet);
			(*origIDs)[dId]=rId;
			dId++;
		}
		localManager->updateDataMat();
		ind++;
	}
	delete[] randInds;
	return 0;
}

int 
EvidenceManager::populateRandIntegers(gsl_rng* r, int* randInds,int size)
{
	double step=1.0/(double)size;
	map<int,int> usedInit;
	for(int i=0;i<size;i++)
	{
		double rVal=gsl_ran_flat(r,0,1);
		int rind=(int)(rVal/step);
		while(usedInit.find(rind)!=usedInit.end())
		{
			rVal=gsl_ran_flat(r,0,1);
			rind=(int)(rVal/step);
		}
		usedInit[rind]=0;
		randInds[i]=rind;
	}
	usedInit.clear();
	return 0;
}


int 
EvidenceManager::populateRandIntegers(gsl_rng* r, int* randInds, INTINTMAP& populateFrom, int size)
{
	double step=1.0/(double)size;
	map<int,int> usedInit;
	map<int,int> temp;
	INTINTMAP_ITER tIter=populateFrom.begin();
	for(int i=0;i<size;i++)
	{
		int tid=tIter->first;
		temp[i]=tid;
		tIter++;
	}
	for(int i=0;i<size;i++)
	{
		double rVal=gsl_ran_flat(r,0,1);
		int rind=(int)(rVal/step);
		while(usedInit.find(rind)!=usedInit.end())
		{
			rVal=gsl_ran_flat(r,0,1);
			rind=(int)(rVal/step);
		}
		usedInit[rind]=0;
		randInds[i]=temp[rind];
	}
	usedInit.clear();
	return 0;
}


int 
EvidenceManager::populateRandIntegers(gsl_rng* r, INTINTMAP& randInds,int size, int subsetsize)
{
	double step=1.0/(double)size;
	for(int i=0;i<subsetsize;i++)
	{
		double rVal=gsl_ran_flat(r,0,1);
		int rind=(int)(rVal/step);
		while(randInds.find(rind)!=randInds.end())
		{
			rVal=gsl_ran_flat(r,0,1);
			rind=(int)(rVal/step);
		}
		randInds[rind]=0;
	}
	return 0;
}

int 
EvidenceManager::populateRandIntegers(gsl_rng* r, vector<int>& randInds,int size, int subsetsize)
{
	double step=1.0/(double)size;
	map<int,int> usedInit;
	for(int i=0;i<subsetsize;i++)
	{
		double rVal=gsl_ran_flat(r,0,1);
		int rind=(int)(rVal/step);
		while(usedInit.find(rind)!=usedInit.end())
		{
			rVal=gsl_ran_flat(r,0,1);
			rind=(int)(rVal/step);
		}
		usedInit[rind]=0;
		randInds.push_back(rind);
	}
	return 0;
}

int
EvidenceManager::estimateCovariance(int uId, int vId, double& ucov, double& vcov, double& uvcov)
{
	gsl_vector_view Du = dataMat->getRowView(uId);
	gsl_vector_view Dv = dataMat->getRowView(vId);
	ucov  = gsl_stats_covariance (Du.vector.data, Du.vector.stride, Du.vector.data, Du.vector.stride, Du.vector.size); 
	vcov  = gsl_stats_covariance (Dv.vector.data, Dv.vector.stride, Dv.vector.data, Dv.vector.stride, Dv.vector.size); 
	uvcov = gsl_stats_covariance (Du.vector.data, Du.vector.stride, Dv.vector.data, Dv.vector.stride, Du.vector.size); 
	return 0;
}

int
EvidenceManager::updateDataMat()
{
	if (dataMat != NULL)
	{
		delete dataMat;
	}
	int nSampels = evidenceSet.size();
	int nGenes   = evidenceSet[0]->size();
	dataMat = new Matrix(nGenes,nSampels);
	for (int j=0;j<nSampels;j++)
	{
		EMAP* evidMap = evidenceSet[j];
		//for (map<int,Evidence*>::iterator itr=evidMap->begin(); itr!=evidMap->end(); itr++)
		for (int i=0;i<evidMap->size();i++)
		{
			//int i=itr->first;
			//Evidence* evid = itr->second;
			Evidence* evid = evidMap->at(i);
			dataMat->setValue(evid->getEvidVal(),i,j);
		}
	}
	return 0;
}
