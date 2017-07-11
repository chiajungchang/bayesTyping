#ifndef INCLUDE_BAYESTYPING_TYPING_H
#define INCLUDE_BAYESTYPING_TYPING_H

#include <vector>
#include <seqan/sequence.h>
#include <omp.h>
#include "options.h"



class Typing
{
private:
	int readlen;
	int reflen;
	int** scoreArr;
	vector<int> candAllele;
private:
	void calMaxHit();
	//void pickCandidateAllele(double filter_threshold);	
	
public:
	vector<PairAlleles> answers;
	vector<int> maxHit; // max scores for each reads gittest
	
public:
	void bayesTyping(int mnoise);
	int readTyping(CharString inTypingFile);
	int writeTyping(CharString outTypingFile);
	Typing(int **inScoreArr,int inReadlen, int inReflen);
	void pickCandidateAllele(double filter_threshold);
};

Typing::Typing(int **inScoreArr,int inReadlen, int inReflen){
	scoreArr = inScoreArr;
	readlen = inReadlen;
	reflen = inReflen;
	maxHit=vector<int>(readlen,0);
}

void Typing::calMaxHit(){ // max scores for each reads
	for(unsigned k=0;k<readlen;k++){
		int	maxscore=0;
		for(int ci=0;ci<candAllele.size();ci++){
			int i=candAllele[ci];
			if(scoreArr[k][i]>maxscore)
				maxscore=scoreArr[k][i];
		}
		if(maxscore>0){
			maxHit[k] = maxscore;
		}
	}
}
void Typing::bayesTyping(int mnoise){
	calMaxHit();	
	vector<int> candRead;
	vector<PairAlleles> tmpanswers;
	
	for(unsigned k=0;k<readlen;k++){
		if(maxHit[k]>0){
			candRead.push_back(k);
		}
	}

	int maxalign=0;
	
	
	//SEQAN_OMP_PRAGMA(parallel for)
	for(unsigned ci=0;ci<candAllele.size()-1;ci++){
		vector<int> difv;
		for(int cj=ci+1;cj<candAllele.size();cj++){
				int sum=0;
				int i = candAllele[ci];
				int j = candAllele[cj];
				difv.clear();
				for(int ck=0;ck<candRead.size();ck++){
					int k=candRead[ck];
					int maxscore=0;
					if(scoreArr[k][i]>scoreArr[k][j]){
						 maxscore=scoreArr[k][i];
					}
					else{
						maxscore=scoreArr[k][j];
					}
					if(maxHit[k]>maxscore)
						difv.push_back(maxHit[k]-maxscore);
					sum+=maxscore;
				}
				sort(difv.begin(),difv.end(),greater<int>());

				for(int k=0;k<mnoise && k<difv.size();k++)
					sum+=difv[k];
				

	//			#pragma omp critical(dataupdate)
				if(sum>=maxalign){ // prepare for the outputs
					maxalign = sum;
					if(mnoise==0)
						tmpanswers.push_back(PairAlleles(i,j,INT_MAX,sum));
					else
						tmpanswers.push_back(PairAlleles(i,j,difv[mnoise-1],sum));
					
				}
					
		}
	}
	for (vector<PairAlleles>::iterator it = tmpanswers.begin() ; it != tmpanswers.end(); ++it){
		if(it->score==maxalign)
			answers.push_back(*it);
	}
}

void Typing::pickCandidateAllele(double filter_threshold){
 
 	vector<int> matchReads(reflen,0); // #reads with max score to a allele

    int maxscore;
    
    for(unsigned k=0;k<readlen;k++){
        maxscore=0;
        for(unsigned i=0;i<reflen;i++){
            if(scoreArr[k][i]>maxscore)
                maxscore=scoreArr[k][i];
        }
        for(unsigned i=0;i<reflen;i++){
            if(scoreArr[k][i]==maxscore)
                matchReads[i]+=1;
        }
    }
	
    for(unsigned i=0;i<reflen;i++){
        if(matchReads[i]>filter_threshold*(double)readlen){
            candAllele.push_back(i);
        }
	}
	
}
int Typing::writeTyping(CharString outTypingFile){
	ofstream ofs (toCString(outTypingFile));
	if (ofs.is_open())
  	{
  		for(unsigned i=0;i<answers.size();i++){
			
			ofs << answers[i].allele1 << "\t" << answers[i].allele2 << "\t" << answers[i].maxDiff << "\t" << answers[i].score << endl;
			
		}				
		ofs.close();
		return 1;
	}
	else{
		return 0;
	}
}
int Typing::readTyping(CharString inTypingFile){
	ifstream ifs;
	ifs.open (toCString(inTypingFile), ifstream::in);
	answers.clear();
	if (ifs.is_open())
  	{
		int allele1,allele2,maxDiff,score;
		while ((ifs >> allele1 >> allele2 >> maxDiff >> score).good()) {
			answers.push_back(PairAlleles(allele1,allele2,maxDiff,score));
		}
		ifs.close();
		return 1;
	}
	else
		return 0;
}

#endif
