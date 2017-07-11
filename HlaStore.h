#include <seqan/basic.h>
#include <seqan/align.h>
#include <seqan/modifier.h>
#include <seqan/arg_parse.h>
#include <seqan/seq_io.h>
#include <seqan/stream.h>
#include <omp.h>
#include <limits.h>
#include <string>
#include <sstream>
#include <seqan/store.h>
#include <seqan/consensus.h>
#include <seqan/random.h>
#include "options.h"
//#include "Typing.h"


using namespace std;
using namespace seqan;


typedef Index<StringSet<Dna5String>, IndexQGram<Shape<Dna, UngappedShape<11> >, OpenAddressing> > TIndex;
typedef Pattern<TIndex, Swift<SwiftSemiGlobal> > TPattern;
typedef Finder<Dna5String, Swift<SwiftSemiGlobal> > TFinder;
typedef Align<Dna5String, ArrayGaps> TAlign;

struct BestHit
{
	int bestAllele;
	int readid;
	int score;
	int score2;

};


class HlaStore{


private:
	StringSet<CharString> readids;
	StringSet<Dna5String> readseqs;
	StringSet<CharString> refids;
	StringSet<Dna5String> refseqs;
	Align<Dna5String> refalign;
	std::set<unsigned> sampledIds;
	int* alleleGroup;
	int* assignedGene;
	//int** scoreArr;
	bool apxMatch;
	double EPSILON;

public:
	int** scoreArr;

private:

	Dna5 getRevCompl(Dna5 const & nucleotide);
	Dna5String getRevComplString(Dna5String genome);
	void approMatch(TPattern &pattern,Dna5String ref,int *maxGeneScore,int i);


public:
	HlaStore(CharString refFile, CharString readsFile, bool inApxMatch);
	~HlaStore();
	void assignGene(Options & options);
	void outSplitCCS(Options & options);
	void downSampling();
	void randSampling(int downsamplingSize);
	int readSampling(CharString inSamplingFile);
	int writeSampling(CharString outSamplingFile);
	void pwAlign(Options & options);
	int readAlign(CharString file);
	int writeAlign(CharString file);
	int readSize();
	int refSize();
	void assignGroup(CharString inClass);
	void printPair(vector<PairAlleles> answers,CharString outPrefix, int maxNumAns);
	void onlyTwo(PairAlleles answer);
	void filterNPrint(CharString outPrefix, double bestHitP);

};


bool sortBestHit(const BestHit &lhs, const BestHit &rhs){ return lhs.score > rhs.score; };

void HlaStore::filterNPrint(CharString outPrefix, double bestHitP){
	int bestHit = length(readids) * bestHitP;

	vector<BestHit> bh(length(readids));
	for(unsigned k=0;k<length(readids);k++){
		bh[k].score=0;
		bh[k].readid=k;
		for(unsigned i=0;i<length(refids);i++){
			if(scoreArr[k][i]>bh[k].score){
				bh[k].score2 = bh[k].score;
				bh[k].score=scoreArr[k][i];
				bh[k].bestAllele = i;
			}
			else{
				bh[k].score2=scoreArr[k][i];
			}
		}
	}
	sort(bh.begin(),bh.end(),sortBestHit);

	vector<int> reads0;
    vector<int> reads1;
    pair<double,double> paircount;
	stringstream resultname;
    resultname << outPrefix << "_result.txt";
    ofstream ofs (resultname.str().c_str());
	if (!ofs.is_open()){
		cerr << "Cannot open output file: " << resultname << endl;
		return;
	}
	int score=0;
	for(unsigned k=0;k<length(readids) && k< bestHit;k++){

		score+=bh[k].score;
		if(bh[k].score == bh[k].score2){
			cout << bh[k].score << endl;
			paircount.first+=0.5;
			paircount.second+=0.5;
		}
		else{
			if(bh[k].bestAllele==0){
				paircount.first+=1;
				reads0.push_back(bh[k].readid);
			}
			else{

				paircount.second+=1;
				reads1.push_back(bh[k].readid);
			}
		}
	}
	ofs << "allele1\tallele2\t#allele1_reads\t#allele2_reads\tscore" << endl;
	ofs << refids[0] << "\t" << refids[1] << "\t" << paircount.first << "\t" << paircount.second << "\t" << score << endl;
	ofs.close();
	stringstream ss1;
		stringstream ss2;
		stringstream ssref1;
		stringstream ssref2;
		ss1 << outPrefix  << "_allele1.fa";
		ss2 << outPrefix  << "_allele2.fa";
		ssref1 << outPrefix  << "_ref1.fa";
		ssref2 << outPrefix << "_ref2.fa";
		SeqFileOut falign1(ss1.str().c_str());
		SeqFileOut falign2(ss2.str().c_str());
		SeqFileOut fref1(ssref1.str().c_str());
		SeqFileOut fref2(ssref2.str().c_str());
		for(vector<int>::iterator it = reads0.begin() ; it != reads0.end(); ++it){
			writeRecord(falign1,readids[*it],readseqs[*it] );
		}
		for(vector<int>::iterator it = reads1.begin() ; it != reads1.end(); ++it){
			writeRecord(falign2,readids[*it],readseqs[*it] );
		}
		writeRecord(fref1,refids[0],refseqs[0] );
		writeRecord(fref2,refids[1],refseqs[1] );

}

void HlaStore::printPair(vector<PairAlleles> answers,CharString outPrefix, int maxNumAns){
	vector<int> maxHit(length(readids),0);

	for(unsigned k=0;k<length(readids);k++){
		int	maxscore=0;
		for(unsigned i=0;i<length(refids);i++){
			if(scoreArr[k][i]>maxscore)
				maxscore=scoreArr[k][i];
		}
		if(maxscore>0){
			maxHit[k] = maxscore;
		}
	}
	vector<int> reads0;
    vector<int> reads1;
    pair<double,double> paircount;
    stringstream resultname;
    resultname << outPrefix << "_result.txt";
    ofstream ofs (resultname.str().c_str());
	if (!ofs.is_open()){
		cerr << "Cannot open output file: " << resultname << endl;
		return;
	}

    ofs << "allele1\tallele2\t#allele1_reads\t#allele2_reads\tscore" << endl;

	for(unsigned a=0;a<answers.size() && a<maxNumAns;a++){
		int i = answers[a].allele1;
		int j = answers[a].allele2;
		paircount.first=0;
		paircount.second=0;
		reads0.clear();
		reads1.clear();
		for(int k=0;k<length(readids);k++){
			if(maxHit[k]<=0)
				continue;
			if(scoreArr[k][i]>scoreArr[k][j]){
				if((maxHit[k]-scoreArr[k][i])<=answers[a].maxDiff){
					paircount.first+=1;
					reads0.push_back(k);
				}
			}
			else{
				if((maxHit[k]-scoreArr[k][j])<=answers[a].maxDiff){
					if(scoreArr[k][i]==scoreArr[k][j]){
						paircount.first+=0.5;
						paircount.second+=0.5;
                     }
					else{
						paircount.second+=1;
						reads1.push_back(k);
					}
				}

			}
		}
		ofs << refids[i] << "\t" << refids[j] << "\t" << reads0.size() << "\t" << reads1.size() << "\t" << answers[a].score << endl;

		stringstream ss1;
		stringstream ss2;
		stringstream ssref1;
		stringstream ssref2;
		ss1 << outPrefix << (a+1) << "_allele1.fa";
		ss2 << outPrefix << (a+1) << "_allele2.fa";
		ssref1 << outPrefix << (a+1) << "_ref1.fa";
		ssref2 << outPrefix << (a+1) << "_ref2.fa";
		SeqFileOut falign1(ss1.str().c_str());
		SeqFileOut falign2(ss2.str().c_str());
		SeqFileOut fref1(ssref1.str().c_str());
		SeqFileOut fref2(ssref2.str().c_str());
		for(vector<int>::iterator it = reads0.begin() ; it != reads0.end(); ++it){
			writeRecord(falign1,readids[*it],readseqs[*it] );
		}
		for(vector<int>::iterator it = reads1.begin() ; it != reads1.end(); ++it){
			writeRecord(falign2,readids[*it],readseqs[*it] );
		}
		writeRecord(fref1,refids[answers[a].allele1],refseqs[answers[a].allele1] );
		writeRecord(fref2,refids[answers[a].allele2],refseqs[answers[a].allele2] );

	}
	ofs.close();
}

void HlaStore::onlyTwo(PairAlleles answer){

	StringSet<CharString> tmpids;
	StringSet<Dna5String> tmprefs;
	appendValue(tmpids,refids[answer.allele1]);
    appendValue(tmprefs,refseqs[answer.allele1]);
    appendValue(tmpids,refids[answer.allele2]);
    appendValue(tmprefs,refseqs[answer.allele2]);
    refids=tmpids;
    refseqs=tmprefs;
	alleleGroup[0] = 1;
	alleleGroup[1] = 2;

}

HlaStore::HlaStore(CharString refFile, CharString readsFile, bool inApxMatch){



	SeqFileIn refseqFileIn(toCString(refFile));


	SeqFileIn readseqFileIn(toCString(readsFile));
	readRecords(readids, readseqs, readseqFileIn);

	apxMatch = inApxMatch;

	if(apxMatch){
        readRecords(refids, refseqs, refseqFileIn);
		assignedGene = new int[length(readseqs)];
	}
	else{

        typedef Iterator<StringSet<CharString >, Standard>::Type TCIter;
        typedef Row< Align<Dna5String> >::Type Trow;
        StringSet<CharString> tmpset;
        readRecords(refids, tmpset, refseqFileIn);

        int alignLength = length(*(begin(tmpset, Standard())));
        resize(rows(refalign),length(tmpset));
        int j=0;
        for (TCIter it = begin(tmpset, Standard()); it != end(tmpset, Standard()); ++it){
            string s(toCString(*it));
            if(s.length()!=alignLength){
                cerr << "Alignments  should have the same length" << j << " " << s.length() << " " << alignLength << endl;
            }
            string result;
            remove_copy(s.begin(), s.end(), back_inserter(result), '-');
            appendValue(refseqs, result.c_str());
            assignSource(row(refalign,j),result.c_str());
            Trow & rowj = row(refalign,j);
            int sindex=0;
            for (string::iterator sit=s.begin(); sit!=s.end(); ++sit){
                if((*sit)=='-'){

                    insertGap(rowj,sindex);
                }
                sindex++;
            }
            j++;
        }

		alleleGroup = new int[length(refseqs)];
		scoreArr = new int*[length(readseqs)];
		for(int i = 0; i < length(readseqs); ++i)
   		 	scoreArr[i] = new int[length(refseqs)];
	}


}

HlaStore::~HlaStore(){
	if(apxMatch){
		delete assignedGene;
	}
	else{
		delete alleleGroup;
		delete[] scoreArr;
	}

}

int HlaStore::readSize(){
	return length(readseqs);
}
int HlaStore::refSize(){
	return length(refseqs);
}

void HlaStore::assignGroup(CharString inClass){
	int curGroup = 1;
	CharString curGroupName;
	if(inClass=="1"){
		for(int i=0;i<length(refids);i++){
			alleleGroup[i] = curGroup;
		}
	}
	else if(inClass=="2"){
		for(unsigned j=0;j<length(refids[0]);j++){
			if(refids[0][j]==':'){
				curGroupName = prefix(refids[0],j);
				break;
			}
		}
		alleleGroup[0]=curGroup;
		for(int i=1;i<length(refids);i++){
			for(unsigned j=0;j<length(refids[i]);j++){
				if(refids[i][j]==':'){
					if(prefix(refids[i],j)!=curGroupName){
						curGroupName = prefix(refids[i],j);
						curGroup++;
					}
					break;
				}
			}
			alleleGroup[i]=curGroup;
		}

	}
}

Dna5 HlaStore::getRevCompl(Dna5 const & nucleotide)
{
    if (nucleotide == (Dna5)'A')
        return (Dna5)'T';
    if (nucleotide == (Dna5)'T')
        return (Dna5)'A';
    if (nucleotide == (Dna5)'C')
        return (Dna5)'G';

	if (nucleotide == (Dna5)'G')
        return (Dna5)'C';
    return (Dna5)'N';
}

Dna5String HlaStore::getRevComplString(Dna5String genome){
	Dna5String revComplGenome;
    resize(revComplGenome, length(genome));

    for (unsigned i = 0; i < length(genome); ++i)
    {
        revComplGenome[length(genome) - 1 - i] = getRevCompl(genome[i]);
    }
    return revComplGenome;
}


int HlaStore::readAlign(CharString file){
	ifstream ifs;
	ifs.open (toCString(file), ifstream::in);
	if (ifs.is_open())
  	{
		int i,j,score;
		while (ifs.good()) {
			ifs >> i >> j >> score;
			scoreArr[i][j]=score;
		}
		ifs.close();
		return 1;
	}
	else
		return 0;
}
int HlaStore::writeAlign(CharString file){
	ofstream ofs (toCString(file));
	if (ofs.is_open())
  	{
		 for(int i=0;i<length(readids);i++)
			for(int j=0;j<length(refids);j++)
				ofs << i << "\t" << j << "\t" << scoreArr[i][j] << endl;

		ofs.close();
		return 1;
	}
	else{
		return 0;
	}

}


void HlaStore::approMatch(TPattern &pattern,Dna5String ref,int *maxGeneScore,int i){
	TFinder finder(ref);
	while (find(finder, pattern, EPSILON))
    {
            // Verify match.

			Finder<Dna5String> verifyFinder(ref);
			if(beginPosition(finder)>=0)
            	setPosition(verifyFinder, beginPosition(finder));
            int j=position(pattern).i1;
           	Pattern<Dna5String, MyersUkkonen> verifyPattern(readseqs[j]);
            unsigned readLength = length(readseqs[j]);
            int minScore = -static_cast<int>(EPSILON * readLength * 2.0);
            setScoreLimit(verifyPattern, minScore);
            int maxs=-100000;
    		while (find(verifyFinder, verifyPattern))
    		{
    			if(getScore(verifyPattern)>maxs)
    				maxs=getScore(verifyPattern);
       		 }
   			 if(maxs>maxGeneScore[j]){

   			 	maxGeneScore[j]=maxs;
   			 	assignedGene[j]=i;
   			 }
     }
}

void HlaStore::assignGene(Options & options){
	EPSILON = options.EPSILON;
	FragmentStore<> fragStore;
	int* maxGeneScore = new int[length(readseqs)];
	for(int i=0;i<length(readseqs);++i)
	{
		assignedGene[i]=-1;
		maxGeneScore[i]=-100000;
	}

	TIndex index(readseqs);
	TPattern pattern(index);
	for(int i=0;i<length(refseqs);++i)
	{
		approMatch(pattern,refseqs[i],maxGeneScore,i);
		approMatch(pattern,getRevComplString(refseqs[i]),maxGeneScore,i);

	}
	if(options.downsampling){
		for(int i=0;i<length(refseqs);++i)
		{
			priority_queue<int> mypq;
			int threshold;
			for(int j=0;j<length(readseqs);++j)
			{
				if(assignedGene[j]==i)
					mypq.push(maxGeneScore[j]);

			}
			int num = options.downsamplingSize;
			while(!mypq.empty() && num>0){
				threshold = mypq.top();
				mypq.pop();
				num--;
			}


			for(int j=0;j<length(readseqs);++j)
			{
				if(assignedGene[j]==i && maxGeneScore[j]>=threshold)
					num--;
				if((assignedGene[j]==i && maxGeneScore[j]<threshold)||num<0)
					assignedGene[j]=-1;

			}
			mypq = priority_queue <int>();

		}

	}
	delete maxGeneScore;
}

void HlaStore::outSplitCCS(Options & options){
	vector<CharString> genes;
	int* refidx = new int[length(refids)];

	CharString curGroupName;
	curGroupName = refids[0];
	genes.push_back(curGroupName);
	int j=0;
	refidx[0] = j;

	for(int i=1;i<length(refids);++i)
	{
		if(curGroupName!=refids[i]){
			curGroupName = refids[i];
			genes.push_back(curGroupName);
			j++;
		}
		refidx[i]=j;
		cout << i << " " << j << endl;
	}
	for(int i=0;i<length(readids);++i)
	{
		if(assignedGene[i]>=0)
			assignedGene[i]=refidx[assignedGene[i]];
	}
	for(int i=0;i<genes.size();i++){
		stringstream ss1;
		ss1 << options.outPrefix << genes[i] << ".fa";
		SeqFileOut ffa(ss1.str().c_str());
		for(j=0 ; j<length(readids); j++){
			if(assignedGene[j]==i)
				writeRecord(ffa,readids[j],readseqs[j] );
		}
	}
}
void HlaStore::randSampling(int downsamplingSize){

	if(downsamplingSize > length(readids)){
		cout << "Request to sample more reads than there actually are! Use all reads." << endl;
		for(unsigned i=0;i<length(readids);i++){
			sampledIds.insert(i);
		}
		return;
	}

	Rng<MersenneTwister> rng(816);
    Pdf<Uniform<unsigned> > pdf(0, length(readids) - 1);

    while (sampledIds.size() < downsamplingSize)
	{
		unsigned x = pickRandomNumber(rng, pdf);
		sampledIds.insert(x);
	}


}
void HlaStore::downSampling(){
	StringSet<CharString> tmpids;
	StringSet<Dna5String> tmpreads;
	std::set<unsigned>::iterator it = sampledIds.begin();
    for (; it != sampledIds.end(); ++it){
    	appendValue(tmpids,readids[*it]);
    	appendValue(tmpreads,readseqs[*it]);
    }
    readids=tmpids;
    readseqs=tmpreads;

}
int HlaStore::readSampling(CharString inSamplingFile){
	ifstream ifs;
	ifs.open (toCString(inSamplingFile), ifstream::in);
	if (ifs.is_open())
  	{
		unsigned x;
		while (ifs.good()) {
			ifs >> x;
			sampledIds.insert(x);
		}
		ifs.close();
		return 1;
	}
	else
		return 0;
}
int HlaStore::writeSampling(CharString outSamplingFile){
	ofstream ofs (toCString(outSamplingFile));
	if (ofs.is_open())
  	{
  		std::set<unsigned>::iterator it = sampledIds.begin();
    	for (; it != sampledIds.end(); ++it){
				ofs << *it << endl;
		}

		ofs.close();
		return 1;
	}
	else{
		return 0;
	}
}

int recalAlignScore(TAlign align, int mismatch_score,int gap_score)
{
    typedef Row<TAlign>::Type TRow;
    typedef Iterator<TRow>::Type TRowIterator;
    TRowIterator it0 = begin(row(align,0)),
                it1=begin(row(align,1)),
                itEnd=end(row(align,0));
    int score=0;
    for(;it0!=itEnd;it0++,it1++)
    {

        if(isGap(it0) || isGap(it1))
        {
            score += gap_score;

        }
        else if(*it0 != *it1)
        {
          //  cout << "mismatch:" << *it0 << " " << *it1 << endl;
            score += mismatch_score;
        }

    }
    return score;

}
void HlaStore::pwAlign(Options & options){



	SEQAN_OMP_PRAGMA(parallel for)
	for(int i=0;i<length(readseqs);++i)
	{
		if(options.verbose){
			cout << "." ;
			cout.flush();
		}

		int *tscores= new int[length(refseqs)];
		Dna5String matchString;
		int boundRation=1;
		int curGroup=0;
		int groupj=0;
		int basescore=0;
		for(int j=0;j<length(refseqs);++j)
		{
        //    if(j==237)
         //       j=j+1;
			int score;
			int lbound;
			int ubound;
			int reScore;
			if(alleleGroup[j]!=curGroup){
				curGroup=alleleGroup[j];
				groupj=j;
				TAlign falign;
				resize(rows(falign), 2);
				assignSource(row(falign, 0), refseqs[j]);
    			assignSource(row(falign, 1), readseqs[i]);
    			int fscore =localAlignment(falign, Score<int, Simple>(options.match, options.mismatch, options.affineExtend,options.affineOpen));

    			matchString=getRevComplString(readseqs[i]);
    			TAlign ralign;
				resize(rows(ralign), 2);
				assignSource(row(ralign, 0), refseqs[j]);
    			assignSource(row(ralign, 1), matchString);
    			int rscore = localAlignment(ralign, Score<int, Simple>(options.match, options.mismatch, options.affineExtend,options.affineOpen));

    			if(fscore>rscore){
    				score=fscore;
    				matchString=infix(readseqs[i], beginPosition(row(falign,1)), endPosition(row(falign,1)))  ;

                    lbound = beginPosition(row(falign,0));
                    ubound = endPosition(row(falign,0));

					if(options.test)
   		 				cout << falign << endl;
    			}
    			else{
    				score=rscore;
    				matchString=infix(matchString, beginPosition(row(ralign,1)), endPosition(row(ralign,1)))  ;

    				 lbound = beginPosition(row(ralign,0));
                    ubound = endPosition(row(ralign,0));

					if(options.test)
    					cout << ralign << endl;

    			}
				if(options.test)
					cout << "Score:" << score << endl;
				if(lbound<0)
					lbound=0;
				if(ubound>length(refseqs[j]))
					ubound=length(refseqs[j]);
				if(options.test){

					cout << "Trimmed " << infix(refseqs[j],lbound,ubound) << endl;
					cout << "Original Ref " << refseqs[j] << endl;
				}
				score = globalAlignmentScore(infix(refseqs[j],lbound,ubound),matchString, Score<int, Simple>(options.match, options.mismatch, options.affineExtend,options.affineOpen), AlignConfig<false, false, false, false>()); //AlignConfig<true,false, false, true>

    			for(boundRation=2;boundRation<options.maxBand;boundRation=boundRation*3/2){
    				int bandscore = globalAlignmentScore(infix(refseqs[j],lbound,ubound),matchString, Score<int, Simple>(options.match, options.mismatch, options.affineExtend,options.affineOpen), AlignConfig<false, false, false, false>(),-1*boundRation*options.frame,boundRation*options.frame);
					if(bandscore==score){
						break;
					}
    			}



    			if(boundRation>=options.maxBand){
    				score=0;
    			}
    			else{
                    TAlign align;
                    resize(rows(align), 2);
                    assignSource(row(align, 0), infix(refseqs[j],lbound,ubound));
                    assignSource(row(align, 1), matchString);
                    globalAlignment(align, Score<int, Simple>(options.match, options.mismatch, options.affineExtend,options.affineOpen), AlignConfig<false, false, false, false>(),-1*boundRation*options.frame,boundRation*options.frame);
                    reScore = recalAlignScore(align,options.mismatch,options.affineOpen);
                    basescore=score-reScore;

    			}
			}
			else{
				if(boundRation>=options.maxBand){
    				score=0;
    			}
    			else{
                    int curlbound=toSourcePosition(row(refalign,j),toViewPosition(row(refalign,groupj),lbound));
                    int curubound=toSourcePosition(row(refalign,j),toViewPosition(row(refalign,groupj),ubound));
                    if(options.test){
                        cout << refids[j] << endl;
                        cout << "Trimmed " << infix(refseqs[j],curlbound,curubound) << endl;

                    }
                    TAlign align;
                    resize(rows(align), 2);
                    assignSource(row(align, 0), infix(refseqs[j],curlbound,curubound));
                    assignSource(row(align, 1), matchString);
					score = globalAlignment(align, Score<int, Simple>(options.match, options.mismatch, options.affineExtend,options.affineOpen), AlignConfig<false, false, false, false>(),-1*boundRation*options.frame,boundRation*options.frame);
                    if(score>0)
                        score = basescore + recalAlignScore(align,options.mismatch,options.affineOpen);
				}

			}
			tscores[j]=score;

    	}

    	for(int j=0;j<length(refseqs);++j)
    			scoreArr[i][j]=tscores[j];

    	delete[] tscores;

	}
}


