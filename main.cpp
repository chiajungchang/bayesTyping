#include <seqan/arg_parse.h>
#include <seqan/seq_io.h>
#include <iostream>
#include <fstream>
#include <omp.h>
#include <string>
#include <seqan/store.h>
#include "bayesTyping.h"
#include "options.h"


using namespace std;
using namespace seqan;


int main(int argc, const char* argv[])
{
    // Parse the command line.
    Options options;;
    ArgumentParser::ParseResult res = parseCommandLine(options, argc, argv);
    if (res != ArgumentParser::PARSE_OK)
        return res == ArgumentParser::PARSE_ERROR;

    HlaStore *hlastore;
    Typing *typing;


    if(!options.runAll){
		if(options.split){
			if(options.verbose)
				cout << "Begin splitting the CCS file..." << endl;
			hlastore = new HlaStore(options.refFile, options.readsFile, true);
    		hlastore->assignGene(options);
    		if(options.verbose)
				cout << "Writing the CCS files with prefix: " << options.outPrefix << endl;
    		hlastore->outSplitCCS(options);
    		delete hlastore;

    		if(options.verbose)
				cout << "End splitting the CCS file." << options.outPrefix << endl;
    		return 0;
    	}


    	//downsampling
    	hlastore = new HlaStore(options.refFile, options.readsFile, false);

    	if(options.downsampling){
    		if(options.inSamplingFile!=""){
    			if(options.verbose)
					cout << "Reading sampling file from " << options.inSamplingFile << endl;
				hlastore->readSampling(options.inSamplingFile);
    		}
    		else{
    			if(options.verbose){
					cout << "Sampling..." << endl;
					cout << "Original size: " << hlastore->readSize()  << endl;
					cout << "Sample size: " << options.downsamplingSize << endl;
				}
				hlastore->randSampling(options.downsamplingSize);

    		}
    		hlastore->downSampling();
    		if(options.outSamplingFile!=""){

    			if(options.verbose)
					cout << "Output sampling file to " << options.outSamplingFile << endl;
				hlastore->writeSampling(options.outSamplingFile);
    		}
    	}

    	//alignment

    	if(options.inAlignFile!=""){
    		if(options.verbose)
					cout << "Reading alignment file from " << options.inAlignFile << endl;
			if(!hlastore->readAlign(options.inAlignFile)){
				cerr << "Cannot read alignment file: " << options.inAlignFile << endl;
				return 0;
			}
    	}
    	else{
    		if(options.verbose){
					cout << "Pairwise alignment...  It might take longer time."  << endl;
					cout << "#reads: " << hlastore->readSize() << endl;
					cout << "#alleles: " << hlastore->refSize() << endl;
			}
			hlastore->assignGroup(options.inClass);
			hlastore->pwAlign(options);

			if(options.verbose){
					cout << "Pairwise alignments are finished."  << endl;

			}

    	}
    	if(options.outAlignFile!=""){

    		if(options.verbose)
					cout << "Output alignment file to " << options.outAlignFile << endl;
    		if(!hlastore->writeAlign(options.outAlignFile)){
				cerr << "Cannot write to alignment file: " << options.outAlignFile << endl;
				return 0;
			}
    	}

    	//bayesTyping
    	typing = new Typing(hlastore->scoreArr,hlastore->readSize(),hlastore->refSize());

    	if(options.inTypingFile!=""){
    		if(options.verbose)
					cout << "Reading bayesTyping results from " << options.inTypingFile << endl;

			if(!typing->readTyping(options.inTypingFile)){
				cerr << "Cannot read typing file:" << options.inTypingFile << endl;
				return 0;
			}
    	}
    	else{
    		if(options.verbose)
					cout << "Running the BayesTyping method."  << endl;
			typing->pickCandidateAllele(options.filter_threshold);
			typing->bayesTyping(options.noise*(double)hlastore->readSize());


			if(options.verbose)
					cout << "BayesTyping isFinished"  << endl;

    	}
    	if(options.outTypingFile!=""){

    		if(options.verbose)
					cout << "Output bayesTyping results to " << options.outTypingFile << endl;
    		if(!typing->writeTyping(options.outTypingFile)){
				cerr << "Cannot write to typing file: " << options.outTypingFile << endl;
				return 0;
			}
    	}



    	if(options.realignment){
    		if(options.verbose)
					cout << "Begin re-alignment..."  << endl;


			for(unsigned i=0;i< typing->answers.size() && i<options.maxNumAnswer;i++){
				std::stringstream ss;
				ss << i+1;
				CharString str = ss.str();
				delete(hlastore);
				hlastore = new HlaStore(options.refFile, options.readsFile, false);
				hlastore->onlyTwo(typing->answers[i]);
				if(options.inreAlignFile!=""){

					CharString inreAlignFile=options.inreAlignFile;
					append(inreAlignFile,str);
					if(options.verbose)
						cout << "Reading realignment file: " << inreAlignFile << endl;

					if(!hlastore->readAlign(inreAlignFile)){
						cerr << "Cannot read realignment " << inreAlignFile << endl;
						return 0;
					}
				}
				else{

					hlastore->pwAlign(options);
				}
				if(options.outreAlignFile!=""){
					CharString outreAlignFile=options.outreAlignFile;
					append(outreAlignFile,str);
					if(options.verbose)
						cout << "Writing realignment file: " << outreAlignFile << endl;
					if(!hlastore->writeAlign(outreAlignFile)){
						cerr << "Cannot write to realignment file: " << outreAlignFile << endl;
						return 0;
					}

				}

				CharString outPrefix=options.outPrefix;
				append(outPrefix,str);
				hlastore->filterNPrint(outPrefix,options.bestHit);
			}


			if(options.verbose)
					cout << "Begin filtering..."  << endl;

    	}
    	else{
    		hlastore->printPair(typing->answers,options.outPrefix,options.maxNumAnswer);
    	}



    }
    else{
    	if(options.verbose)
					cout << "Running All... " << endl;

    }



    if(options.verbose)
					cout << "End of the program."  << endl;

	return res;
}
