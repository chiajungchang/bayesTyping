#ifndef INCLUDE_BAYESTYPING_OPTIONS_H
#define INCLUDE_BAYESTYPING_OPTIONS_H

#include <seqan/arg_parse.h>
#include <seqan/sequence.h>
using namespace seqan;

#define VERSION "0.05"

struct PairAlleles
{
	int allele1;
	int allele2;
	int maxDiff;
	int score;

	PairAlleles(int in1,int in2,int in3, int in4){
		allele1 = in1;
		allele2 = in2;
		maxDiff = in3;
		score = in4;
	}
	bool operator<(const PairAlleles& rhs) const
    {
        return score < rhs.score;
    }

};

struct Options
{

    CharString readsFile; // name of the ccs read file
    CharString refFile;  // name of the allele file

    //Down Sampling
    bool downsampling;
    int downsamplingSize;
    CharString inSamplingFile;
    CharString outSamplingFile;

    //Align
    int match; //match score
    int mismatch; //mismatch score
    int affineOpen; //affineOpen
    int affineExtend; //affineExtend
    CharString inAlignFile;
    CharString outAlignFile;
    int frame; //border of the alignment
    int maxBand; // hidden
    double EPSILON;

    //bayestype
    CharString inClass;
    double noise;
    double filter_threshold;  //alleles with best-match reads < filter_threshold*#reads are filtered to accelerate
    CharString inTypingFile;
    CharString outTypingFile;


    //realign
    bool realignment;
    double bestHit;
    int maxNumAnswer;
    CharString inreAlignFile;
    CharString outreAlignFile;
    CharString outPrefix; //prefix to output allele1.fa allele2.fa

 //   CharString outPair; //candidate pairs and #reads
//    CharString inPrefix; //prefix to output allele1.fa allele2.fa


   //Function
   bool runAll;
   CharString configure;
   bool split;
   bool verbose;
   bool test;








    Options() : frame(30),filter_threshold(0.01),verbose(false),maxBand(60),test(false),downsampling(false),realignment(false),EPSILON(0.08),runAll(false),inClass(1)
    {}
};

seqan::ArgumentParser::ParseResult
parseCommandLine(Options & options, int argc, char const ** argv)
{
    seqan::ArgumentParser parser("bayesTyping");
    setShortDescription(parser, "");
	setVersion(parser, VERSION);
	addUsageLine(parser, "[\\fIOPTIONS\\fP] \\fB-s\\fP \\fIFILE\\fP  \\fB-r\\fP \\fIFILE\\fP");
    setCategory(parser, "HLA Typing");
    addDescription(parser,
                   "Align reads to HLA references as well as HLA typing");

	addSection(parser, "Required");
    addOption(parser, ArgParseOption("r", "reads", "CCS reads from the targeted gene of the sample.", ArgParseOption::INPUT_FILE, "IN"));
    setRequired(parser, "reads");
    setValidValues(parser, "reads", getFileExtensions(Fasta()));

    addOption(parser, ArgParseOption("g", "refgene", "HLA artificial genomic sequences of the target gene.", ArgParseOption::INPUT_FILE, "IN"));
    setRequired(parser, "refgene");
    setValidValues(parser, "refgene", getFileExtensions(Fasta()));


    addSection(parser, "Downsampling Options");
    addOption(parser, ArgParseOption("ds", "downsampling", "downsampling or not"));
    addOption(parser, ArgParseOption("sd", "downsamplingSize", "downsampling size", seqan::ArgParseOption::INTEGER, "INT"));
    setDefaultValue(parser, "downsamplingSize", "200");
    addOption(parser, ArgParseOption("id", "inSamplingFile", "input of sampling: IDs (from zero) of sampling sequences", ArgParseOption::INPUT_FILE, "IN"));
    addOption(parser, ArgParseOption("od", "outSamplingFile", "output of sampling: IDs (from zero) of sampling sequences", ArgParseOption::OUTPUT_FILE, "OUT"));


    addSection(parser, "Align Options");
     addOption(parser, ArgParseOption("m", "match", "Match score", seqan::ArgParseOption::INTEGER, "INT"));
    setDefaultValue(parser, "match", "5");
     addOption(parser, ArgParseOption("mi", "mismatch", "Mismatch penalty", seqan::ArgParseOption::INTEGER, "INT"));
    setDefaultValue(parser, "mismatch", "-6");
     addOption(parser, ArgParseOption("ao", "affineOpen", "AffineOpen penalty", seqan::ArgParseOption::INTEGER, "INT"));
    setDefaultValue(parser, "affineOpen", "-10");
     addOption(parser, ArgParseOption("ae", "affineExtend", "Match score", seqan::ArgParseOption::INTEGER, "INT"));
    setDefaultValue(parser, "affineExtend", "-3");
    addOption(parser, ArgParseOption("ia", "inAlignFile", "Saved alignment file to skip pairwise alignment", ArgParseOption::INPUT_FILE, "IN"));
    addOption(parser, ArgParseOption("oa", "outAlignFile", "Save alignment file", ArgParseOption::OUTPUT_FILE, "OUT"));

    addSection(parser, "Bayes Typing Options");
    addOption(parser, ArgParseOption("c", "class", "HLA Gene Class: either 1 or 2", seqan::ArgParseOption::STRING, "STRING"));
    setValidValues(parser, "class", "1 2");
    setDefaultValue(parser, "class", "1");
    addOption(parser, ArgParseOption("n", "noise", "Noise Level.", seqan::ArgParseOption::DOUBLE, "DOUBLE"));
    setDefaultValue(parser, "noise", "0.1");
    addOption(parser, ArgParseOption("ft", "filter_threshold", "alleles with best-match reads < filter_threshold*#reads are filtered to accelerate", seqan::ArgParseOption::DOUBLE, "DOUBLE"));
    setDefaultValue(parser, "filter_threshold", "0.01");
    addOption(parser, ArgParseOption("ib", "inTypingFile", "input of BayesTyping file: pairs of alleles, numbers, scores", ArgParseOption::INPUT_FILE, "IN"));
    addOption(parser, ArgParseOption("ob", "outTypingFile", "output of BayesTyping file: pairs of alleles, numbers, scores", ArgParseOption::OUTPUT_FILE, "OUT"));


	addSection(parser, "Re-alignment and filter options");
	addOption(parser, ArgParseOption("rf", "realignment", "realignment and filter"));
	addOption(parser, ArgParseOption("bh", "bestHit", "number of best hit", seqan::ArgParseOption::DOUBLE, "INT"));
    setDefaultValue(parser, "bestHit", "0.9");
    addOption(parser, ArgParseOption("ma", "maxNumAnswer", "maximum number of answers", seqan::ArgParseOption::INTEGER, "INT"));
    setDefaultValue(parser, "maxNumAnswer", "3");
    addOption(parser, ArgParseOption("ir", "inreAlignFile", "Saved realignment file to skip pairwise realignment", ArgParseOption::INPUT_FILE, "IN"));
    addOption(parser, ArgParseOption("or", "outreAlignFile", "Save realignment file", ArgParseOption::OUTPUT_FILE, "OUT"));

    addOption(parser, ArgParseOption("op", "outPrefix", "Prefix of output", seqan::ArgParseOption::STRING, "STRING"));

     addSection(parser, "Others");
	addOption(parser, ArgParseOption("all", "runAll", "run the entire pipeline"));
	addOption(parser, ArgParseOption("cf", "configure", "configure file for runAll", ArgParseOption::OUTPUT_FILE, "INT"));
    setDefaultValue(parser, "configure", "configure.txt");
	addOption(parser, ArgParseOption("sp", "split", "split reads by assigned genes"));
     addOption(parser, ArgParseOption("v", "verbose", "verbose"));
    addOption(parser, ArgParseOption("t", "test", "test"));


  //  addSection(parser, "Intermediate");
   // addOption(parser, ArgParseOption("ip", "inPrefix", "Prefix for allele1.fa, allele2.fa and noise.fa as input", ArgParseOption::STRING, "STRING"));


  //  addSection(parser, "Output");

    //	addOption(parser, ArgParseOption("o", "outPair", "Output filename. List the pair of predicted alleles and the number of supporting reads", ArgParseOption::OUTPUT_FILE, "OUT"));
 //   setDefaultValue(parser, "outPair", "out.txt");




	ArgumentParser::ParseResult res = seqan::parse(parser, argc, argv);
	getOptionValue(options.readsFile, parser, "reads");
	getOptionValue(options.refFile, parser, "refgene");

	options.downsampling = isSet(parser,"downsampling");
	getOptionValue(options.downsamplingSize, parser, "downsamplingSize");
	getOptionValue(options.inSamplingFile, parser, "inSamplingFile");
	getOptionValue(options.outSamplingFile, parser, "outSamplingFile");

	getOptionValue(options.match, parser, "match");
	getOptionValue(options.mismatch, parser, "mismatch");
	getOptionValue(options.affineExtend, parser, "affineExtend");
	getOptionValue(options.affineOpen, parser, "affineOpen");
    getOptionValue(options.inAlignFile, parser, "inAlignFile");
    getOptionValue(options.outAlignFile, parser, "outAlignFile");

	getOptionValue(options.inClass, parser, "class");
	getOptionValue(options.noise, parser, "noise");
	getOptionValue(options.filter_threshold, parser, "filter_threshold");
	getOptionValue(options.inTypingFile, parser, "inTypingFile");
    getOptionValue(options.outTypingFile, parser, "outTypingFile");

	options.realignment = isSet(parser,"realignment");
	getOptionValue(options.bestHit, parser, "bestHit");
	getOptionValue(options.maxNumAnswer, parser, "maxNumAnswer");
    getOptionValue(options.inreAlignFile, parser, "inreAlignFile");
    getOptionValue(options.outreAlignFile, parser, "outreAlignFile");
	 getOptionValue(options.outPrefix, parser, "outPrefix");


   getOptionValue(options.configure, parser, "configure");
   options.verbose = isSet(parser,"verbose");
	options.test = isSet(parser,"test");
	options.split = isSet(parser,"split");
	options.runAll = isSet(parser,"runAll");


	return res;
}

#endif
