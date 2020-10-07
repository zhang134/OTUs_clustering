#include "mcc.h"

using namespace std;
Options options;
SequenceDB seq_db;
#pragma warning(disable:4996)
int main(int argc, const char *argv[])
{
	//int argc = 5;
	if (argc < 5) {
		print_usage();
		return 0;
	}
	string db_in;
	string genome_in;
	string db_out;
	string command;
	//clock_t start, finish;
	//start = clock();
	//double totaltime;
	float identify;
	//int argc = 8;
	int max_length = 0, max_banded_width = 0, kmerLength;
	//int **matrix; // score matrix
	//int **track; // track score
	//const char *argv[] = { "./mcc", "-i" ,"D:\\OTU数据\\前100条.txt", "-o" ,"V6SL_sim.txt" };
	//const char *argv[] = { "./smsAlign","-locate","-i","D:\\OTU拼接\\回国后文章\\数据\\ecoli真实数据\\ecoliCLRreads.fa", "-lib", "ecoili_raw_lib.txt","-o" ,"ecoili_raw_pos.txt"}; 
	//const char *argv[] = { "./smsAlignPositionLocate","-i","D:\\OTU拼接\\回国后文章\\数据\\E.coli_genome_NCBI.4641652.fa.npbss_simulated_CLR50X.fa", "-j","D:\\OTU拼接\\回国后文章\\数据\\E.coli_genome_NCBI.4641652.Normalize.fasta", "-T","5","-lib" ,"genomeLib.txt", "-o", "sequencePosition.txt" };
	//const char *argv[] = { "./smsAlign","-align","-i","D:\\OTU拼接\\回国后文章\\数据\\模拟数据\\E.coli_genome_NCBI.4641652.fa.npbss_simulated_CLR50X.fa", "-j","D:\\OTU拼接\\回国后文章\\数据\\模拟数据\\E.coli_genome_NCBI.4641652.Normalize.fasta", "-o", "aligned.txt","-pos","seqPos.txt"};
	bool flag = options.SetOptions(argc, argv);
	options.Validate();
	db_in = options.input;   //options member
	if (flag) {
			printf("Reading sequence file...");
			seq_db.Read(db_in.c_str(), options); // read the sequence file;
			printf("    Done.\n\n");
			seq_db.SequenceStatistic(options);
			//seq_db.buildSequenceSimilarityMatrix();
			seq_db.buildSequenceSimilaritySquareMatrix();
			// calculate the similarity
			printf("\nCalculating sequence identified...\n");
			seq_db.seq_align_global_my_multi_threads(options);
			// write to file
			printf("\nWriting sequence similarity to file...");
			seq_db.WriteAlignedSimilarityToFile(options);
			printf("    Done.\n\n");
		}
	else {
		bomb_error("Invalid option detected, exit!");
	}
	printf("Finished at: %s\n\n", getTime().c_str());
	//system("pause");
	return 0;
}
