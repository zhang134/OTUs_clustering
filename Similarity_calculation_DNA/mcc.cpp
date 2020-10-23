#include "mcc.h"
#include <assert.h>
#include <limits.h>
#include <algorithm>
#include <math.h>
#include <string>
#include <cstdio>
#include <sstream>


#ifndef NO_OPENMP

#include<omp.h>

#define WITH_OPENMP " (+OpenMP)"

#else

#define WITH_OPENMP ""

#define omp_set_num_threads(T) (T = T)
#define omp_get_thread_num() 0

#endif
int min(int x, int y) {
	return (x>y ? y : x);
}
int max(int x, int y) {
	return (x>y ? x : y);
}

//int aa2idx[] = { 0, 2, 4, 3, 6, 13,7, 8, 9,20,11,10,12, 2,20,14,
//5, 1,15,16,20,19,17,20,18, 6 };
// idx for  A  B  C  D  E  F  G  H  I  J  K  L  M  N  O  P
//          Q  R  S  T  U  V  W  X  Y  Z
// so  aa2idx[ X - 'A'] => idx_of_X, eg aa2idx['A' - 'A'] => 0,
// and aa2idx['M'-'A'] => 12

int aa2idx_ACGT[] = { 0, -1, 1, -1, -1, -1, 2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, 3 };
// index for          A   B  C   D   E   F  G   H   I   J   K   L   M   N   O   P   Q   R   S  T
// so  aa2idx_ACGT[ X - 'A'] => idx_of_X
// eg aa2idx_ACGT['A' - 'A'] => 0, aa2idx_ACGT['T' - 'A'] => 3. 

// kmer weight                            1  2   3   4   5     6      7      8      9      10
unsigned int NAAN_array[15] = { 1, 4, 16, 64, 256, 1024, 4096, 16384, 65536, 262144,
1048576, 4194304, 16777216, 67108864, 268435456 };
//                         11       12        13        14         15


// node for sort
struct node
{
	unsigned int data;
	unsigned int index_in_genome;
	unsigned int index_in_sequence;
};

int comp(const void *a, const void *b) { // ascending order for struct
	return (*(struct node *)a).data> (*(struct node *)b).data ? 1 : -1;
}

struct TempFile
{
	FILE *file;
	char buf[512];

	TempFile(const char *dir = NULL) {
		int len = dir ? strlen(dir) : 0;
		assert(len < 400);
		buf[0] = 0;
		if (len) {
			strcat(buf, dir);
			if (buf[len - 1] != '/' && buf[len - 1] != '\\') {
				buf[len] = '/';
				len += 1;
			}
		}
		strcat(buf, "cdhit.temp.");
		len += 11;
		sprintf(buf + len, "%p", this);
		file = fopen(buf, "w+");
	}
	~TempFile() {
		if (file) {
			fclose(file);
			remove(buf);
		}
	}
};

struct TempFiles
{
	NVector<TempFile*> files;

	~TempFiles() { Clear(); }

	void Clear() {
		int i;
#pragma omp critical
		{
			for (i = 0; i<files.size; i++) if (files[i]) delete files[i];
			files.Clear();
		}
	}
};

const char *temp_dir = "";
TempFiles temp_files;

FILE* OpenTempFile(const char *dir = NULL)
{
	TempFile *file = new TempFile(dir);
#pragma omp critical
	{
		temp_files.files.Append(file);
	}
	return file->file;
}
static void CleanUpTempFiles()
{
	if (temp_files.files.Size()) printf("Clean up temporary files ...\n");
	temp_files.Clear();
}

void bomb_error(const char *message)
{
	fprintf(stderr, "\nFatal Error:\n%s\nProgram halted !!\n\n", message);
	temp_files.Clear();
	exit(1);
} // END void bomb_error

string getTime() {
	time_t timep;
	time(&timep);
	char tmp[64];
	strftime(tmp, sizeof(tmp), "%Y/%d/%m %H:%M:%S", localtime(&timep));
	return tmp;
}


int print_usage() {
	cout << "*****************************************************************************" << endl;
	cout << "Calculate the MCC index for OTU clustering methods using globle alignmnet." << endl;
	cout << "=============== Usage example: ================" << endl;
	cout << "./mcc -i seq.fa -c cluster.txt" << endl << endl;
	cout << "-----required:" << endl;
	cout << "-i      input sequence fasta file (e.g. -i seq.fa)." << endl;
	cout << "-o      output file name (e.g. -o sim.txt)." << endl;
	//cout << "-c      cluster file by clustering methods (e.g. -c cluster.txt)." << endl;
	cout << " \n----optional:" << endl;
	cout << "-T      the threads number (default: your computer has)." << endl;
	cout << "-ms     match score in the sequence alignment (default: 2)." << endl;
	cout << "-ss     substitution (mismatch) score in the sequence alignment (default: -2)." << endl;
	cout << "-gs     gap score in the sequence alignment (default: -2)." << endl;
	cout << "-h      print this usage (or --help)." << endl;
	return 0;
}


bool Options::SetOptions(int argc, const char *argv[])
{
	int i, n;
	char date[100];
	strcpy(date, __DATE__);
	n = strlen(date);
	for (i = 1; i < n; i++) {
		if (date[i - 1] == ' ' && date[i] == ' ')
			date[i] = '0';
	}
	printf("================================================================\n");
	printf("Program: MCC, V" VERSION WITH_OPENMP ", %s, " __TIME__ "\n", date);
	printf("Your command:");
	//n = 9;
	for (i = 0; i<argc; i++) {
		printf(" %s", argv[i]);
	}
	printf("\n\n");
	time_t tm = time(NULL);
	printf("Started: %s", ctime(&tm));
	printf("================================================================\n");
	printf("                            Output                              \n");
	printf("----------------------------------------------------------------\n");
	for (i = 1; i + 1<argc; i += 2) {
		//printf("\n argv[%d]= %s", i, argv[i]);
		//printf("\n argv[%d]= %s", i + 1, argv[i + 1]);
		bool dddd = SetOption(argv[i], argv[i + 1]);
		if (dddd == 0) return false;
	}
	if (i < argc) return false;

	atexit(CleanUpTempFiles);
	return true;
}

bool Options::SetOption(const char *flag, const char *value)
{
	bool ffff = SetOptionCommon(flag, value);
	return ffff;
}

bool Options::SetOptionCommon(const char *flag, const char *value)
{
	int intval = atoi(value);
	if (strcmp(flag, "-i") == 0) input = value;
	else if (strcmp(flag, "-o") == 0) output = value;
	else if (strcmp(flag, "-M") == 0) max_memory = atol(value) * 1000000;
	//else if (strcmp(flag, "-b") == 0) banded_width = intval;
	else if (strcmp(flag, "-p") == 0) print = intval;
	//else if (strcmp(flag, "-tmp") == 0) temp_dir = value;
	else if (strcmp(flag, "-ms") == 0) match_score = intval;
	else if (strcmp(flag, "-ss") == 0) mismatch_score = intval;
	else if (strcmp(flag, "-gs") == 0) gap_score = intval;
	else if (strcmp(flag, "-h") == 0 || strcmp(flag, "--help") == 0) {
		print_usage();
		exit(1);
	}
	else if (strcmp(flag, "-T") == 0) {
#ifndef NO_OPENMP
		int cpu = omp_get_num_procs();
		threads = intval;
		if (threads > cpu) {
			threads = cpu;
			printf("Warning: total number of CPUs in your system is %i\n", cpu);
		}
		if (threads == 0) {
			threads = cpu;
			printf("total number of CPUs in your system is %i\n", cpu);
		}
		if (threads != intval) printf("Actual number of CPUs to be used: %i\n\n", threads);
#else
		printf("Option -T is ignored: multi-threading with OpenMP is NOT enabled!\n");
#endif
	}
	else return false;
	return true;
}

void Options::Validate()
{
#ifndef NO_OPENMP
	int cpu = omp_get_num_procs();
	threads = cpu;
#else
	printf("Option -T is ignored: multi-threading with OpenMP is NOT enabled!\n");
#endif
}

void Options::Print()
{
	//printf("query_file = %i\n", input);
	printf("threads = %i\n", threads);
	printf("match_score = %i\n", match_score);
	printf("mismatch_score = %i\n", mismatch_score);
	printf("gap_score = %i\n", gap_score);
	printf("print = %i\n", print);
}

void Sequence::Clear()
{
	if (data) delete[] data;
	/* do not set size to zero here, it is need for writing output */
	bufsize = 0;
	data = NULL;
}
Sequence::Sequence()
{
	memset(this, 0, sizeof(Sequence));
	nMatch = 0;
	nSubsitute = 0;
	nDelete = 0;
	nInsert = 0;
}

Sequence::Sequence(const Sequence & other)
{
	int i;
	//printf( "new: %p  %p\n", this, & other );
	memcpy(this, &other, sizeof(Sequence));
	if (other.data) {
		size = bufsize = other.size;
		data = new char[size + 1];
		//printf( "data: %p  %p\n", data, other.data );
		data[size] = 0;
		memcpy(data, other.data, size);
		//for (i=0; i<size; i++) data[i] = other.data[i];
	}
	if (other.identifier) {
		int len = strlen(other.identifier);
		identifier = new char[len + 1];
		memcpy(identifier, other.identifier, len);
		identifier[len] = 0;
	}
}
Sequence::~Sequence()
{
	//printf( "delete: %p\n", this );
	if (data) delete[] data;
	if (identifier) delete[] identifier;
}

void Sequence::operator=(const char *s)
{
	size = 0; // avoid copying;
	Resize(strlen(s));
	strcpy(data, s);
}
void Sequence::operator+=(const char *s)
{
	int i, m = size, n = strlen(s);
	Reserve(m + n);
	memcpy(data + m, s, n);
}
void Sequence::Resize(int n)
{
	int i, m = size < n ? size : n;
	size = n;
	if (size != bufsize) {
		char *old = data;
		bufsize = size;
		data = new char[bufsize + 1];
		if (data == NULL) bomb_error("Memory");
		if (old) {
			memcpy(data, old, m);
			delete[]old;
		}
		if (size) data[size] = 0;
	}
}
void Sequence::Reserve(int n)
{
	int i, m = size < n ? size : n;
	size = n;
	if (size > bufsize) {
		char *old = data;
		bufsize = size + size / 5 + 1;
		data = new char[bufsize + 1];
		if (data == NULL) bomb_error("Memory");
		if (old) {
			memcpy(data, old, m);
			delete[]old;
		}
	}
	if (size) data[size] = 0;
}

void Sequence::Swap(Sequence & other)
{
	Sequence tmp;
	memcpy(&tmp, this, sizeof(Sequence));
	memcpy(this, &other, sizeof(Sequence));
	memcpy(&other, &tmp, sizeof(Sequence));
	memset(&tmp, 0, sizeof(Sequence));
}
int Sequence::Format()
{
	int i, j = 0, m = 0;
	while (size && isspace(data[size - 1])) size--;
	if (size && data[size - 1] == '*') size--;
	if (size) data[size] = 0;
	for (i = 0; i<size; i++) {
		char ch = data[i];
		m += !(isalpha(ch) | isspace(ch));
	}
	if (m) return m;
	for (i = 0; i<size; i++) {
		char ch = data[i];
		if (isalpha(ch)) data[j++] = toupper(ch);
	}
	data[j] = 0;
	size = j;
	return 0;
}

void Sequence::SwapIn()
{
	if (data) return;
	if (swap == NULL) bomb_error("Can not swap in sequence");
	Resize(size);
	fseek(swap, offset, SEEK_SET);
	if (fread(data, 1, size, swap) == 0) bomb_error("Can not swap in sequence");
	data[size] = 0;
}
void Sequence::SwapOut()
{
	if (swap && data) {
		delete[] data;
		bufsize = 0;
		data = NULL;
	}
}

void SequenceDB::Read(const char *file, const Options & options)
{
	Sequence one;
	Sequence dummy;
	Sequence des;
	Sequence *last = NULL;
	FILE *swap = NULL;
	FILE *fin = fopen(file, "r");
	char *buffer = NULL;
	char *res = NULL;
	size_t swap_size = 0;
	if (fin == NULL) bomb_error("Failed to open the sequence file");
	Clear();
	dummy.swap = swap;
	buffer = new char[MAX_LINE_SIZE + 1];

	while (!feof(fin) || one.size) { /* do not break when the last sequence is not handled */
		buffer[0] = '>';
		if ((res = fgets(buffer, MAX_LINE_SIZE, fin)) == NULL && one.size == 0)
			break;
		if (buffer[0] == '+') {
			int len = strlen(buffer);
			int len2 = len;
			while (len2 && buffer[len2 - 1] != '\n') {
				if ((res = fgets(buffer, MAX_LINE_SIZE, fin)) == NULL) break;
				len2 = strlen(buffer);
				len += len2;
			}
			one.des_length2 = len;
			dummy.des_length2 = len;
			fseek(fin, one.size, SEEK_CUR);
		}
		else if (buffer[0] == '>' || buffer[0] == '@' || (res == NULL && one.size)) {
			if (one.size) { // write previous record
				one.dat_length = dummy.dat_length = one.size;
				if (one.identifier == NULL || one.Format()) {
					printf("Warning: from file \"%s\",\n", file);
					printf("Discarding invalid sequence or sequence without header and description!\n\n");
					if (one.identifier) printf("%s\n", one.identifier);
					printf("%s\n", one.data);
					one.size = 0;
				}
				one.index = dummy.index = sequences.size();
				if (one.size > 0) {
					if (swap) {
						swap_size += one.size;
						// so that size of file < MAX_BIN_SWAP about 2GB
						if (swap_size >= MAX_BIN_SWAP) {
							dummy.swap = swap = OpenTempFile(temp_dir);
							swap_size = one.size;
						}
						dummy.size = one.size;
						dummy.offset = ftell(swap);
						dummy.des_length = one.des_length;
						sequences.Append(new Sequence(dummy));
						//one.ConvertBases();
						fwrite(one.data, 1, one.size, swap);
					}
					else {
						//printf( "==================\n" );
						sequences.Append(new Sequence(one));
						//printf( "------------------\n" );
						//if( sequences.size() > 10 ) break;
					}
					//if( sequences.size() >= 10000 ) break;
				}
			}
			one.size = 0;
			one.des_length2 = 0;

			int len = strlen(buffer);
			int len2 = len;
			des.size = 0;
			des += buffer;
			while (len2 && buffer[len2 - 1] != '\n') {
				if ((res = fgets(buffer, MAX_LINE_SIZE, fin)) == NULL) break;
				des += buffer;
				len2 = strlen(buffer);
				len += len2;
			}
			size_t offset = ftell(fin);
			one.des_begin = dummy.des_begin = offset - len;
			one.des_length = dummy.des_length = len;

			int i = 0;
			if (des.data[i] == '>' || des.data[i] == '@' || des.data[i] == '+')
				i += 1;
			if (des.data[i] == ' ' || des.data[i] == '\t')
				i += 1;
			//if (options.des_len && options.des_len < des.size) des.size = options.des_len;
			while (i < des.size && !isspace(des.data[i]))
				i += 1;
			des.data[i] = 0;
			one.identifier = dummy.identifier = des.data;
		}
		else {
			one += buffer;
		}
	}
#if 0
	int i, n = 0;
	for (i = 0; i<sequences.size(); i++) n += sequences[i].bufsize + 4;
	cout << n << "\t" << sequences.capacity() * sizeof(Sequence) << endl;
	int i;
	scanf("%i", &i);
#endif
	one.identifier = dummy.identifier = NULL;
	delete[] buffer;
	fclose(fin);
}
void SequenceDB::SequenceStatistic(Options & options)
{
	int i, j, k, len;
	int N = sequences.size();
	total_letter = 0; // total bases
	total_desc = 0;   // total letters of sequence headers
	max_len = 0;
	mean_len = 0.0;
	min_len = (size_t)-1;
	for (i = 0; i<N; i++) {
		Sequence *seq = sequences[i];
		len = seq->size;
		mean_len += seq->size;
		total_letter += len;
		if (len > max_len)
			max_len = len;
		if (len < min_len)
			min_len = len;
		if (seq->swap == NULL)
			//			seq->ConvertBases();
			if (seq->identifier)
				total_desc += strlen(seq->identifier);
	}
	mean_len = mean_len / N;
	cout << "  Total number: " << N << endl;
	cout << "  Longest:      " << max_len << endl;
	cout << "  Shortest:     " << min_len << endl;
	cout << "  Mean length:  " << mean_len << endl;
	cout << "  Total bases:  " << total_letter << endl;
	cout << endl;
	// END change all the NR_seq to iseq

	//len_n50 = (max_len + min_len) / 2; // will be properly set, if sort is true;

}// END sort_seqs_divide_segs

void SequenceDB::buildSequenceSimilarityMatrix() {
	// build sequence similarity matrix: 3 columns
	// seq0 seq1 0.98
	int N = sequences.size();
	int row_num = (N  * (N - 1) )/ 2;
	int memory_matrix = row_num * 3 * 4 / (1024 * 1024);
	cout << "Similarity matrix needs memory: ~" << memory_matrix << "Mb" << endl;
	simlarity_matrix = new float*[row_num + 1];
	for (int i = 0; i < row_num; i++) {
		simlarity_matrix[i] = new float[3] ();
	}
}

void SequenceDB::buildSequenceSimilaritySquareMatrix() {
	int N = sequences.size();
	simlarity_matrix = new float*[N];
	for (int i = 0; i < N; i++) {
                simlarity_matrix[i] = new float[N] ();
        }

}

void SequenceDB::WriteAlignedSimilarityToFile(Options & options) {
	ofstream outfile2;
	const char* filename = options.output.data();
	outfile2.open(filename, ios::out);
	//outfile2.open(options.output, ios::app);
	if (!outfile2.is_open())
		bomb_error("Open output file failure");
	int N = sequences.size();
	int row_num = (N  * (N - 1)) / 2;
	//outfile2 << "SeqID\tSeqID\tIdentify" << endl;
	for (int i = 0; i < N; i++) {
		for (int j = 0; j < N; j++){
			 outfile2 << simlarity_matrix[i][j] << "\t";
		}
		outfile2 << endl;
		//if (simlarity_matrix[i][2] < 0.90)
		//	continue;
		//outfile2 << (int)simlarity_matrix[i][0] << "\t" << (int)simlarity_matrix[i][1] << "\t" << simlarity_matrix[i][2] << endl;
	}
	outfile2.close();
}

void SequenceDB::seq_align_global_my_multi_threads(const Options & options) {
	int tid = 0;
	float p = 0.0, p0 = 0.0, sim;
	int seq_num = sequences.size();
	int *row_num = new int[seq_num]();
	row_num[0] = 0;
	for (int ii = 1; ii < seq_num; ii++) {
		row_num[ii] = row_num[ii - 1] + seq_num - ii;
	}
	for (int i = 0; i < seq_num; i++) {
#pragma omp parallel for schedule( dynamic, 1 )
		for (int j = i + 1; j < seq_num; j++) {
			sim = seq_align_global_my(sequences[i], sequences[j], i, j, options); 
			//int index_row = row_num[i] + j - i - 1;
			//simlarity_matrix[i][j] = sim;//(float)i;
			//simlarity_matrix[j][i] = sim;//(float)j;
			//simlarity_matrix[index_row][2] = sim;
			//printf("matrix[%f][%f] = %f,index = %d\n", simlarity_matrix[index_row][0], simlarity_matrix[index_row][1], simlarity_matrix[index_row][2], index_row);
		}
		 //printf("%d-th in sequence: %d, in genome: %d\n", i, j);																//printf("The %d-th sequence position is located\n", i);
		tid = omp_get_thread_num();
		if (omp_get_thread_num() == 0) {
			p = (100.0*i) / sequences.size();
			if (p > p0 + 1E-1) {
				printf("\r%5.1f%%   %8d-th sequence", p, i); //printf("The %d-th sequence is located", i);
				p0 = p;
			}
		}
		fflush(stdout);
	}
}


float SequenceDB::seq_align_global_my(Sequence *seq1, Sequence *seq2, int id1, int id2, const Options & options) {
	int len1 = seq1->size;
	int len2 = seq2->size;
	int gapScore = options.gap_score;
	int matchScore = options.match_score;
	int mismatchScore = options.mismatch_score;
	int **matrix;
	int **track;
	float simm;
	string seq1_aligned = "";
	string seq2_aligned = "";
	string middle = "";
	int from_up, from_diag, from_left;
	int match_num = 0;
	//  up 1, diag 2, left 3
	matrix = new int*[len1 + 1]; //设置行 或直接double **data=new double*[m]; 一个指针指向一个指针数组。  
	track = new int*[len1 + 1];
	for (int j = 0; j <= len1; j++) {
		matrix[j] = new int[len2 + 1] ();        //这个指针数组的每个指针元素又指向一个数组
		track[j] = new int[len2 + 1] ();
	}
	// initialize the first column
	for (int i = 0; i <= len1; i++) {
		matrix[i][0] = i * gapScore;
		track[i][0] = 1;
	}
	// initialize the first row
	for (int j = 0; j <= len2; j++) {
		matrix[0][j] = j * gapScore;
		track[0][j] = 3;
	}

	for (int i = 1; i <= len1; i++) {
		for (int j = 1; j <= len2; j++) {
			if (seq1->data[i - 1] == seq2->data[j - 1])
				from_diag = matrix[i - 1][j - 1] + matchScore;
			else
				from_diag = matrix[i - 1][j - 1] + mismatchScore;
			from_up = matrix[i - 1][j] + gapScore;
			from_left = matrix[i][j - 1] + gapScore;
			if (from_diag >= from_up && from_diag >= from_left) {
				matrix[i][j] = from_diag;
				track[i][j] = 2;//  up 1, diag 2, left 3
			}

			else if (from_up >= from_diag && from_up >= from_left) {
				matrix[i][j] = from_up;
				track[i][j] = 1;//  up 1, diag 2, left 3
			}
			else {
				matrix[i][j] = from_left;
				track[i][j] = 3;//  up 1, diag 2, left 3
			}

		}
	}

	//string haha = "AXNXCXXGXXXXXXXXT"; // A 0, C 4, G 7, T 16, N 2
	int i = len1;
	int j = len2;
	for (; i > 0 && j > 0;) {
		if (track[i][j] == 2) {//  up 1, diag 2, left 3
			seq1_aligned += seq1->data[i - 1];
			seq2_aligned += seq1->data[j - 1];
			if (seq1->data[i - 1] == seq2->data[j - 1]) {
				middle += "|";
				match_num++;
			}
			else
				middle += " ";
			i--;
			j--;
		}
		else if (track[i][j] == 1) {//  up 1, diag 2, left 3
			seq1_aligned += seq1->data[i - 1];
			i--;
			seq2_aligned += "-";//haha.at((int)seq2->data[j]);
			middle += " ";
		}
		else {
			seq1_aligned += "-";//haha.at((int)seq1->data[i - 1]);
			seq2_aligned += seq2->data[j - 1];
			j--;
			middle += " ";
		}
	}
	if (seq1_aligned.size() != seq2_aligned.size()) {
		bomb_error("Sequence aligned error: two aligned sequences lengths are not equal!");
	}
	reverse(seq1_aligned.begin(), seq1_aligned.end());
	reverse(seq2_aligned.begin(), seq2_aligned.end());
	reverse(middle.begin(), middle.end());
#if 0
	cout << seq1_aligned << endl;
	cout << middle << endl;
	cout << seq2_aligned << endl;
#endif

	// release memory
	for (int i = 0; i <= len1; i++) {
		delete[] matrix[i]; //先撤销指针元素所指向的数组
		delete[] track[i];
	}
	delete[] matrix;
	delete[] track;

	simm =  match_num / (float)seq1_aligned.size();
	simlarity_matrix[id1][id2] = simm;
	simlarity_matrix[id2][id1] = simm;
	return simm;
}
