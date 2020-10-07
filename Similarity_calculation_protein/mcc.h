#pragma once
#include<iostream>
#include<fstream>
#include<iomanip>
#include<cstdlib>
#include<stdio.h>
#include<string.h>
#include<ctype.h>
#include<stdint.h>
#include<time.h>
#include<valarray>
#include<vector>
#include<map>

#define VERSION  "1.0"
#define MAX_LINE_SIZE 300000
#define MAX_BIN_SWAP 2E9

using namespace std;

template<class TYPE>
class Vector : public vector<TYPE>
{
public:
	Vector() : vector<TYPE>() {}
	Vector(size_t size) : vector<TYPE>(size) {}
	Vector(size_t size, const TYPE & deft) : vector<TYPE>(size, deft) {}

	void Append(const TYPE & item) {
		size_t n = this->size();
		if (n + 1 >= this->capacity()) this->reserve(n + n / 5 + 1);
		this->push_back(item);
	}
	int size()const { return (int)vector<TYPE>::size(); }
};

template<class TYPE>
class NVector
{
public:
	TYPE * items;
	int     size;
	int     capacity;

	NVector() { size = capacity = 0; items = NULL; }
	NVector(int n, const TYPE & v = TYPE()) {
		size = capacity = 0; items = NULL;
		Resize(n, v);
	}
	NVector(const NVector & other) {
		size = capacity = 0; items = NULL;
		if (other.items) {
			Resize(other.size);
			memcpy(items, other.items, other.size * sizeof(TYPE));
		}
	}

	~NVector() { if (items) free(items); }

	int  Size()const { return size; }
	void Clear() {
		if (items)
			free(items);
		size = capacity = 0; items = NULL;
	}

	void Resize(int n, const TYPE & value = TYPE()) {
		if (n == size && capacity > 0) return;
		int i;
		// When resize() is called, probably this is the intended size,
		// and will not be changed frequently.
		if (n != capacity) {
			capacity = n;
			items = (TYPE*)realloc(items, capacity * sizeof(TYPE));
		}
		for (i = size; i<n; i++)
			items[i] = value;
		size = n;
	}
	void Append(const TYPE & item) {
		if (size + 1 >= capacity) {
			capacity = size + size / 5 + 1;
			items = (TYPE*)realloc(items, capacity * sizeof(TYPE));
		}
		items[size] = item;
		size++;
	}

	TYPE& operator[](const int i) {
		//if( i <0 or i >= size ) printf( "out of range\n" );
		return items[i];
	}
	TYPE& operator[](const int i)const {
		//if( i <0 or i >= size ) printf( "out of range\n" );
		return items[i];
	}
};
typedef NVector<int>      VectorInt;
typedef Vector<VectorInt> MatrixInt;

void bomb_error(const char *message);
int print_usage();
string getTime();

struct Options
{
	size_t  max_memory; // -M: 400,000,000 in bytes
	int     print;
	int     threads;
	int		match_score;
	int		mismatch_score;//substitution
	int		gap_score;
	size_t  max_entries;
	size_t  max_sequences;
	size_t  mem_limit;

	string  input;
	string  output;

	Options() {
		max_memory = 800000000;
		print = 0;
		threads = 1;
		max_entries = 0;
		max_sequences = 1 << 20;
		mem_limit = 100000000;
		match_score = 2;
		mismatch_score = -2;
		gap_score = -2;
	};

	bool SetOptionCommon(const char *flag, const char *value);
	bool SetOption(const char *flag, const char *value);
	bool SetOptions(int argc, const char *argv[]);
	void Validate();
	void Print();
};

struct Sequence {
	// real sequence, if it is not stored swap file:
	char *data;
	int   size; // length of the sequence
	int   bufsize;
	double  identify; // the identify with genome
	FILE *swap;
	// stream offset of the sequence:
	int   offset;
	// stream offset of the description string in the database:
	size_t   des_begin;
	// length of the description:
	int   des_length;
	// length of the description in quality score part:
	int   des_length2;
	// length of data in fasta file, including line wrapping:
	int   dat_length;
	char *identifier; // sequence header
	int id_sorted;// sorted index of the sequence
	int   index;  // index of the sequence in the original database:
	short state;
	int   coverage[4];
	int nMatch;
	int nSubsitute;
	int nDelete;
	int nInsert;
	int nAlignedBase;
	Sequence();
	Sequence(const Sequence & other);
	~Sequence();
	void Clear();
	void operator=(const char *s);
	void operator+=(const char *s);
	void Resize(int n);
	void Reserve(int n);
	void Swap(Sequence & other);
	int Format();
	//void ConvertBases();
	void SwapIn();
	void SwapOut();
};

class SequenceDB
{
public:

	Vector<Sequence*>  sequences;
	long long total_letter;
	long long total_desc;
	size_t max_len;
	float mean_len;
	size_t min_len;
	float **simlarity_matrix;
	void Clear() {
		for (int i = 0; i < sequences.size(); i++) delete sequences[i];
		sequences.clear(); //rep_seqs.clear();
	}

	SequenceDB() {
		total_letter = 0;
		total_desc = 0;
		min_len = 0;
		max_len = 0;
		//len_n50 = 0;
		mean_len = 0.0;
	}
	~SequenceDB() { Clear(); }

	void Read(const char *file, const Options & options);
	void SequenceStatistic(Options & options);
	void buildSequenceSimilarityMatrix();
	void buildSequenceSimilaritySquareMatrix();
	void WriteAlignedSimilarityToFile(Options & options);
	float seq_align_global_my(Sequence *seq1, Sequence *seq2, int id1, int id2, const Options & options);
	float seq_align_protein_semi_global(Sequence *seq1, Sequence *seq2, int id1, int id2, const Options & options);
	void seq_align_global_my_multi_threads(const Options & options);

};
