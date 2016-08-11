#ifndef LARGEVIS_H
#define LARGEVIS_H

#include "largeVis.h"

#include <stdlib.h>
#include <math.h>
#include <vector>

#include <pthread.h>

typedef float visreal;

struct arg_struct{
	void *ptr;
	int id;
	arg_struct(void *x, int y) :ptr(x), id(y){}
};

class LargeVis{
private:
	arma::sp_mat*			edges;
	long long N, D, n_samples, n_threads, M, edge_count_actual;
	visreal initial_rho, gamma, perplexity;
	visreal *vis;
	std::vector<int> *knn_vec;
	long long E, *head;
	std::vector<long long> next, reverse;
	std::vector<int> edge_from, edge_to;
	std::vector<visreal> edge_weight;
	int *neg_table;
	long long neg_size;
	long long *alias;
	visreal *prob;
	std::uniform_real_distribution<double> rnd;
	std::mt19937_64 mt;
	void init_alias_table();
	long long sample_aE(visreal rand_value1, visreal rand_value2);

	void compute_similarity_thread(int id);
	static void *compute_similarity_thread_caller(void *arg);
	void search_reverse_thread(int id);
	static void *search_reverse_thread_caller(void *arg);
	void init_neg_table();
	void visualize_thread(int id);
	static void *visualize_thread_caller(void *arg);
public:
	LargeVis(int D, long long N, long long E, arma::sp_mat* edges);
	void compute_similarity(visreal perp);
	arma::fmat visualize(visreal gamm, long long n_samp, long long n_neg);


};

#endif