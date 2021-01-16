
// [[Rcpp::plugins(cpp11)]]
#include <RcppArmadillo.h>
#include <stack>
#include <iostream>
#include <cmath>
#include <vector>
#include <stack>
#include <algorithm>

// [[Rcpp::depends(RcppArmadillo)]]

using namespace std;



// declare functions to use later
/* Equal probability sampling; without-replacement case */
// pulled from R-2.0.1/src/main/random.c lines 352-364
// slightly modified by KQ 2/21/2005
// x: n array of original indices running from 0 to (n-1)
// y: k array of samled indices
// k: length of y (must be <= n)
// n: length of x (must be >= k)
// revise this function to work with vectors
static void SampleNoReplace(const int k, int n, vector<int>& y, vector<int>& x) {
	for (int i = 0; i < n; i++)
		x[i] = i;
	for (int i = 0; i < k; i++) {
		double u = R::runif(0, 1);//runif(1,0,1) has a NumericVecor output, so I need to index to get the double number out
		int j = static_cast<int>(n * u);
		y[i] = x[j];
		x[j] = x[--n];
	}
}




// draws 1 element from a double array x with prob weights in array p
// n is length of x and p 
int ProbSamp(vector<int>& x, double* p, const int& n) {
	const double u = R::runif(0, 1);
	double cumprob = p[0];
	for (int i = 0; i < (n - 1); ++i) {
		if (u <= cumprob) {
			return x[i];
		}
		else {
			cumprob += p[i + 1];
		}
	}
	return x[n - 1];
}
// draws 1 element from a double array x with prob weights in array p
// n is length of x and p 
double ProbSamp(vector<double>& x, double* p, const int& n) {
	const double u = R::runif(0, 1);
	double cumprob = p[0];
	for (int i = 0; i < (n - 1); ++i) {
		if (u <= cumprob) {
			return x[i];
		}
		else {
			cumprob += p[i + 1];
		}
	}
	return x[n - 1];
}


// beta-binomial density at y, sample size, n, and parameters a, and b
double dbetabinom(const int& y, const int& n,
	const double& a, const double& b,
	const bool log = false) {
	/*
	cout << "in dbetabinom" << endl;
	cout << "y = " << y << endl;
	cout << "n = " << n  << endl;
	cout << "a = " << a << endl;
	cout << "b = " << b  << endl << endl;
	*/
	const double term1 = std::lgamma(n + 1.0) - std::lgamma(y + 1.0) -
		std::lgamma(n - y + 1.0);
	const double term2 = std::lgamma(a + y) + std::lgamma(n + b - y) - std::lgamma(a + b + n);
	const double term3 = std::lgamma(a + b) - std::lgamma(a) - std::lgamma(b);
	if (log) {
		return(term1 + term2 + term3);
	}
	else {
		return(std::exp(term1 + term2 + term3));
	}
}


// product beta binomial density 
// takes arrays of a and b parameters so most appropriate 
// for calculating the probability of a new cluster 
// assumes the binomial sample size is always 1
// y: J-array of 0-1 data
// a: J-array of success counts
// b: J-array of failure counts
// J: length of arrays

// log scale
double prodbetabinomG0(const int* y, const double* a, const double* b,
	const int& J) {
	double funval = 0.0;
	for (int j = 0; j < J; ++j) {
		if (y[j] == 0 || y[j] == 1) {
			/*
			cout << "in G0" << endl;
			cout << "y[j] = " << y[j] << endl;
			cout << "a[j] = " << a[j] << endl;
			cout << "b[j] = " << b[j] << endl;
			*/
			const double holder = dbetabinom(y[j], 1, a[j], b[j], true);
			funval += holder;
		}
	}
	return funval;
}


// product beta binomial density 
// takes arrays of a and b as well as additional counts as parameters
// so most appropriate 
// for calculating the probability of y being in a particular existing 
// cluster 
// assumes the binomial sample size is always 1
// y: J-array of 0-1 data
// a: J-array of success counts
// b: J-array of failure counts
// n1: J-vector of success counts
// n0: J-vector of failure counts
// J: length of arrays

//log scale
double prodbetabinomH(const int* y, const double* a, const double* b,
	const int* n1,
	const int* n0,
	const int& J) {

	double funval = 0.0;
	for (int j = 0; j < J; ++j) {

		if (y[j] == 0 || y[j] == 1) {
			/*
			cout << "in H" << endl;
			cout << "y[j] = " << y[j] << endl;
			cout << "a[j] = " << a[j] << endl;
			cout << "b[j] = " << b[j] << endl;
			cout << "n1[j] = " << n1[j] << endl;
			cout << "n0[j] = " << n0[j] << endl << endl;
			*/
			const double holder = dbetabinom(y[j], 1, a[j] + n1[j], b[j] + n0[j], true);
			funval += holder;
		}
	}
	return funval;
}

// product beta binomial density 
// takes arrays of a and b as well as additional counts as parameters
// this function differs from prodbetabinomH in that it subtracts the
// indicators in y[i] from n1 and n0 as appropriate and in this sense
// it is most useful for calculating the probability of y being in its
// current cluster
// assumes the binomial sample size is always 1
// y: J-array of 0-1 data
// a: J-array of success counts
// b: J-array of failure counts
// n1: J-vector of success counts
// n0: J-vector of failure counts
// J: length of arrays
// t: the time priod
//log scale
double prodbetabinomH_not_i(const int* y, const double* a, const double* b,
	const int* n1,
	const int* n0,
	const int& J) {

	int n1_not_i[J];
	int n0_not_i[J];

	for (int j = 0; j < J; ++j) {
		if (y[j] == 0) {
			n0_not_i[j] = n0[j] - 1;
			n1_not_i[j] = n1[j];
		}
		else if (y[j] == 1) {
			n1_not_i[j] = n1[j] - 1;
			n0_not_i[j] = n0[j];
		}
	}

	double funval = 0.0;
	for (int j = 0; j < J; ++j) {
		if (y[j] == 0 || y[j] == 1) {
			/*
			cout << "in H not i" << endl;
			cout << "y[j] = " << y[j] << endl;
			cout << "a[j] = " << a[j] << endl;
			cout << "b[j] = " << b[j] << endl;
			cout << "n1[j] = " << n1[j] << endl;
			cout << "n0[j] = " << n0[j]  << endl;
			cout << "n0_not_i[j] = " << n0_not_i[j] << endl;
			cout << "n1_not_i[j] = " << n1_not_i[j] << endl << endl;
			*/
			const double holder = dbetabinom(y[j], 1, a[j] + n1_not_i[j],
				b[j] + n0_not_i[j], true);
			funval += holder;
		}
	}
	return funval;
}

double alpha_logLik(const double& alpha, const int& n, const int& k){
	return k*std::log(alpha)+std::lgamma(alpha)-std::lgamma(alpha+n);
}

//created on Jan 20, 2020
//This is based on a new model. I put a log normal prior directly on alpha.
//Then log(\alpha) has a normal distribution in every period
//I use dyanmic linear model across period to smooth it.  


// set seed: this guarantees reproducibiltiy, and it's called in MCMCddpmbb() function

void set_seed(double seed) {
    Rcpp::Environment base_env("package:base");
    Rcpp::Function set_seed_r = base_env["set.seed"];
    set_seed_r(std::floor(std::fabs(seed)));
}



//' Dynamic Dirichlet Process Mixture Beta-Binomial Model to Indentify Voting Coalitions in Roll Call Vote Data of Multiple Periods.
//'
//' @export
// [[Rcpp::export]]
Rcpp::List MCMCddpmbb(arma::cube X, arma::Mat<int> dataStructure,  
	double lambda1_value, double lambda0_value, arma::Col<double> alpha_start, 
	arma::Col<double> tau, double W = 1, double m0 = 0, double H0 = 1, double r0 = 1, 
	double s0 = 1, double V_start = 1, int sams_iter=3, int mcmc=1000, 
	int burnin=0, int thin=1, int verbose = 100, int seed=1234, int fix_W = 1, double tau_W = 0.05, double r1 = 0.1, double s1 = 1) {
	
	//change argument names to fit the function body
	//change tau to tune_value
	arma::Col<double> tune_value(tau.n_rows);
	for(unsigned int i=0; i<tune_value.n_rows; ++i){
		tune_value(i)=tau(i);
	}
	
	//change  H0 to C0
	double C0=H0;
	
	
	//set random seed
	set_seed(seed);   
	// MCMC-related quantities
	const int tot_iter = burnin + mcmc;
	const int nsamp = mcmc / thin;
	// the number of periods
	const int T = dataStructure.n_rows;
	int I[T]; // number of individuals for each period
	for (int i = 0; i < T; ++i) {
		I[i] = dataStructure(i, 0);
	}
	int J[T]; // number of items for each period
	for (int i = 0; i < T; ++i) {
		J[i] = dataStructure(i, 1);
	}

	int max_I = *std::max_element(I, I + T);//I can directly do pointer mathematics on the array
	int max_J = *std::max_element(J, J + T);

	int y[T][max_I][max_J];    // the observed binary responses: 3-D array 
	// (0, 1, anything else is missing)
	for (int t = 0; t < T; ++t) {
		for (int i = 0; i < max_I; ++i) {
			for (int j = 0; j < max_J; ++j) {
				y[t][i][j] = X(j, i, t);
				//arma::cube index-rule: the 3rd position index a slice matrix, 1st and 2nd positions index matrix postion
				//indexing starts with 0
			}
		}
	}





	// starting values
	double alpha[T];
	double alpha_log[T];
	for (int t = 0; t < T; ++t) {
		alpha[t] = alpha_start(t);
		alpha_log[t] = std::log(alpha[t]);
	}

	
	//MH_accept monitors the acceptance of new z[t] for each Metropolis-Hastings step for each time period
	int MH_accept[T];
	for (int t = 0; t < T; ++t) {
		MH_accept[t] = 0;
	}
	
	
	//tune parameter for Metroplis-Hastings step
	double tune[T];
	for (int t = 0; t < T; ++t) {
		tune[t] = tune_value(t);
	}

	//the DLM mean parameter gamma[t], the length should be T
	double gamma[T];
	for (int t = 0; t < T; ++t) {
		gamma[t] = 0;
	}
	//The rest of the DLM paramters containers
	double a[T];
	double R[T];
	double f[T];
	double Q[T];
	double e[T];
	double A[T];
	double m[T];
	double C[T];
	for (int t = 0; t < T; ++t) {
		a[t]=0;
		R[t]=0;
		f[t]=0;
		Q[t]=0;
		e[t]=0;
		A[t]=0;
		m[t]=0;
		C[t]=0;
	}

	// the DLM variance for the innovation equation
	double V = V_start;
	//store containers for V, gamma
	arma::Col<double> V_store_data(nsamp);
	arma::mat gamma_store_data(nsamp, T);
	gamma_store_data.fill(0.0);
	
	//store containers for W
	arma::Col<double> W_store_data(nsamp);
	
	
	//MH_accept_W monitors the acceptance of new W for each Metropolis-Hastings step for each time period
	int MH_accept_W = 0;
	
	
	double lambda1[T][max_J]; // pseudo counts of y=1 for beta base distribution
	double lambda0[T][max_J]; // pseudo counts of y=0 for beta base distribution
	
	// For now, I don't sample lambda
	for (int t = 0; t < T; ++t) {
		for (int j = 0; j < J[t]; ++j) {
			lambda1[t][j] = lambda1_value;
			lambda0[t][j] = lambda0_value;
		}
	}
	

	
	arma::mat alpha_store_data(nsamp,T);
	alpha_store_data.fill(0.0);
	

	// storage arrays
	arma::Cube<int> c_store(nsamp, max_I,T);
	c_store.fill(-100);
	

	arma::Mat<int> nclust_store_data(nsamp, T);
	nclust_store_data.fill(0);

	
	



	// initial values for cluster variables and sufficient statistics
	vector<vector<int>> cr;  // gives the row index of the cluster each 
	for (int t = 0; t < T; ++t) {
		vector<int> holder;
		holder.reserve(I[t] + 1); // in period t, obs belongs to 
		for (int i = 0; i < I[t]; ++i) {
			// initially put everyone in a separate cluster for each time priod
			holder.push_back(i);
		}
		cr.push_back(holder);
	}
	
	// cru holds the unique cluster rows
	vector<vector<int>> cru = cr;
	for (int i = 0; i < (int)cru.size(); ++i) {
		sort(cru[i].begin(), cru[i].end());
		cru[i].erase(unique(cru[i].begin(), cru[i].end()),
			cru[i].end());
	}
	

	// K[t] is the number of distinct clusters for each period t before any sampling
	int K[T];
	for (int t=0; t < T; ++t) {
		K[t] = cru[t].size();
	}
	// stack of available cluster rows 
	vector<stack<int>> castack;
	for (int t = 0; t < T; ++t) {
		stack<int> holder;
		for (int k = K[t]; k < (I[t] + 1); ++k) {
			holder.push(k);
		}
		castack.push_back(holder);
	}
	/*for (int k = K; k < (I + 1); ++k) {
		castack.push(k);
	}*/


	int n_0[T][max_I + 1][max_J]; // holds number of y_kj == 0 for cluster c(2nd position) for item j (3rd position) in period t
	int n_1[T][max_I + 1][max_J]; // holds number of y_kj == 1 for cluster c(2nd position) for item j (3rd position) in period t
	// fill with 0s
	for (int t = 0; t < T; ++t) {
		for (int i = 0; i < (max_I + 1); ++i) {
			for (int j = 0; j < max_J; ++j) {
				n_0[t][i][j] = 0;
				n_1[t][i][j] = 0;
			}
		}
	}
	
	// fill with counts based on y
	for (int t = 0; t < T; ++t) {
		for (int kk = 0; kk < K[t]; ++kk) {
			int k = cru[t][kk];
			for (int j = 0; j < J[t]; ++j) {
				for (int i = 0; i < I[t]; ++i) {
					if (cr[t][i] == k) {
						if (y[t][i][j] == 0) {
							++n_0[t][k][j];
						}
						if (y[t][i][j] == 1) {
							++n_1[t][k][j];
						}
					}
				}
			}
		}
	}
	




	// number of elements in each cluster 
	int n_in_cluster[T][max_I + 1];
	// fill with 0s
	for (int t = 0; t < T; ++t) {
		for (int i = 0; i < (max_I + 1); ++i) {
			n_in_cluster[t][i] = 0;
		}
	}
	
	// put correct data in
	for (int t = 0; t < T; ++t) {
		for (int kk = 0; kk < K[t]; ++kk) {
			int k = cru[t][kk];
			for (int i = 0; i < I[t]; ++i) {
				if (cr[t][i] == k) {
					++n_in_cluster[t][k];
				}
			}
		}
	}
	

	//I may  need to initilize them over and over again
	//arrays for permutations
	vector<vector<int>> I_array;
	vector<vector<int>> I_inds_perm;
	vector<vector<int>> Iseq_array;//I use Iseq_array to assign new cluster label to an individual. Therefore, the one more length (position) is for a potential new cluster
	for (int t = 0; t < T; ++t) {
		vector<int> holder_I_array;
		vector<int> holder_I_inds_perm;
		for (int i = 0; i < I[t]; ++i) {
			holder_I_array.push_back(0);
			holder_I_inds_perm.push_back(0);
		}
		I_array.push_back(holder_I_array);
		I_inds_perm.push_back(holder_I_inds_perm);
	}

	for (int t = 0; t < T; ++t) {
		vector<int> holder_Iseq_array;
		for (int i = 0; i < I[t] + 1; ++i) {
			holder_Iseq_array.push_back(i); //for each period Iseq_array is one unit longer than I_array and I_inds_perm
		}
		Iseq_array.push_back(holder_Iseq_array);
	}

	/*
	//original code for this part
	int I_inds_perm[I];
	int Iseq_array[I + 1];
	for (int i = 0; i < (I + 1); ++i) {
		Iseq_array[i] = i;
	}*/


	//sampling starts

	arma::Col<int> merge_split_accepts(T);
	merge_split_accepts.fill(0);
	int merge_split_count = 0; // merge_split is always done for all periods at the same time, so one counter is good
	int count = 0; //use this count to store draws, one count is one sample index
	
	arma::Col<int> n_split(T);
	n_split.fill(0);
	arma::Col<int> n_merge(T);
	n_merge.fill(0);



	for (int iter = 0; iter < tot_iter; ++iter) {
		// alpha log scale for new cluster
		
		if (iter % sams_iter == 0) {
			++merge_split_count;
			for (int t = 0; t < T; ++t) {
				// randomly select indices i1 and i2 from {0,1,...,I[t]}
				int i1 = static_cast<int>(std::floor(R::runif(0, 1) * I[t]));
				int i2 = static_cast<int>(std::floor(R::runif(0, 1) * I[t]));
				while (i1 == i2) {
					i2 = static_cast<int>(std::floor(R::runif(0, 1) * I[t]));
				}

				// if cr[t][i1] == cr[t][i2] attempt a split
				if (cr[t][i1] == cr[t][i2]) {
					// find all obs that are in the same cluster as i1 and i2
					// remove i1 and i2 from this set of obs and call this set S
					vector<int> Svec;
					Svec.reserve(I[t]);
					int n0_S[J[t]];
					int n1_S[J[t]];
					for (int j = 0; j < J[t]; ++j) {
						n0_S[j] = 0;
						n1_S[j] = 0;
					}
					for (int i = 0; i < I[t]; ++i) {
						if (cr[t][i] == cr[t][i1]) {
							if (i != i1 && i != i2) {
								Svec.push_back(i);
							}
						}
					}
					int NS = Svec.size(); // number of elements in S
					/*
					cout << "Svec = " << endl;
					printVector(Svec);
					*/
					// form initial singleton sets S_i1 and S_i2
					// S_i1 quantities, now only record i1's info in S_i1, will update more later
					vector<int> S_i1;
					S_i1.reserve(I[t]);
					S_i1.push_back(i1);
					int N_Si1 = S_i1.size();
					int n0_Si1[J[t]];
					int n1_Si1[J[t]];
					for (int j = 0; j < J[t]; ++j) {
						if (y[t][i1][j] == 0) {
							n0_Si1[j] = 1;
							n1_Si1[j] = 0;
							++n0_S[j];
						}
						else if (y[t][i1][j] == 1) {
							n0_Si1[j] = 0;
							n1_Si1[j] = 1;
							++n1_S[j];
						}
						else {
							n0_Si1[j] = 0;
							n1_Si1[j] = 0;
						}
					}

					// S_i2 quantities, now only record i2's info in S_i2, will update more later
					vector<int> S_i2;
					S_i2.reserve(I[t]);
					S_i2.push_back(i2);
					int N_Si2 = S_i2.size();
					int n0_Si2[J[t]];
					int n1_Si2[J[t]];
					for (int j = 0; j < J[t]; ++j) {
						if (y[t][i2][j] == 0) {
							n0_Si2[j] = 1;
							n1_Si2[j] = 0;
							++n0_S[j];
						}
						else if (y[t][i2][j] == 1) {
							n0_Si2[j] = 0;
							n1_Si2[j] = 1;
							++n1_S[j];
						}
						else {
							n0_Si2[j] = 0;
							n1_Si2[j] = 0;
						}
					}
					/*
					cout << "initial S_i1 = " << endl;
					printVector(S_i1);
					cout << "initial S_i2 = " << endl;
					printVector(S_i2);
					*/

					// permute the indices of the obs in S 
					SampleNoReplace(NS, NS, I_inds_perm[t], I_array[t]);

					// quantities needed for the MH acceptance ratio
					double like_term1 = prodbetabinomG0(y[t][i1], lambda1[t], lambda0[t], J[t]);//log scale
					double like_term2 = prodbetabinomG0(y[t][i2], lambda1[t], lambda0[t], J[t]);//log scale
					double like_term3 = like_term1;//log scale
					like_term3 = like_term3 + prodbetabinomH(y[t][i2], lambda1[t],
						lambda0[t], n1_Si1, n0_Si1, J[t]);//log scale

					  // add obs in S to either S_i1 or S_i2 in order of permutation
					  // with probs given in Dahl eq 12
					int Si2_indic[NS]; // 0 if obs in Si1, 1 if obs in Si2
					double trans_prob = 0.0; // transition probability //log scale  std::log(1.0)
					for (int iiprime = 0; iiprime < NS; ++iiprime) {
						int iprime = Svec[I_inds_perm[t][iiprime]];
						//iprime is one individual in S that are not i1 or i2

						like_term3 = like_term3 + prodbetabinomH(y[t][iprime], lambda1[t], lambda0[t], n1_S, n0_S, J[t]); //log scale

						// calculate prob of being put in cluster i1
						double term1 = std::log(N_Si1) + prodbetabinomH(y[t][iprime],
							lambda1[t], lambda0[t],
							n1_Si1, n0_Si1, J[t]); //log scale
					  //Rcpp::Rcout<<"term1 is " <<term1<<'\n';
						double term2 = std::log(N_Si2) + prodbetabinomH(y[t][iprime],
							lambda1[t], lambda0[t],
							n1_Si2, n0_Si2, J[t]); //log scale
					   //Rcpp::Rcout<<"term2 is " <<term2<<'\n';
						double term_max = std::max(term1, term2);
						term1 = term1 - term_max;
						term2 = term2 - term_max;
						double denom = std::exp(term1) + std::exp(term2);//back to decimal scale to compute the weight sum
						const double prob_i1 = std::exp(term1) / denom;
						//Rcpp::Rcout<<"prob_i1 is " <<prob_i1<<'\n';

						if (R::runif(0, 1) < prob_i1) { // put obs iprime in S_i1
							S_i1.push_back(iprime);
							++N_Si1;
							Si2_indic[I_inds_perm[t][iiprime]] = 0;
							trans_prob = std::log(prob_i1) + trans_prob; //log scale
							like_term1 = like_term1 + prodbetabinomH(y[t][iprime], lambda1[t],
								lambda0[t], n1_Si1, n0_Si1, J[t]); //log scale
							for (int j = 0; j < J[t]; ++j) {
								if (y[t][iprime][j] == 0) {
									++n0_Si1[j];
									++n0_S[j];
								}
								else if (y[t][iprime][j] == 1) {
									++n1_Si1[j];
									++n1_S[j];
								}
							}
						}
						else { // put obs iprime in S_i2
							S_i2.push_back(iprime);
							++N_Si2;
							Si2_indic[I_inds_perm[t][iiprime]] = 1;
							trans_prob = std::log(1.0 - prob_i1) + trans_prob;//log scale
							like_term2 = like_term2 + prodbetabinomH(y[t][iprime], lambda1[t],
								lambda0[t], n1_Si2, n0_Si2, J[t]); //log scale
							for (int j = 0; j < J[t]; ++j) {
								if (y[t][iprime][j] == 0) {
									++n0_Si2[j];
									++n0_S[j];
								}
								else if (y[t][iprime][j] == 1) {
									++n1_Si2[j];
									++n1_S[j];
								}
							}
						}

					} // end iiprime loop // finishing temprorarily splitting the set S, then decide to accept or reject
					/*
					cout << "candidate S_i1 = " << endl;
					printVector(S_i1);
					cout << "candidate S_i2 = " << endl;
					printVector(S_i2);
					*/

					// compute MH ratio and accept with this prob 	  
					const double like_ratio = like_term1 + like_term2 - like_term3;//log scale

					/*Rcpp::Rcout << "(log scale)  split like_term1 is " << like_term1 << '\n';
					Rcpp::Rcout << "(log scale)  split like_term2 is " << like_term2 << '\n';
					Rcpp::Rcout << "(log scale)  split like_term3 is " << like_term3 << '\n';*/
					const double prior_ratio = alpha_log[t] + std::lgamma(N_Si1) + std::lgamma(N_Si2) - std::lgamma(N_Si1 + N_Si2);//log scale

					const double trans_ratio = std::log(1.0) - trans_prob;//log scale equation(3.13) in Jain and Neal
					// because 
					// Pr(merge) = 1
					/*Rcpp::Rcout << "(log scale)  split like_ratio is " << like_ratio << '\n';
					Rcpp::Rcout << "(log scale)  split prior_ratio is " << prior_ratio << '\n';
					Rcpp::Rcout << "(log scale)  split trans_ratio is " << trans_ratio << '\n';*/
					const double MH_ratio = std::exp(like_ratio + prior_ratio + trans_ratio);//back to decimal scale
					/*Rcpp::Rcout << "(demical scale) split MH_ratio is " << MH_ratio << '\n';*/
					if (R::runif(0, 1) < MH_ratio) {
						// accept split (obs in Si1 keep same index as before split
						//               obs in Si2 get index at top of castack)
						++merge_split_accepts[t];
						++n_split[t];
						// do bookkeeping if a split is made

						// adjust cr
						cr[t][i2] = castack[t].top();
						for (int iiprime = 0; iiprime < NS; ++iiprime) {
							int iprime = Svec[I_inds_perm[t][iiprime]];
							if (Si2_indic[I_inds_perm[t][iiprime]] == 1) {
								cr[t][iprime] = castack[t].top();
							}
						}
						/*
						cout << "new cr after bookkeeping = " << endl;
						printVector(cr);
						*/

						// adjust n_in_cluster
						n_in_cluster[t][cr[t][i1]] = N_Si1;
						n_in_cluster[t][cr[t][i2]] = N_Si2;

						// adjust n_0 and n_1
						for (int j = 0; j < J[t]; ++j) {
							n_0[t][cr[t][i1]][j] = n0_Si1[j];
							n_1[t][cr[t][i1]][j] = n1_Si1[j];
							n_0[t][cr[t][i2]][j] = n0_Si2[j];
							n_1[t][cr[t][i2]][j] = n1_Si2[j];
						}

						// adjust cru
						cru[t] = cr[t];
						sort(cru[t].begin(), cru[t].end());
						cru[t].erase(unique(cru[t].begin(), cru[t].end()), cru[t].end());

						// adjust castack
						castack[t].pop();

						// adjust K
						++K[t];
						/*
						cout << "n_in_cluster after bookkeeping = " << endl;
						for (int k=0; k<K; ++k){
						  cout << n_in_cluster[cru[k]] << " " ;
						}
						*/

					}// merge step for period t ends here
				}
				else { // if cr[t][i1] != cr[t][i2] attempt a merge

				  // let S_i1 and S_i2 denote the clusters containing 
				  // i1 and i2 respectively 	  
				  // attempt to create a merged component S = S_i1 union S_i2
					vector<int> Svec; // Svec does not contain i1 and i2 at the moment
					Svec.reserve(I[t]);
					int n0_S[J[t]];
					int n1_S[J[t]];
					for (int j = 0; j < J[t]; ++j) {
						n0_S[j] = 0;
						n1_S[j] = 0;
					}
					for (int i = 0; i < I[t]; ++i) {
						if (cr[t][i] == cr[t][i1] || cr[t][i] == cr[t][i2]) {
							if (i != i1 && i != i2) {
								Svec.push_back(i);
							}
						}
					}
					int NS = Svec.size(); // number of elements in S

					// form initial singleton sets S_i1 and S_i2
					// S_i1 quantities
					vector<int> S_i1;
					S_i1.reserve(I[t]);
					S_i1.push_back(i1);
					int N_Si1 = S_i1.size();
					int n0_Si1[J[t]];
					int n1_Si1[J[t]];
					for (int j = 0; j < J[t]; ++j) {
						if (y[t][i1][j] == 0) {
							n0_Si1[j] = 1;
							n1_Si1[j] = 0;
							++n0_S[j];
						}
						else if (y[t][i1][j] == 1) {
							n0_Si1[j] = 0;
							n1_Si1[j] = 1;
							++n1_S[j];
						}
						else {
							n0_Si1[j] = 0;
							n1_Si1[j] = 0;
						}
					}

					// S_i2 quantities
					vector<int> S_i2;
					S_i2.reserve(I[t]);
					S_i2.push_back(i2);
					int N_Si2 = S_i2.size();
					int n0_Si2[J[t]];
					int n1_Si2[J[t]];
					for (int j = 0; j < J[t]; ++j) {
						if (y[t][i2][j] == 0) {
							n0_Si2[j] = 1;
							n1_Si2[j] = 0;
							++n0_S[j];
						}
						else if (y[t][i2][j] == 1) {
							n0_Si2[j] = 0;
							n1_Si2[j] = 1;
							++n1_S[j];
						}
						else {
							n0_Si2[j] = 0;
							n1_Si2[j] = 0;
						}
					}
					/*
					cout << "initial S_i1 = " << endl;
					printVector(S_i1);
					cout << "initial S_i2 = " << endl;
					printVector(S_i2);
					*/

					// permute the indices of the obs in S 
					SampleNoReplace(NS, NS, I_inds_perm[t], I_array[t]);

					// quantities needed for the MH acceptance ratio
					double like_term1 = prodbetabinomG0(y[t][i1], lambda1[t], lambda0[t], J[t]);//log scale
					double like_term2 = prodbetabinomG0(y[t][i2], lambda1[t], lambda0[t], J[t]);//log scale
					double like_term3 = like_term1;//log scale
					like_term3 = like_term1 + prodbetabinomH(y[t][i2], lambda1[t],
						lambda0[t], n1_Si1, n0_Si1, J[t]);//log scale

					// add obs in S to either S_i1 or S_i2 in order of permutation
					// with probs given in Dahl eq 12
					int Si2_indic[NS]; // 0 if obs in Si1 1 if obs in Si2
					double trans_prob = 0.0; // transition probability //log scale log(1.0)=0.0
					for (int iiprime = 0; iiprime < NS; ++iiprime) {
						int iprime = Svec[I_inds_perm[t][iiprime]];


						like_term3 = like_term3 + prodbetabinomH(y[t][iprime], lambda1[t],
							lambda0[t], n1_S, n0_S, J[t]);//log scale

						// calculate prob of being in cluster i1
						double term1 = std::log(N_Si1) +
							prodbetabinomH(y[t][iprime],
								lambda1[t], lambda0[t],
								n1_Si1, n0_Si1, J[t]);
						//Rcpp::Rcout<<"term1 is " <<term1<<'\n';
						double term2 = std::log(N_Si2) +
							prodbetabinomH(y[t][iprime],
								lambda1[t], lambda0[t],
								n1_Si2, n0_Si2, J[t]);
						// Rcpp::Rcout<<"term2 is " <<term2<<'\n';
						double term_max = std::max(term1, term2);
						term1 = term1 - term_max;
						term2 = term2 - term_max;
						double denom = std::exp(term1) + std::exp(term2);
						const double prob_i1 = std::exp(term1) / denom;
						//Rcpp::Rcout<<"prob_i1 is " <<prob_i1<<'\n';

						if (cr[t][iprime] == cr[t][i1]) { // put obs iprime in S_i1 // this is imaginary split, so that I can compute the transition ratio
							S_i1.push_back(iprime);
							++N_Si1;
							Si2_indic[I_inds_perm[t][iiprime]] = 0;
							trans_prob = std::log(prob_i1) + trans_prob;//log scale
							like_term1 = like_term1 + prodbetabinomH(y[t][iprime], lambda1[t],
								lambda0[t], n1_Si1,
								n0_Si1, J[t]);
							for (int j = 0; j < J[t]; ++j) {
								if (y[t][iprime][j] == 0) {
									++n0_Si1[j];
									++n0_S[j];
								}
								else if (y[t][iprime][j] == 1) {
									++n1_Si1[j];
									++n1_S[j];
								}
							}
						}
						else { // put obs iprime in S_i2 // this is imaginary split, so that I can compute the transition ratio
							S_i2.push_back(iprime);
							++N_Si2;
							Si2_indic[I_inds_perm[t][iiprime]] = 1;
							trans_prob = std::log(1.0 - prob_i1) + trans_prob;
							like_term2 = like_term2 + prodbetabinomH(y[t][iprime], lambda1[t],
								lambda0[t], n1_Si2,
								n0_Si2, J[t]);//log scale
							for (int j = 0; j < J[t]; ++j) {
								if (y[t][iprime][j] == 0) {
									++n0_Si2[j];
									++n0_S[j];
								}
								else if (y[t][iprime][j] == 1) {
									++n1_Si2[j];
									++n1_S[j];
								}
							}
						}

					} // end iiprime loop
					/*
					cout << "candidate S_i1 = " << endl;
					printVector(S_i1);
					cout << "candidate S_i2 = " << endl;
					printVector(S_i2);
					*/

					// compute MH ratio and accept with this prob 	  
					const double like_ratio = like_term3 - like_term1 - like_term2;
					/*Rcpp::Rcout << "(log scale) merge like_term3 is " << like_term3 << '\n';
					Rcpp::Rcout << "(log scale) merge like_term1 is " << like_term1 << '\n';
					Rcpp::Rcout << "(log scale) merge like_term2 is " << like_term2 << '\n';*/
					const double prior_ratio = -alpha_log[t] - std::lgamma(N_Si1) - std::lgamma(N_Si2) + std::lgamma(N_Si1 + N_Si2);//log scale




					const double trans_ratio = trans_prob;
					/*Rcpp::Rcout << "(log scale) merge like_ratio is " << like_ratio << '\n';
					Rcpp::Rcout << "(log scale) merge prior_ratio is " << prior_ratio << '\n';
					Rcpp::Rcout << "(log scale) merge trans_ratio is " << trans_ratio << '\n';*/
					const double MH_ratio = std::exp(like_ratio + prior_ratio + trans_ratio);//back to decimal scale
					/*Rcpp::Rcout << "(decimal scale) merge MH_ratio is " << MH_ratio << '\n';*/
					if (R::runif(0, 1) < MH_ratio) {
						// accept merge (all obs get cr[i1])
						++merge_split_accepts[t];
						++n_merge[t];
						// do bookkeeping if a merge is made

						// adjust cr
						int holder = cr[t][i2];
						cr[t][i2] = cr[t][i1];
						for (int iiprime = 0; iiprime < NS; ++iiprime) {
							int iprime = Svec[I_inds_perm[t][iiprime]];
							if (Si2_indic[I_inds_perm[t][iiprime]] == 1) {
								cr[t][iprime] = cr[t][i1];
							}
						}
						/*
						cout << "new cr after bookkeeping = " << endl;
						printVector(cr);
						*/

						// adjust n_in_cluster
						n_in_cluster[t][cr[t][i1]] = N_Si1 + N_Si2;
						n_in_cluster[t][holder] = 0;

						// adjust n_0 and n_1
						for (int j = 0; j < J[t]; ++j) {
							n_0[t][cr[t][i1]][j] = n0_Si1[j] + n0_Si2[j];
							n_1[t][cr[t][i1]][j] = n1_Si1[j] + n1_Si2[j];
						}

						// adjust cru
						cru[t] = cr[t];
						sort(cru[t].begin(), cru[t].end());
						cru[t].erase(unique(cru[t].begin(), cru[t].end()), cru[t].end());

						// adjust castack
						castack[t].push(holder);

						// adjust K
						--K[t];
						/*
						cout << "n_in_cluster after bookkeeping = " << endl;
						for (int k=0; k<K; ++k){
						  cout << n_in_cluster[cru[k]] << " " ;
						}
						*/

					}

				}
			}
			
		} //split-merge step for all periods ends here
		else {
			for (int t = 0; t < T; ++t) {
				SampleNoReplace(I[t], I[t], I_inds_perm[t], I_array[t]);
				for (int ii = 0; ii < I[t]; ++ii) {
					const int i = I_inds_perm[t][ii];

					double numer[K[t] + 1];   // numerator for cluster weights, log scale
					double denom = 0.0; // denominator for cluster weights, log scare
					double numer_max;    // rescaling factor
					double clusterweights[K[t] + 1]; // prob of being in each cluster
					// last element is prob of new cluster

					//I use Log-sum-exp Trick.
					//I always use the new cluster's weight on log scale as the correction term
					//This way the new cluster weight will always be exp(0)=1. Thus, the denominator will not be zero.
					numer_max = alpha_log[t] + prodbetabinomG0(y[t][i], lambda1[t], lambda0[t], J[t]); // current maximum 
					numer[K[t]] = 1.0;//exp(numer_max-numer_max)=exp(0)=1

					for (int kk = 0; kk < K[t]; ++kk) {
						int k = cru[t][kk];
						if (k != cr[t][i]) {
							numer[kk] = n_in_cluster[t][k] *
								std::exp(prodbetabinomH(y[t][i], lambda1[t], lambda0[t],
									n_1[t][k], n_0[t][k], J[t]) - numer_max); //log scale
						}
						else if (n_in_cluster[t][k] > 1) {
							numer[kk] = (n_in_cluster[t][k] - 1.0)*
								std::exp(prodbetabinomH(y[t][i], lambda1[t], lambda0[t],
									n_1[t][k], n_0[t][k], J[t]) - numer_max); //log scale
						}
						else {
							numer[kk] = 0.0;
						}
						// calculate the current maximum numerator
					   // numer_max = hpmax(numer[kk], numer_max);
					}// end k loop

					// rescale numer and form denom
					for (int kk = 0; kk < (K[t] + 1); ++kk) {
						//numer[kk] = hpdiv(numer[kk], numer_max);
						denom += numer[kk];
					}
					// put numer and denom together to form weights
					for (int kk = 0; kk < (K[t] + 1); ++kk) {
						clusterweights[kk] = numer[kk] / denom;
						//Rcpp::Rcout << "peirod " << t << " cluster " << kk << " prob is " << clusterweights[kk];
					}

					// sample clusterindexp[i]
					// cluster_i is in {0, 1, ..., K}
					// this needs to be mapped back onto the appropriate row number
					
					const int cluster_i = ProbSamp(Iseq_array[t], clusterweights,
						K[t] + 1);
					//Rcpp::Rcout << "period " << t << " unit " << i << " 'new cluster is " << cluster_i << endl;
					// cr_i is an element of {cru, top of available castack} 
					int cr_i;
					if (cluster_i < K[t]) { // cluster_i is an existing cluster
						cr_i = cru[t][cluster_i];
					}
					else {
						cr_i = castack[t].top();
					}

					// do the necessary bookkeeping
					//
					// cases: 
					//        1) obs i stays in same cluster 
					//           (cr_i == cr[i])
					//
					//        2) obs i moves to other existing cluster and obs i's 
					//           old cluster remains non-empty
					//           (cr_i != cr[i] &
					//            cr_i != castack.top() & 
					//            n_in_cluster[cr[i]] > 1)
					//
					//        3) obs i moves to other existing cluster and obs i's 
					//           old cluster becomes empty
					//           (cr_i != cr[i] &
					//            cr_i != castack.top() & 
					//            n_in_cluster[cr[i]] == 1)
					//
					//        4) obs i forms a new cluster and obs i's 
					//           old cluster remains non-empty
					//           (cr_i == castack.top() & 
					//            n_in_cluster[cr[i]] > 1)
					//
					//        5) obs i forms a new cluster and obs i's 
					//           old cluster becomes empty
					//           (cr_i == castack.top() & 
					//            n_in_cluster[cr[i]] == 1)

					// case 2
					if (cr_i != cr[t][i] && cr_i != castack[t].top() &&
						n_in_cluster[t][cr[t][i]] > 1) {

						// adjust n_0 and n_1
						for (int j = 0; j < J[t]; ++j) {
							if (y[t][i][j] == 0) {
								--n_0[t][cr[t][i]][j];
								++n_0[t][cr_i][j];
							}
							if (y[t][i][j] == 1) {
								--n_1[t][cr[t][i]][j];
								++n_1[t][cr_i][j];
							}
						}

						// adjust n_in_cluster
						--n_in_cluster[t][cr[t][i]];
						++n_in_cluster[t][cr_i];

						// adjust cr
						//const int holder = cr[t][i]; don't need this holder, because the castask is unchanged in this case
						cr[t][i] = cr_i;

						// adjust cru
						cru[t] = cr[t];
						sort(cru[t].begin(), cru[t].end());
						cru[t].erase(unique(cru[t].begin(), cru[t].end()), cru[t].end());

						// adjust castack (not necessary in this case castack unchanged)

						// adjust K (not necessary in this case-- K remains the same)

					} // end case 2 if
					// case 3
					if (cr_i != cr[t][i] && cr_i != castack[t].top() &&
						n_in_cluster[t][cr[t][i]] == 1) {

						// adjust n_0 and n_1
						for (int j = 0; j < J[t]; ++j) {
							if (y[t][i][j] == 0) {
								--n_0[t][cr[t][i]][j];
								++n_0[t][cr_i][j];
							}
							if (y[t][i][j] == 1) {
								--n_1[t][cr[t][i]][j];
								++n_1[t][cr_i][j];
							}
						}

						// adjust n_in_cluster
						--n_in_cluster[t][cr[t][i]];
						++n_in_cluster[t][cr_i];

						// adjust cr
						const int holder = cr[t][i];
						cr[t][i] = cr_i;

						// adjust cru
						cru[t] = cr[t];
						sort(cru[t].begin(), cru[t].end());
						cru[t].erase(unique(cru[t].begin(), cru[t].end()), cru[t].end());

						// adjust castack
						castack[t].push(holder);

						// adjust K
						--K[t];

					} // end case 3 if
					// case 4 
					else if (cr_i == castack[t].top() && n_in_cluster[t][cr[t][i]] > 1) {

						// adjust n_0 and n_1
						for (int j = 0; j < J[t]; ++j) {
							if (y[t][i][j] == 0) {
								--n_0[t][cr[t][i]][j];
								++n_0[t][cr_i][j];
							}
							if (y[t][i][j] == 1) {
								--n_1[t][cr[t][i]][j];
								++n_1[t][cr_i][j];
							}
						}

						// adjust n_in_cluster
						--n_in_cluster[t][cr[t][i]];
						++n_in_cluster[t][cr_i];

						// adjust cr
						//const int holder = cr[t][i]; no need for holder in this case either
						cr[t][i] = cr_i;

						// adjust cru
						cru[t] = cr[t];
						sort(cru[t].begin(), cru[t].end());
						cru[t].erase(unique(cru[t].begin(), cru[t].end()), cru[t].end());

						// adjust castack
						castack[t].pop();

						// adjust K
						++K[t];

					}// end case 4 if
					// case 5 
					else if (cr_i == castack[t].top() && n_in_cluster[t][cr[t][i]] == 1) {

						// adjust n_0 and n_1
						for (int j = 0; j < J[t]; ++j) {
							if (y[t][i][j] == 0) {
								--n_0[t][cr[t][i]][j];
								++n_0[t][cr_i][j];
							}
							if (y[t][i][j] == 1) {
								--n_1[t][cr[t][i]][j];
								++n_1[t][cr_i][j];
							}
						}

						// adjust n_in_cluster
						--n_in_cluster[t][cr[t][i]];
						++n_in_cluster[t][cr_i];

						// adjust cr
						const int holder = cr[t][i];
						cr[t][i] = cr_i;

						// adjust cru
						cru[t] = cr[t];
						sort(cru[t].begin(), cru[t].end());
						cru[t].erase(unique(cru[t].begin(), cru[t].end()), cru[t].end());

						// adjust castack
						castack[t].pop();
						castack[t].push(holder);

						// adjust K (not necessary in this case-- K remains the same)

					}// end case 5 if

				} // end i loop, to draw updated clusters for every individual
				//Rcpp::Rcout << "period " << t << " has number of clusters " << K[t] << std::endl;
			}//period t loop ends
			

		} // a new draw news here, including drawing new cluster or split-merge



	  // if iter is not evenly divisible by *sams_iter use Gibbs
	  // permute indices of individuals for better sampling behavior



	  

		// sample alpha: alpha has a log normal distribution. log(alpha) has a normal distribution
		// we use a MH random walk algorithm to sample alpha_t 
		//Rcpp::rgamma is parameterized by shape and scale, so I use 1/rate=scale
		
		for (int t = 0; t < T; ++t) {
			double step = R::runif(-tune[t], tune[t]);
			double new_alpha = alpha[t] + step;
			//new_alpha has to be positive
			while(new_alpha<=0){
				step = R::runif(-tune[t], tune[t]);
				new_alpha = alpha[t] + step;
			}
			// The likelihood for alpha_t \propto \alpha_t^{K[t]}\frac{\Gamma(alpha_t)}{\Gamma(alpha_t)+I[t]}
			// The prior for alpha_t \propto \frac{1}{\alpha_t} \exp{\frac{(\ln alpha_t - gamma_t)}{2 V}}, which is a log-normal distribution.
			// double alpha_logLik(const double& alpha, const int& n, const int& k){
						// return k*std::log(alpha)+std::lgamma(alpha)-std::lgamma(alpha+n);
					// }
			
			// The Rcpp sugar functions are meant for vector type arguments like Rcpp::NumericVector. 
			// For scalar arguments you can use the functions in the R namespace:
			double ratio = alpha_logLik(new_alpha, I[t], K[t])-alpha_logLik(alpha[t], I[t], K[t])
			+R::dlnorm(new_alpha, gamma[t],std::sqrt(V),true)
			-R::dlnorm(alpha[t], gamma[t],std::sqrt(V),true);
			//Rcpp::Rcout << "The MH ratio is " << std::exp(ratio) << std::endl;
			if (R::runif(0, 1) < std::exp(ratio)) {
				alpha[t] = new_alpha;
				alpha_log[t] = std::log(alpha[t]);
				MH_accept[t]+=1;
			}
		}
		
		
		
		//use DLM to smooth all gamma[t]s
		//alpha_log[t] \sim N(\gamma[t], V)
		//\gamma[t] \sim N(\gamma[t-1], W)
        // known parameter values
		
		//## run the Kalman filter forward through time
		a[0] = m0;
		R[0] = C0 + W;
		f[0] = a[0];
		Q[0] = R[0] + V;
		e[0] = alpha_log[0] - f[0];
		A[0] = R[0] / Q[0];
		m[0] = a[0] + A[0] * e[0];
		C[0] = R[0] - A[0] * Q[0] * A[0];
		for (int i = 1; i < T; ++i) { // C++index starts with 0
			a[i] = m[i - 1];
			R[i] = C[i - 1] + W;
			f[i] = a[i];
			Q[i] = R[i] + V;
			e[i] = alpha_log[i] - f[i];
			A[i] = R[i] / Q[i];
			m[i] = a[i] + A[i] * e[i];
			C[i] = R[i] - A[i] * Q[i] * A[i];
		}

		// sample the last gamma, gamma[T-1]
		gamma[T - 1] = R::rnorm(m[T - 1], std::sqrt(C[T - 1])); //C++ index starts with 0
		
	    //sample gamma backwards through time
		double B = 0;
		double h = 0;
		double H = 0;
		for (int i = T - 2; i >= 0; --i) { //C++ index starts with 0
			B = C[i] / R[i + 1];
			h = m[i] + B * (gamma[i + 1] - a[i + 1]);
			H = C[i] - B * R[i + 1] * B;
			gamma[i] = R::rnorm(h, std::sqrt(H));
		}

       //sample new V (sigma2)
		double residual[T];
		for (int t = 0; t < T; ++t) {
			residual[t] = alpha_log[t] - gamma[t];
		}
		double sse = 0;
		for (int t = 0; t < T; ++t) {
			sse += std::pow(residual[t], 2.0);
		}
		
	   /* sigma2 < -rinvgamma(1, (alpha + T) / 2, (beta + sse) / 2), paramterized by shape and rate
	    V < -sigma2*/
		//V prior inverse gamma(r0/2,s0/2)
		//V posterior is inverse gamma((r0 + T) / 2,  (s0 + sse)/2),(s0 + sse)/2 is rate
		// which is equivalent to 1 / Rcpp::rgamma(1, (r0 + T) / 2, 2 / (s0 + sse))[0]; I checked this in R
		V = 1 / R::rgamma((r0 + T) / 2, 2 / (s0 + sse));
		//}
		
		if(fix_W != 1){
			
			// sample W			
			
			double step_W = R::runif(-tau_W, tau_W);
			double new_W =  W + step_W;
			//new_W has to be positive
			while(new_W<=0){
				step_W = R::runif(-tau_W, tau_W);
				new_W =  W + step_W;
			}
			// The likelihood for W \propto \sum_{t=1}^{t=T} \phi(gamma_{t}-gamma{t-1}; 0, W), where \phi(., 0, W) is the PDF of a normal distribution with mean 0 and variance W. 
			// W \sim Gamma(shape = r1, rate=s1), so s1^r1 is on the numerator.
			//Rcpp::rgamma is parameterized by shape and scale, so I use 1/rate=scale
			
			double W_logLik_old = 0;
			double W_logLik_new = 0;
			for(int t = 1; t < T; ++t){
				
				W_logLik_old += R::dnorm(gamma[t]-gamma[t-1],0,std::sqrt(W),true);
				W_logLik_new += R::dnorm(gamma[t]-gamma[t-1],0,std::sqrt(new_W),true);
			}
			
			// The Rcpp sugar functions are meant for vector type arguments like Rcpp::NumericVector. 
			// For scalar arguments you can use the functions in the R namespace:
			double W_ratio = W_logLik_new-W_logLik_old+R::dgamma(new_W, r1,1/s1,true)-R::dgamma(W, r1,1/s1,true);
			//Rcpp::Rcout << "The MH ratio is " << std::exp(ratio) << std::endl;
			if (R::runif(0, 1) < std::exp(W_ratio)) {
				W = new_W;
				MH_accept_W+=1;
			}
			
		}
		

		if (iter >= burnin && ((iter % thin) == 0)) {
			for (int t = 0; t < T; ++t) {
				nclust_store_data(count, t) = K[t];
				alpha_store_data(count, t) = alpha[t];
				gamma_store_data(count, t) = gamma[t];
				// b0_store(count,t) = b0[t];
				
				for (int i = 0; i < I[t]; ++i) {
					c_store(count, i, t) = cr[t][i];
				}
			}
			V_store_data(count)=V;
			W_store_data(count)=W;
	

			++count;
		}

		if ((iter+1) % verbose == 0) {
			Rcpp::Rcout << "iteration = "<<iter+1 << endl;
			Rcpp::checkUserInterrupt();
		}

	} // end an entire iteration


    //report the HM acceptance rate for each period
	//For every iteration, not matter with split-merge or a normal iteration, we still do 
	arma::Col<double> accept_rate(T);
	for (int t = 0; t < T; ++t) {
		accept_rate(t) = (double)MH_accept[t] / tot_iter;
		
	}
	double accept_rate_W = (double)MH_accept_W/tot_iter;

	///////Passing output to R
	return Rcpp::List::create(Rcpp::Named("cluster_label_store") = c_store,
		Rcpp::Named("n_cluster_store_data") = nclust_store_data,
		Rcpp::Named("merge_split_accept") = merge_split_accepts,
		Rcpp::Named("merge_split_count") = merge_split_count,
		Rcpp::Named("n_split") = n_split,
		Rcpp::Named("n_merge") = n_merge,
		Rcpp::Named("alpha_store_data") = alpha_store_data,
		Rcpp::Named("gamma_store_data") = gamma_store_data,
		Rcpp::Named("V_store_data") = V_store_data,
		Rcpp::Named("W_store_data") = W_store_data,
		Rcpp::Named("alpha_accept_rate")=accept_rate,
		Rcpp::Named("W_accept_rate")=accept_rate_W);
		
}



/*** R
print("C++ file compiled")
*/
