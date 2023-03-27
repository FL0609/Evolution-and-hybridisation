#include "fisher.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <ctime>
#include <cmath>
#include <vector>
#include <algorithm>
#include <chrono>
#include <random>
#include <array>
using namespace std;

// Constants:
#define PIE 3.14159265359
#define EXITFLAG -1

// Hardcoded variables (need to be changed by hand and code recompiled):
#define ALPHA 0.5 // overall strength of selection
#define BN_LENGTH 25 // Expected number of neutral fixations before population size increase in the bottleneck model (Ev=-1)
#define BN_DEPTH 10 // Population declines by this factor with the bottleneck model (Ev=-1)
#define CYCLELENGTH 100000 // how many generations does it take for the moving optimum to complete a cycle (Ev=1)
#define UPDATELENGTH 10000 // Report to the screen after this many generations

extern MTRand rnd;

/*****************************************************************************************************************************************
 *****************************************************************************************************************************************
 *****************************************************************************************************************************************
 ******** FUNCTION DEFINITIONS
 *****************************************************************************************************************************************
 *****************************************************************************************************************************************
 *****************************************************************************************************************************************/

/*************************************************************************
 ******** 	Return the position of the optimum at generation "gen"
 ******** 	Ev=-1: always fixed at 0 (the ancestral trait value)
 ******** 	Ev=0: always fixed at 0 (the ancestral trait value)
 ******** 	Ev=2: always fixed at 1
 ******** 	Ev=1: oscillating pattern over "CYCLELENGTH" generations
 *************************************************************************/
double get_opt(int gen, int Ev, int cyclength, double kv)
{
	double tmp;

	if(Ev==0 || Ev== -1) // Fixed optimum at the origin
	{
		return(0.0);
	}
	if(Ev==2) // Optimum fixed at a unit length
	{
		return(1.0);
	}
	if(Ev==1)
	{
		tmp = 0.5 + sin(2.0*PIE*gen/cyclength-PIE/2)*0.5;
		if(tmp>0.0)
			return(pow(tmp,2.0/kv));
		if(tmp<0.0)
			return(-pow(-tmp,2.0/kv));
		if(tmp==0.0)
			return(0.0);
	}
	exit(1);
}

/*****************************************************************************************************************************************
 *****************************************************************************************************************************************
 *****************************************************************************************************************************************
 ******** The main recursion
 *****************************************************************************************************************************************
 *****************************************************************************************************************************************
 *****************************************************************************************************************************************

Individual-based simulation of evolutionary divergence for a single population under Fisher's geometrical model,
in a range of population genetic regimes

Input parameters:

** Population genetics parameters
Nv: size of the population (post-bottleneck if implemented); Number of diploid individuals
Lv: genome map length in Morgans (mean nb of cross-overs per meiosis); set to -1 for free recombination
Uv: mutation rate per diploid genome

** Fitness landscape
kv: curvature of fitness function
nv: number of phenotypic traits

** Distribution of mutation effects
Mv: Flag for mutational effect model: 1 "bottom up"; 0 "top down"
sv: The mean selection coefficient for random mutations in an optimal genotype
qv: mean phenotypic dominance effect (set to 1/2 in published work)
Fv: inflation of variance of distribution of phenotypic dominance (set to -1 for a uniform distribution)

** Demographic/environmental change
Ev: Flag for environmental change (i.e. position of the optimum):
	-1 (optimum is fixed to match ancestral state + initial bottleneck)
	0  (optimum is fixed to match ancestral state)
	1  (oscillating over time)
	2  (fixed at 1, away from ancestral state)

** Simulation parameters
Dv: Stop simulation after Dv fixations have occurred.

 *****************************************************************************************************************************************
 *****************************************************************************************************************************************
 *****************************************************************************************************************************************/
void recursion( int Nv, double Lv, double Uv,
                double kv, double sv, int nv,
				double Fv, double qv, int Mv,
				int Ev, int Dv, int repv)
{
	// declare variables:
	int par1, par2, nb, mut, indmut, index, N_1, fourN;
	double w, d, sz2, mutsize, pos, Wmax, mydom;
	bool withrec;
	long long int maxgen;

   	long long int gen = 0;
	int nb_seg = 0; //number of segregating mutations
	int nb_fix = 0;	//number of fixed variants

	double kd2 = kv / 2.0; // half of the fitness landscape curvature
	withrec = Lv == -1 ? false : true; // Flag for free recombination
	bool mono = true; //is the population monomorphic? skip to mutation step

	// compute various quantities (applying mainly to the bottleneck model: Ev=-1)
	int genNup = (Ev==-1) ? BN_LENGTH / Uv : 0;
	int Nsmall = Nv / BN_DEPTH; //population size before increase
	int N_1small = Nsmall - 1;
	int fourNsmall = 4*Nsmall;
	int Nlarge = Nv; //population size after increase
	int N_1large = Nlarge - 1;
	int fourNlarge = 4*Nlarge;
	double optval = get_opt(((gen) % CYCLELENGTH),Ev,CYCLELENGTH,kv);

	// set up population of diploid invididuals
    //chr * xxx = new chr [fourNlarge]; //store all genotypes in the population
	chr * pop = new chr [fourNlarge]; //store all genotypes in the population
	chr * temp = new chr [fourNlarge]; // used to generate next generation
	chr * cp;

	// "Wtot" will hold the fitness of each individual:
	double * Wtot = new double [Nlarge];

	//set up lists to store information on mutations
	vector<double> mutline((2 + nv), 0.0); //store mutational effects of variants for each trait. elements are: position, mutational effects, allele count
	vector<vector <double> > mutations(0, mutline);
	vector<double> domline(nv, 0.0);
	vector<vector <double> > delta(0, domline); //store dominance coefficients of variants for each trait
	vector<int> muttimes; //store generation at which variant arises
	vector<int> fixtimes; //store generation at which variant fixes
	vector<int> fixed; //store indices of fixed variants
	vector<int> lost; //store indices of lost variants
	vector<int> homo; //store indices of homozygous variants
	vector<int> hetero; // store indices of heterozygous variants
	vector<double>currentmutation(nv,0.0);
	vector<double>wildtype(nv,0.0);

    vector<int> c = {0, 1, 2, 3};
	//reserve capacity for vectors that will grow
	//delta.reserve(100 * nv);
	//mutations.reserve(100 * nv);

	//record duration of simulation
	time_t start, finish;
	//struct tm *ptr;
	start = time(0);


	//set up output file
	string fileName;
	stringstream nameF;
	nameF << "res_N" << Nv <<
					"_L" << Lv << "_U" << Uv <<
					"_k" << kv << "_s" << sv << "_n" << nv << "_F" << Fv << "_q" << qv << "_M" << Mv <<
					"_E" << Ev << "_rep" << repv << ".txt";
	nameF >> fileName;

	//make header output file
	ofstream fout;
	fout.open(fileName);
        fout << "position\tmut_time\tfix_time\topt_dist\t";
        for (int i = 0; i < nv; i++)
        {
                fout << "mut_effect" << i << "\t" << "dom" << i << "\t";
        }
	fout << endl;

	// set parameter for the distribution of mutation effect sizes
	if(Mv) //Multivariate normal
		mutsize = get_mutsize_MVN(nv, kv, sv, ALPHA);
	else //Exponential
		mutsize = get_mutsize_EXP(2.0, sv, ALPHA); // r_bar

	//set population size before increase
	Nv = Nsmall;
    N_1 = N_1small;
    fourN = fourNsmall;

	//   ---------------------------   START RECURSION    ---------------------------
	// loop through generations:
	while(gen != EXITFLAG)
	{

		if(gen == genNup) //in the generation of population size increase, recompute fitnesses
			mono = 0;


		//   (1) ---------------------------   COMPUTE FITNESS    ---------------------------
		if (mono == 0) //if the population is monomorphic, skip to mutation step. if not: compute fitnesses
		{
			Wmax = 0; //maximum fitness in the population
			optval = get_opt(((gen) % CYCLELENGTH),Ev,CYCLELENGTH,kv); //distance to optimum

			for (int j = 0; j < Nv; j++)  // for each individual
			{
				nb = 4 * j; //(nb = index of 1st chromosome individual j, nb+1 = index of 2nd chromosome individual j)
				sz2 = 0; //summed squared distance to optimum
				d = 0; //distance to optimum

				sort(pop[nb].sel.begin(), pop[nb].sel.end());
				sort(pop[nb + 1].sel.begin(), pop[nb + 1].sel.end());
                sort(pop[nb + 2].sel.begin(), pop[nb + 2].sel.end());
                sort(pop[nb + 3].sel.begin(), pop[nb + 3].sel.end());
                
                /* Print the mutations of each chromosome at each generaiton:
                
                cout << "For chromosome " << nb << ", at generation " << gen << ", the mutations are" << endl;
                for (int i = 0; i < pop[nb].sel.size(); i++)
               {
                    cout <<  pop[nb].sel[i] << endl;
               }
               cout << "For chromosome " << nb+1 << ", at generation " << gen << ", the mutations are" << endl;
                for (int i = 0; i < pop[nb+1].sel.size(); i++)
               {
                    cout <<  pop[nb+1].sel[i] << endl;
               } */
                   
 
				//remove fixed mutations from genotype
				for (unsigned int i = 0; i < fixed.size(); i++)
				{
					if (find(pop[nb].sel.begin(), pop[nb].sel.end(), fixed[i]) == pop[nb].sel.end())
						cout << "'fixed' mut " << fixed[i] << " absent in chrom 1 ind " << j << endl;
					if (find(pop[nb + 1].sel.begin(), pop[nb + 1].sel.end(), fixed[i]) == pop[nb + 1].sel.end())
						cout << "'fixed' mut " << fixed[i] << " absent in chrom 2 ind " << j << endl;
                    if (find(pop[nb + 2].sel.begin(), pop[nb + 2].sel.end(), fixed[i]) == pop[nb + 1].sel.end())
						cout << "'fixed' mut " << fixed[i] << " absent in chrom 3 ind " << j << endl;
                    if (find(pop[nb + 3].sel.begin(), pop[nb + 3].sel.end(), fixed[i]) == pop[nb + 1].sel.end())
						cout << "'fixed' mut " << fixed[i] << " absent in chrom 4 ind " << j << endl;
					pop[nb].sel.erase(remove(pop[nb].sel.begin(), pop[nb].sel.end(), fixed[i]), pop[nb].sel.end());
					pop[nb + 1].sel.erase(remove(pop[nb + 1].sel.begin(), pop[nb + 1].sel.end(), fixed[i]), pop[nb + 1].sel.end());	
                   pop[nb + 2].sel.erase(remove(pop[nb + 2].sel.begin(), pop[nb + 2].sel.end(), fixed[i]), pop[nb + 2].sel.end());
                   pop[nb + 3].sel.erase(remove(pop[nb + 3].sel.begin(), pop[nb + 3].sel.end(), fixed[i]), pop[nb + 3].sel.end());					
                }

				//determine homozygous and heterozygous sites
				//set_intersection(pop[nb].sel.begin(), pop[nb].sel.end(),pop[nb + 1].sel.begin(), pop[nb + 1].sel.end(), back_inserter(homo));
				//set_symmetric_difference(pop[nb].sel.begin(), pop[nb].sel.end(),pop[nb + 1].sel.begin(), pop[nb + 1].sel.end(), back_inserter(hetero));
				
                
                
				//compute distance to optimum for each trait
				for (int k = 0; k < nv; k++) //loop through n traits
				{
					d = wildtype[k];
					if (k == 0) //for 1st trait, consider current position of optimum
						d -= optval;

					for (unsigned int i = 0; i < pop[nb].sel.size(); i++) //loop through 1st chromosome
					{
						index = pop[nb].sel[i];
						d += mutations[index][k + 1]*0.25;
					}

					for (unsigned int i = 0; i < pop[nb+1].sel.size(); i++) //loop through 2nd chromosome
					{
						index = pop[nb+1].sel[i] ;
						d += mutations[index][k + 1]*0.25;
					}
                    for (unsigned int i = 0; i < pop[nb+2].sel.size(); i++) //loop through 3rd chromosome
					{
						index = pop[nb+2].sel[i] ;
						d += mutations[index][k + 1]*0.25;
					}
                    for (unsigned int i = 0; i < pop[nb+3].sel.size(); i++) //loop through 4th chromosome
					{
						index = pop[nb+3].sel[i] ;
						d += mutations[index][k + 1]*0.25;
					}

					sz2 += d * d;
				}
             //				homo.clear();
             //				hetero.clear();


				// compute fitness:
				w = getfitness(sz2, kd2, ALPHA);
				Wtot[j] = w;
				if (Wmax < w)
					Wmax = w;
			} //end for loop going through indvidiuals
			// **** bookkeeping ****
			//clear indices of fixed variants
			if (!(fixed.empty()))
			{
				for (unsigned int i = 0; i < fixed.size(); i++)
					lost.push_back(fixed[i]);
				fixed.clear();
			}

			//increase population size
			if (gen == genNup)
			{
				Nv = Nlarge; // in the next for loop we will sample Nv individuals from parents 0 to N_1.
				//By setting Nv to large before the for loop but not N_1, we make sure the increased population of Nlarge offspring are sampled from the current population's Nsmall parents
				fourN = fourNlarge;
			}
			// **** end bookkeeping ****




			//   (2) ---------------------------   SAMPLE NEXT GENERATION    ---------------------------
			for (int ind = 0; ind < Nv; ind++) //loop through offspring
			{
				// sampling first parent:
				do
				{
						par1 = int(rnd.randInt(N_1)); //pick random candidate as parent 1
				} while (rnd.rand()> (Wtot[par1] / Wmax)); //accept candidate with probabiliy proportional to their fitness

				// recombination
				//if (withrec)
				//				rec(temp[2*ind], pop[2*par1], pop[2*par1+1], Lv, mutations, nv);
				//else
                //    {
                        std::random_device rd;
                        std::mt19937 g(rd());
                        std::shuffle(c.begin(),c.end(),g);
//                        cout << "chrom pairing parent " << par1 << " ";
//                        for(unsigned int i=0; i < c.size(); i++)
//                        {
//                            cout << c[i] << ",";
//                        }
//                        cout << endl;
                        
								freerec(temp[4*ind], pop[4*par1+c[0]], pop[4*par1+c[1]]);
                               freerec(temp[4*ind+1], pop[4*par1+c[2]], pop[4*par1+c[3]]);
                             
                  //  }
                //}

				//sample second parent
				do
				{
					par2 = int(rnd.randInt(N_1)); //pick random candidate as parent 1
				} while (rnd.rand()> (Wtot[par2] / Wmax)); //accept candidate with probabiliy proportional to their fitness

				// recombination
				//if (withrec)
				//				rec(temp[2*ind+1], pop[2*par2], pop[2*par2+1], Lv, mutations, nv);
				//else
                //    {
                        random_shuffle(c.begin(), c.end());
//                        cout << "chrom pairing parent " << par2 << " ";
//                        for(int i=0;i<4;i++)
//                        {
//                            cout << c[i] << ",";
//                        }
//                        cout << endl;
                        
//								freerec(temp[4*ind+2], pop[4*par2+0], pop[4*par2+1]);
//                               freerec(temp[4*ind+3], pop[4*par2+2], pop[4*par2+3]);
                               freerec(temp[4*ind+2], pop[4*par2+c[0]], pop[4*par2+c[1]]);
                               freerec(temp[4*ind+3], pop[4*par2+c[2]], pop[4*par2+c[3]]);
                  //  }
                
			}

			// **** bookkeeping ****
			if(gen == genNup) // Now we set N_1 to large as well to make sure all individuals can reproduce from here onwards
				N_1 = N_1large;
			// **** end bookkeeping ****



		} //skip straight to mutation in case of monomorphic population (no need to recompute fitnesses)



		//   (3) ---------------------------   MUTATION    ---------------------------
		mut = int(poisdev(4 * Uv * Nv)); // number of new mutations in entire population

		for (int j = 0; j < mut; j++) //loop through new mutations to distribute over individuals
		{
			indmut = int(rnd.randInt(fourN - 1)); //pick random chromosome to receive mutation
			pos = rnd.rand(); //pick random position in the genome [0,1] ie infinite site model

			// info new mutation
			mutline[0] = pos;
			generatemutation(currentmutation,nv,mutsize,Mv,kv);


			for (int k = 0; k < nv; k++) // generate mutational effect and dominance of new mutation for each trait
			{
				mutline[k + 1] = currentmutation[k];
			}
			mutline[1 + nv] = 0; //this is where the mutation's frequency will go.

			//Generate dominance coefficient(s) of new mutation
			if (Fv == -999) //type of dominance: per mutation
			{
				mydom = getdominance(qv, -1);
				for (int k = 0; k < nv; k++)
				{
					domline[k] = mydom;
				}
			}
			else if (Fv == 123) //dominance coefficients that make large muts recessive and small muts (semi)dominant
		  {
				for (int k = 0; k < nv; k++)
				{
					domline[k] = 2.0 - ( 2.0 / (1.0 + exp(-abs(mutline[k + 1])/mutsize )) );
				}
			}
			else //type of dominance: per trait
			{
				for (int k = 0; k < nv; k++)
				{
					domline[k] = getdominance(qv, Fv);
				}
			}


			if (lost.empty()) //generate new mutation index if none available to recycle
			{
				index = mutations.size();
				delta.push_back(domline);
				mutations.push_back(mutline); //add mutation to list of segregating mutations
				muttimes.push_back(gen); //store generation in which mutation arose
				fixtimes.push_back(0);
			}

			//to keep the keep the size of the mutations list managable, we re-use the position in the mutation list once a mutation is lost
			if (!(lost.empty())) //retrieve unused index if available
			{
				if (lost.size() > 1)
					sort(lost.begin(), lost.end());
				index = lost[0];
				lost.erase(remove(lost.begin(), lost.end(), index), lost.end()); //this index is now no longer available (until mutation is lost)
				mutations[index] = mutline; //add mutation to list of segregating mutations
				delta[index] = domline;
				muttimes[index] = gen; //store generation in which mutation arose
				fixtimes[index] = 0;
			}
			temp[indmut].sel.push_back(index); // add mutation to individual genotype
            
            
               
		}

 
		// **** bookkeeping ****

		// update population:
		cp = pop;
		pop = temp;
		temp = cp;
		// **** end bookkeeping ****





		//   (4) ---------------------------   COMPUTE ALLELE FREQUENCIES    ---------------------------

		if ((mut > 0) || (mono==false)) //if no mutations are segregating (monomorphic population), no need to compute allele frequencies
		{
			mono = true;
			nb_seg = 0;

			//first set all allele counts to 0
			for (unsigned int loc = 0; loc < mutations.size(); loc++)
			{
				mutations[loc][1 + nv] = 0;
			}

			for (int ind = 0; ind < fourN; ind++) //count up alleles of all individuals
			{
				temp[ind].sel.clear();
				if (!(pop[ind].sel.empty()))
				{
					mono = false;
					for (unsigned int loc = 0; loc < pop[ind].sel.size(); loc++)
					{
						index = pop[ind].sel[loc];
						mutations[index][1 + nv] += 1;
					}
				}
			}

			// record fixed, lost and segregating variants and update wildtype
			for (unsigned int loc = 0; loc < mutations.size(); loc++)
			{
				if ((mutations[loc][1 + nv] / fourN) == 1) // fixed variants
				{
					//cout << "mut " << loc << " fixed! al count: " << mutations[loc][1 + nv] << endl;
					fixed.push_back(loc); //add recently fixed mutations to list of fixed mutations
                    //cout << "At generation " << gen << ", we have " << fixed.size() << ", fixed variants, " << ", with indices:";
                    //for (int i=0; i<fixed.size(); i++)
                    //{
                      //cout << fixed[i] << ", ";  
                   // } 
                   // cout << endl;
					if( fixtimes[loc] == 0)
						fixtimes[loc] = gen; //record generation of fixation
					fout << mutations[loc][0] << "\t" << muttimes[loc] << "\t" << fixtimes[loc] << "\t" << optval << "\t"; //report fixation to output file
					for (int k = 0; k < nv; k++)
					{
						wildtype[k] += mutations[loc][k + 1]; //update the wildtype
						fout << mutations[loc][k + 1] << "\t" << delta[loc][k] << "\t";
					}
					fout << endl;
					nb_fix += 1;
				}
				if (mutations[loc][1 + nv] == 0) //lost variants
				{
					if (find(lost.begin(), lost.end(), loc) == lost.end()) //if locus was not yet on the lost list
					{	
                        lost.push_back(loc);
                        //cout << "At generation " << gen << ", number of lost mutations:" << lost.size() << endl;
                    }
				}
				if (double((mutations[loc][1 + nv] / fourN)) > 0 && (mutations[loc][1 + nv] / fourN) < 1) //compute number of segregating variants
					nb_seg += 1;
			}
		}


		if (gen % UPDATELENGTH == 0)    //display progress of simulation
			fprintf(stderr,"After this cycle, at gen %lld, we have fixed %d with %d segregating.\n",gen,nb_fix, nb_seg);

		gen += 1;
		// **** bookkeeping ****

		// stop simulation once D mutations fixed
		if (nb_fix >= Dv)
		{
			maxgen = gen;
			gen = EXITFLAG;
		}


		// **** end bookkeeping ****


	} // end of recursion for loop
	//   ---------------------------   END RECURSION    ---------------------------

	fout.close(); //close output file

	//record duration of simulation
	finish = time(0);
	int temps = int(difftime(finish, start));
	cout << maxgen << " generations have taken " << temps << " seconds\n";

	delete [] temp;
	delete [] Wtot;
	delete [] pop;

}
