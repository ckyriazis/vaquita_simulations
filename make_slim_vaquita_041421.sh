# Script to make SLIM job script


# Set K1, the ancestral carrying capacity for the first epoch
K1=${1}

# Set K2, the ancestral carrying capacity for the second epoch
K2=${2}

# Set bottMin, the minumum number of individuals for the bottleneck
bottMin=${3}

# Set recovMortRate, the stochastic  mortality rate during the recovery
recovMortRate=${4}

# Set burn-in duration
burnIn=${5}

# Set duration of second epoch - leave this hard coded for now
T2=26000


# Make script
cat > vaquita_${K1}K1_${K2}K2_${bottMin}bottMin_${recovMortRate}recovMortRate_${burnIn}burnIn_041421.slim << EOM

initialize() {
	
	initializeSLiMModelType("nonWF");
	initializeSLiMOptions(keepPedigrees=T);
	defineConstant("K1", ${K1}); //ancestral carrying capacity for first epoch
	defineConstant("K2", ${K2}); //ancestral carrying capacity for second epoch	
	defineConstant("sampleSize", 60); //for sampling the pop for summary stats
	defineConstant("count", 0); //counter used for stopping decline
	defineConstant("repr_int", 1); // calving interval - can handle 1 or 2 years
	defineConstant("p_repr", 0.75); // probability of an adult female reproducing when calving interval = 1
	defineConstant("minRepAge", 5); //minimum age for reproduction
	defineConstant("minIndNum", ${bottMin}); //minimum pop size following bottleneck
	defineConstant("bottMortRate1", 0.15); //mortality rate  during initial part of decline 
	defineConstant("bottMortRate2", 0.34); //mortality rate  during second phase of decline
	defineConstant("recovMortRate", ${recovMortRate}); //mortality rate  during recovery
	defineConstant("geneLength", 1760);
	initializeSex("A");
	
	// age-dependent increases in mortality
	defineConstant("L", c(0.4, 0.05, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.1, 0.1, 0.1, 0.1, 0.1, 0.2, 0.2, 0.2, 0.2, 0.2, 1.0));	
	
	initializeMutationRate(5.8e-9);
	defineConstant("h_VstrDel", 0.0);
	defineConstant("h_strDel", 0.01);
	defineConstant("h_modDel", 0.1);
	defineConstant("h_wkDel", 0.4);
	
	//draw deleterious mutations from vaquita DFE
	//and allow for different dominance coefficients for mutations with different s
	//by creating different mutation types (faster than fitness callbacks)	
	//s for moderately and weakly deleterious mutations is divided by the generation time (11.9)
	//to convert selection coefficients from per generation to per year
	

	//very strongly deleterious mutations (s<-0.1)
	initializeMutationType("m1", h_VstrDel, "s", "do x=rgamma(1,-0.025722485,0.1314); while (x >= -0.1); return x;");
	//strongly deleterious mutations (s<-0.01)
	initializeMutationType("m2", h_strDel, "s", "do x=rgamma(1,-0.025722485,0.1314); while (x < -0.1 | x >= -0.01); return x;");
	//moderately deleterious mutations (-0.001 > s >= -0.01)
	initializeMutationType("m3", h_modDel, "s", "do x=rgamma(1,-0.0021615533,0.1314); while (x < -0.01 | x >= -0.001); return x;");
	//weakly deleterious mutations (s >= -0.001)
	initializeMutationType("m4", h_wkDel, "s", "do x=rgamma(1,-0.0021615533,0.1314); while (x < -0.001); return x;");
	//lethal mutations - optionally add some fraction of mutations with s=-1
	//proportion set by initializeGenomicElementType() call below (default is 0.0)
	initializeMutationType("m5", 0.0,"f",-1.0);
	//neutral mutations
	initializeMutationType("m6", 0.5,"f",0.0);
	
	//ratio of different deleterious mutation types taken from vaquita DFE (sum to 100 below)
	//assume ratio of deleterious to neutral muts of 2.31:1 (Huber et al 2017) 
	//giving 100/2.31=43.3 for neutral mutations below
	initializeGenomicElementType("g1", c(m1,m2,m3,m4,m5,m6), c(7.7,20.7,18.4,53.2,0.0,43.3));	


	//number of genes on each autosome from vaquita annotations
	gene_vec=c(1860,1319,1454,789,687,890,655,1131,631,1011,1082,433,609,673,1192,642,436,335,1075,1084,185);
	
	defineConstant("seqLength", sum(gene_vec)*geneLength);
	
	gene_num=sum(gene_vec);
	
	for (i in 1:gene_num){
		initializeGenomicElement(g1, ((i-1)*geneLength)+(i-1), (i*geneLength)+(i-2) );
	}
	
	rates=NULL;
	
	//assume no recombination within genes, a rate of 1e-3 between genes, and free recombination between chroms
	for (i in 1:(size(gene_vec)-1)){
		rates=c(rates, 0, rep(c(1e-3, 0), asInteger(gene_vec[i-1]-1)), 0.5);
	}
	rates=c(rates, 0, rep(c(1e-3, 0), asInteger(gene_vec[size(gene_vec)-1]-1)));
	
	ends=NULL;
	for (i in 1:gene_num){
		ends=c(ends, (i*geneLength)+(i-2), (i*geneLength)+(i-1));
	}
	ends=ends[0:(size(ends)-2)];
	
	initializeRecombinationRate(rates, ends);

}


1:$((${burnIn} + ${T2} + 500)) reproduction() {
	
	males = p1.individuals[p1.individuals.sex=='M'];
	females = p1.individuals[p1.individuals.sex=='F'];
	
	//get reproductive age males and females
	repro_females = females[females.age >= minRepAge];
	repro_males = males[males.age >= minRepAge];
	
	//loop through females and reproduce
	//with one randomly selected male (1 offspring)
	for(mom in repro_females){
		
		// get fitness of mom and use to weight probability of reproduction
		// i.e., if fitness is 0.95, reproduction will be successful 95% of time
		// but if its 0.9 it will be successful 90% of time	
		mom_fitness = p1.cachedFitness(mom.index)/p1.individuals.tagF[mom.index];
		if(runif(1)>mom_fitness){
			next;
		}

		if(repro_males.size() > 0){
			
			// for calving interval of 1 or 1.5
			if(repr_int == 1){
				
				// allow all repr females to mate when repr_int == 1
				//probability of actually reproducing determined by p_repr
				if(runif(1)<p_repr){
					dad = sample(repro_males,1);
					child = p1.addCrossed(mom, dad);
					//use tag to keep track of F_ped for each individual
					//multiply by 100000 to make sure its an int
					//since tag only alows for ints, and tagF is used elsewhere
					//divide out 100000 below when outputting mean F_ped
					child.tag = asInteger(mom.relatedness(dad)/2*100000);
				}
			}
			
			// only allow half of repr females to mate when repr_int ==2
			// select by even or odd index every other year
			if(repr_int == 2){
				
				if(mom.index % 2 == sim.generation % 2){
					dad = sample(repro_males,1);
					child = p1.addCrossed(mom, dad);
					//use tag to keep track of F_ped for each individual as described above
					child.tag = asInteger(mom.relatedness(dad)/2*100000);
				}
			}
		}
	}
	self.active = 0;

}



1 early() {
	cat("gen,popSize,numReprMale,numReprFem,numKilled,meanFitness,meanHet,B,FROH_500kb,FROH_1Mb,F_ped,avgLethal,avgVstrDel,avgStrDel,avgModDel,avgWkDel" + "\n");
	sim.addSubpop("p1", K1);
	p1.individuals.age = rdunif(K1, min=0, max=25);
	
	//use sim.tag to enforce sims to run for 50 years after reaching min pop size
	sim.tag = 50;

	//uncomment these if simulating neutral mutations
        //m3.convertToSubstitution = T;
	//m4.convertToSubstitution = T;
	m6.convertToSubstitution = T;
}



1:$((${burnIn} + ${T2} + 500)) early() {
	
	//use life table to enforce mortality after 25 years
	inds = p1.individuals;
	ages = inds.age;
	mortality = L[ages];
	survival = 1 - mortality;
	inds.fitnessScaling = survival;
	
	p_death = recovMortRate;

	
	//after burn in, initiate mortality due to fishing
	if(sim.generation > $((${burnIn} + ${T2}))  & count <1){
		
		//initial phase of decline - lower mortality rate
		if(p1.individualCount > 250){
			p_death = bottMortRate1;
		}
		
		//second phase of decline - increased mortality rate
		if(p1.individualCount <= 250){
			p_death = bottMortRate2;
		}
	
	}
	
	// use p1.tag to keep track of number of bycatch mortalities
	p1.tag = 0;
	
	//apply p_death to each individual as bernoulli trial
	if(sim.generation > $((${burnIn} + ${T2}))){
		for(i in inds){
			
			// kill off individuals with probability p_death
			if(runif(1) < p_death){
				i.fitnessScaling = 0.0;
				p1.tag = p1.tag+ 1;
			}
		}

		//kill off calves whose mother has died	
		calves = p1.individuals[p1.individuals.age==0];
		
		for(calf in calves){
			parent_IDs = calf.pedigreeParentIDs;
			//catn(parent_IDs);
			for(ID in parent_IDs){
				parent = p1.individuals[p1.individuals.pedigreeID==ID];				
				if(parent.sex == "F"){
					mom = parent;
				}
				if(mom.fitnessScaling==0.0){
					calf.fitnessScaling=0.0;
				}
			}
		}
	}


	
        if(sim.generation <= ${burnIn}){
                K = K1;
        }
        if(sim.generation >${burnIn}){
                K = K2;
        }

        //no fitness increases due to density dependence allowed
        p1.fitnessScaling = min(K /(p1.individualCount * mean(survival)), 1.0);
	
	//use p1.individuals.tagF to keep track of population-level fitness scaling
	//as well as fitness scaling for each individual due to age
	//need quantity to divide out of cachedFitness to get unscaled absolute fitness for the population
	p1.individuals.tagF = p1.individuals.fitnessScaling*p1.fitnessScaling;

}



//track statistics pre-bottleneck every 1000 generations
1:$((${burnIn} + ${T2})) late() {
	if (sim.generation % 1000 == 0) {
		stats = getStats(p1, sampleSize);
		
		males = p1.individuals[p1.individuals.sex=='M'];
		females = p1.individuals[p1.individuals.sex=='F'];
		old_females = females[females.age >= minRepAge];
		old_males = males[males.age >= minRepAge];
		cat(sim.generation + "," + p1.individuals.size() + "," + old_females.size() + "," + old_males.size() + "," + p1.tag + "," + stats + "\n");
	
	}

}




// track statistics after decline every generation and terminate when the population goes to 1 individual 
$((${burnIn} + ${T2} + 1)):$((${burnIn} + ${T2} + 500)) late() {

	//once pop size declines to minIndNum individuals, flip 'switch' to reduce mortality rates to recovMortRate
	if(p1.individualCount <=minIndNum){
		rm("count", removeConstants=T);
		defineConstant("count",1);
	}

	if(p1.individuals.size() < 2){
		stats = c("NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA"); //cant get stats from just one individual
	}
	
	// case when p1 size is less than sample size but greater than 1
	if(p1.individuals.size() < sampleSize & p1.individuals.size() > 1){	
		stats = getStats(p1, p1.individuals.size());
	}
	//case when p1 size is greater than or equal to sample size
	if(p1.individuals.size() >= sampleSize){ 
		stats = getStats(p1, sampleSize);
	}
	
	males = p1.individuals[p1.individuals.sex=='M'];
	females = p1.individuals[p1.individuals.sex=='F'];
	old_females = females[females.age >= minRepAge];
	old_males = males[males.age >= minRepAge];
	
	
	cat(sim.generation + "," + p1.individuals.size() + "," + old_females.size() + "," + old_males.size() + "," + p1.tag + "," + stats + "\n");
	
	
	//end sim if pop goes extinct
	if(p1.individuals.size() < 2){
		sim.simulationFinished();
		cat("The population has gone extinct");
	}
	
	
	// count down from 50 years
	if(count==1){
		sim.tag = sim.tag - 1;
	}
	
	
	// end sim when reach 1
	if(sim.tag == 0){
		sim.simulationFinished();
	}
}



// define function to sample a population for
// mean fitness, heterozygosity, inbreeding load (2B), mean Froh, mean Fped 
// and avg num of mutations of different classes per individual (very str del, str del, mod del, wk del)
function (s) getStats(o pop, i sampSize)
{
	i = sample(pop.individuals, sampSize, F);
	
	m = sortBy(i.genomes.mutations, "position"); //get all mutations in sample
	m_uniq = unique(m); // get rid of redundant muts
	DAF = sapply(m_uniq, "sum(m == applyValue);"); // count number of each mut in pop
	m_uniq_polym = m_uniq[DAF != i.genomes.size()]; //remove fixed mutations
	
	//calculate mean pop heterozygosity
	meanHet = calcHeterozygosity(pop.genomes);  
	
	//initialize vectors used for individual statistics
	ROH_length_sumPerInd_1Mb = c();
	ROH_length_sumPerInd_500Kb = c();
	Num_lethal_muts = c();	
	Num_vStrDel_muts = c();
	Num_strDel_muts = c();
	Num_modDel_muts = c();
	Num_wkDel_muts = c();
	F_ped = c();
	B_pop = c();
	
	for (individual in i) {
		
		indm = sortBy(individual.genomes.mutations, "position");
		indm = indm[match(indm, m_uniq_polym) >= 0];   // Check that individual mutations are not fixed 
		indm_uniq = unique(indm);
		
		genotype = sapply(indm_uniq, "sum(indm == applyValue);");
		
		// tally number of mutations for different classes of selection coefficient per individual
		s = individual.genomes.mutations.selectionCoeff;
		
		Num_lethal_muts = c(Num_lethal_muts, sum(s<=-0.5));
		Num_vStrDel_muts = c(Num_vStrDel_muts, sum(s<=-0.1));
		Num_strDel_muts = c(Num_strDel_muts, sum(s<=-0.01));
		Num_modDel_muts = c(Num_modDel_muts, sum(s<=-0.001 & s > -0.01));
		Num_wkDel_muts = c(Num_wkDel_muts, sum(s<=-0.00001 & s > -0.001));
		
		//code for getting ROHs
		ID_het = (genotype == 1); //outputs T/F for genotypes if they are het or homDer
		ID_homDer = (genotype == 2);
		pos_het = indm_uniq.position[ID_het]; //outputs positions of heterozgoys genotypes
		
		startpos = c(0, pos_het); //adds 0 to beggining of vector of hets
		endpos = c(pos_het, sim.chromosome.lastPosition); //adds last position in genome to vector of hets
		pos_het_diff = endpos - startpos;
		
		//sum for ROHs > 1Mb
		ROH_startpos_1Mb = startpos[pos_het_diff > 1000000]; //filter out startpos that dont correspond to ROH > 1Mb
		ROH_endpos_1Mb = endpos[pos_het_diff > 1000000];
		ROH_length_1Mb = pos_het_diff[pos_het_diff > 1000000]; //vector of ROHs for each individual	
		ROH_length_sum_1Mb = sum(ROH_length_1Mb);
		ROH_length_sumPerInd_1Mb = c(ROH_length_sumPerInd_1Mb, ROH_length_sum_1Mb); // add sum of ROHs for each individual to vector of ROHs for all individuals
		
		//sum for ROHs > 500kb
		ROH_startpos_500Kb = startpos[pos_het_diff > 500000]; //filter out startpos that dont correspond to ROH > 1Mb
		ROH_endpos_500Kb = endpos[pos_het_diff > 500000];
		ROH_length_500Kb = pos_het_diff[pos_het_diff > 500000]; //vector of ROHs for each individual	
		ROH_length_sum_500Kb = sum(ROH_length_500Kb);
		ROH_length_sumPerInd_500Kb = c(ROH_length_sumPerInd_500Kb, ROH_length_sum_500Kb); // add sum of ROHs for each individual to vector of ROHs for all individuals

		//get Fped from individual.tag
		F_ped_ind = individual.tag/100000;
		F_ped = c(F_ped, F_ped_ind);

		//calculate 2B (inbreeding load)
		del_muts = c(individual.genomes.mutationsOfType(m1),individual.genomes.mutationsOfType(m2),individual.genomes.mutationsOfType(m3),individual.genomes.mutationsOfType(m4),individual.genomes.mutationsOfType(m5));
		B_ind = c();

		if (del_muts.length()>0) {
			for(m in del_muts){
				//check if mut is heterozygous
				if(individual.genomes.mutationCountsInGenomes(m)==1){
					
					//protect against case where s < -1 (can happen with gamma DFE)
					s = max(m.selectionCoeff,-1.0);
					//difference in fitness between het and hom is s*(h-1) (1+sh -(1+s))
					B_ind = c(B_ind, s*(m.mutationType.dominanceCoeff-1));
					//catn(m.id + "," +  m.selectionCoeff + "," + m.selectionCoeff*(m.mutationType.dominanceCoeff-1));
				}
			}
			// this is summed rather than multiplied
			//even if fitness is multiplicative (Gao 2015 suggests so)
			B_pop = c(B_pop, sum(B_ind));
		}
		else{
			B_pop = c(B_pop, 0.0);
		}		
	}
	
	
	return(mean(pop.cachedFitness(NULL)/pop.individuals.tagF) + "," + meanHet + "," + mean(B_pop) + "," + mean(ROH_length_sumPerInd_500Kb)/seqLength + "," + mean(ROH_length_sumPerInd_1Mb)/seqLength + "," + mean(F_ped) + "," +  mean(Num_lethal_muts) + "," +  mean(Num_vStrDel_muts) + "," + mean(Num_strDel_muts)+ "," + mean(Num_modDel_muts) + "," + mean(Num_wkDel_muts));
}



EOM
