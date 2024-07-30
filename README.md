# A Particle Markov Chain Monte Carlo Approach to Inference in Surface Kinetics

### Intro
In this paper, we propose a PMCMC sampler to infer the kinetics of the adsorption and desorption of carbon monoxide (CO) on a palladium (Pd(111)) crystal. 
Although it is a system of only 2 molecular species, its complications arise from several possible adsorption sites, and incredibly fast adsorption occurence. 
We use fast-scanning temporal infrared spectroscopy (IR) data, from which we obtain a time-series of the area under the spectrum and take that area as the
input to our PMCMC sampler. 

### Constraints and Modeling
The system comes with several physical constraints which we incorporate as mathematical constraints in our proposed modeling and approach. 
We propose a model of 4 regions which naturally occur based on physical system and experimental design.

### Parameters of Interest
From the sampler, we directly infer: <br/> 
$\theta$ - Time series of the coverage of CO on Pd(111) - interpreted as the hidden states in a State-Space Model. <br/> 
$k_1, k_2$ - Adosrption constants in Region 1 and 2. <br/> 
$k_3, k_4$ - Desorption constants in Region 3 and 4. <br/> 
$\sigma_A$ - Observation Noise of the system. <br/> 

Indirectly from the sampler:
$Ea_1, Ea_2, Ea_3, Ea_4$ - Activation energies in each region. <br/> 
$A_1, A_2, A_3, A_4$ - Pre-exponential factors in each region. <br/> 

### Validity of Results
The kinetic parameters $Ea_1, Ea_4$ and $A_4$ have been obtained independently from several experiments in the literature. 
Our results match these closely. Additionally, the relative obtained values ($Ea_2, Ea_3, A_1, A_2, A_3$) reflect physically meaningful results.

## CODE

Note: Unfortunately at this time we do not publish the data we use (we will in the near future). However, we do provide the structure of the code for the readers interested in the implementation. <br/> 

### Main script
main_PMCMC.m - the main script which runs the proposed PMCMC sampler for a user-specified temperature. <br/> 
approach1.m - a script that reproduces the estimation of the Ea and lnA for all regions based on Approach 1 described in the paper. <br/>
approach2.m - a script that reproduces the estimation of the Ea and lnA for all regions based on Approach 2 described in the paper. <br/>


### Functions
monte_carlo/pf_chem.m - the particle filter to estimate the coverages (hidden states) from the data (area). <br/> 
util/beta_random.m - a function that samples from the Beta distribution with a specified mean, and a fixed alpha prior. <br/> 
util/compute_weights.m - a function that computes the weights inside the particle filter. <br/> 


monte_carlo/MH14.m - the metropolis-hastings algorithm for the parameters in region 1 and 4. <br/> 
monte_carlo/MH23.m - the metropolis-hastings algorithm for the parameters in region 2 and 3. <br/> 

util/compute_Ea.m - a function to compute the activation energy and pre-exponential factor for given rate parameter and temperature inputs. <br/>

util/plotting.m - a script for formatting plots for the paper.




 
