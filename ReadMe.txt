Author: Matthieu Nadini. 03/24/2020
For questions, please contact matthieu.nadini@gmail.com

###########################################################

Sample code for the article "Epidemic spreading and vaccination strategies in a urban-like environment"
- It is about our two-dimensional agent-based model that encapsulates the typical core-peripheral structure of many urban environments
- It considers the SIR disease process only, with few modifications it can be adapted to the SIS as well


###########################################################################
The codes requires either python 2.x or python 3. and c++11 with the boost library

First, you need to modify the run_code.py and type where your boost library is
Then, in order to run the code, type 
- python run_code.py


###############################################################################
In Input folder 

- Parameters.txt
	It contains all input parameters in the model
	You may change the value of the input parameters as you like.
  
- Meaning_parameters.txt
	It explains the meaning of all input parameters
  
- N_runs.txt
	Number of independent runs the algorithm will run with the same parameters
  
##############################################################################
In Output folder 

All output files have the following name:
- "Lambda_over_mu_"+NumberToString(int(100.*pjump))+"_"+NumberToString(int(100.*noise))+"_"+NumberToString(run_i)+".txt"
	NumberToString(int(100.*pjump)) multiply the input parameter pjump per 100 and convert it into a string
	NumberToString(int(100.*noise)) multiply the input parameter noise per 100 and convert it into a string
	NumberToString(int(run_i)) convert the number of indipendent run into a string
	
	Example with pjump = 0.5, noise = 0.3, and run_i = 3:
	"Lambda_over_mu_50_30_3.txt"
	
All output file contains the following info:
- first column represents the timestamp
- second column represents the fraction of recovered
