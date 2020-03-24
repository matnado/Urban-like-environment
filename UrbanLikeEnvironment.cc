//To run the programm
//g++ -std=c++11 -o SIR SIR.cc -L/home/mat/Desktop/Thesis/boost/ -I/home/mat/Desktop/Thesis/boost/

#include <iostream>
#include <fstream>
#include <sstream>
#include <stdlib.h> 
#include <math.h>
#include <time.h>
#include <random>
#include <boost/random.hpp>
#include <boost/generator_iterator.hpp>
#include <boost/format.hpp>
#include <boost/math/constants/constants.hpp>
#include <vector>
#include <random>
#include <map>
#include <algorithm>
#include <set>
#define PI 3.14159265

using namespace std; 

class Network{
	public:
		int N;						//agents in the network
		int n;						//Number of agents in each location
	    	double D;					//dimension of the square
	    	
	    	double L;					//Number of the locations
	    	double Rcmin;					//Minimum radius locations							
	    	double Rcmax;					//Maximum radius locations
	    	double gamma;					//Exponent distribution of the locations
	    	
	    	vector<double> Rc;				//radius of the locations circle
	    	vector<double> x_c;				//x variable, center of the locations -> according to the algorithm
	    	vector<double> y_c;				//y variable, center of the locations -> according to the algorithm
	    	
	    	vector<vector<double>> agents_assigned_to_c;	//Labels 1: locations ID, Labels 2: agents ID in the location
	    	vector<vector<string>> agents_in_c;		//verify whether a specific node is in the location of not. "yes" node is in the location, "no" it is outside
	    	vector<int> count_agents_in_c;
	    	
	    	double p_in;					//probability of being in the location
	    	int T_steady_state;				//Time steps at the steady state
	    	double eta; 					//white noise
	    	double delta_t; 				//discrete time interval
	    	
	    	vector<double> posx;				//x position of the agents
	    	vector<double> posy;				//y position of the agents
	    	
	    	vector<double> posr; 				//position agents with respect to the center of the location
	    	
	    	vector<double> radius; 				//area of influence, when interactions may come
	    	double omega;					//exponent power law distribution of the radius of interactions
	    	
	    	double rmin; 					//minimum area of influence, when interactions may come
	    	double rmax; 					//maximum area of influence, when interactions may come
	    	vector<double> v; 				//velocities
	    	
	    	double vexp; 					//exploration velocities that agents have when they are outside the location 
	    	double c;					//constant in the exponential distribution --> it is a multiple of vexp!!! Net.c=1./(Net.vexp*Net.c)
	    	vector<double> densities_locations;		//density locations.. the velocity inside the locations is inversely proportional to it
	    	
	    	vector<double> theta; 				//angle that determines the node movement
	    	
	    	double pjump;
	    	double noise;
	    	
	    	vector<vector<double>> temporal_interactions; 	//temporal interactions
	    	vector<vector<double>> strengthin;		//individual strength inside the location
	    	
	    	double Lambda; 					//infection rate
	    	double mu; 					//recovery rate 
	    	double perc_Infected; 				//percentage of initial infected agents
	    	
	    	map<string, int> possible_states;		//in SIR: S=0; I=1; R=2
		vector<int> state;				//each node has one of the possible states
	    	
		set<int> infected_agents;			//infected agents
};



long seed;

void srand_un(void) {
	seed=(long)time(NULL);
}

//draw locations' radius from the power law distribution
double draw_radius_agents_from_power_law (Network &Net, const double &y){
	double Ermin = pow(Net.rmin, 1.-Net.omega);
	double Ermax = pow(Net.rmax, 1.-Net.omega);
	
	return pow(y*(Ermax-Ermin) + Ermin, 1./(1.-Net.omega));
}

//draw locations' radius from the power law distribution
double draw_radius_locations_from_power_law (Network &Net, const double &y){
	double Ermin = pow(Net.Rcmin, 1.-Net.gamma);
	double Ermax = pow(Net.Rcmax, 1.-Net.gamma);
	
	return pow(y*(Ermax-Ermin) + Ermin, 1./(1.-Net.gamma));
}

//read important parameters
void read_parameters(Network &Net, const string &input_path, const string &Name_file){
	string line;
	ifstream myfile;
	myfile.open(input_path+Name_file);
	
	if (myfile.is_open()){
		while (true){
   			getline(myfile, line);
   			getline(myfile, line);
   			istringstream ss(line);
 			ss >> Net.N >> Net.D >> Net.L >> Net.Rcmin >> Net.Rcmax >> Net.rmin >> Net.rmax >> Net.omega >> Net.gamma >> Net.vexp >> Net.Lambda >> Net.mu >> Net.perc_Infected >> Net.c >> Net.T_steady_state >> Net.pjump >> Net.noise >> Net.delta_t;
 			break;
		}
		myfile.close();    
	}
	else
		cout<<"file not opened"<<endl;
	
	Net.n = Net.N/Net.L;
	if (int(Net.n)!=Net.n)
		cout<<"N/L should be integer"<<endl;
	
	Net.eta = 360.; //the angle is 360 in total
	Net.c=1./(Net.vexp*Net.c);// 1./sqrt(Net.D/2.-Net.Rc);//2./(Net.D/2.-Net.Rc);//
	Net.p_in = 1./(1.+ (Net.pjump)/(1.-exp(Net.c)));
}

void initialize_states_SIR(Network &Net){
	Net.possible_states["S"] = 0;
	Net.possible_states["I"] = 1;
	Net.possible_states["R"] = 2;
}

//this function check whether 
int check_overlapping(vector<double> &Rc, vector<double> &x_c, vector<double> &y_c, double &suppx, double &suppy, double &Rlast){
	
	for (int j=0; j<x_c.size(); j++){ 
		if (sqrt(pow(x_c[j]-suppx, 2.)+ pow(y_c[j]-suppy, 2.)) < Rc[j]+Rlast){
			return 1;
		}
	}
	return 0;
}

//create the core-peripheral structure
void place_locations_in_space(Network &Net, boost::uniform_01<boost::mt19937> &zeroone){
	double R_av = 0.;
	for (int i=0; i<Net.L; i++){
		R_av += Net.Rc[i];
	}
	R_av/=Net.L;
		
	double Ring_min = Net.Rc[0]+Net.Rcmin, Ring_max = Ring_min+R_av; 
	double angle, pos, suppx, suppy, supp, is_empty, i, p_center_x, p_center_y;
	
	p_center_x = Net.x_c[0];
	p_center_y = Net.y_c[0];
	is_empty = 1;
	i = 1;
	
	while(i< Net.L){
		angle = (zeroone()*Net.eta)*PI/180.;
		pos = zeroone()*(Ring_max-Ring_min);
		suppx = p_center_x+(pos+Ring_min)*cos(angle);
		suppy = p_center_y+(pos+Ring_min)*sin(angle);
		supp = check_overlapping(Net.Rc, Net.x_c, Net.y_c, suppx, suppy, Net.Rc[i]);
	
		while(supp == 1){
			angle = (zeroone()*Net.eta)*PI/180.;
			pos = zeroone()*(Ring_max-Ring_min);
			suppx = p_center_x+(pos+Ring_min)*cos(angle);
			suppy = p_center_y+(pos+Ring_min)*sin(angle);
			supp = check_overlapping(Net.Rc, Net.x_c, Net.y_c, suppx, suppy, Net.Rc[i]);
			is_empty +=1;
			if (is_empty > 99){
				is_empty = 1;
				Ring_min = Ring_max;//+Rmin;
				Ring_max = Ring_min+R_av;
			}
		}
		is_empty = 1;
		Net.x_c.push_back(suppx);
		Net.y_c.push_back(suppy);
		i++;
	}	
}

//Set initial parameters, create the urban-like enviroment, and initialize the system itself
void initialize_network(Network &Net, const string &input_path, const  string &Name_file, const string &output_path, boost::uniform_01<boost::mt19937> &zeroone){
	read_parameters(Net, input_path, Name_file);
	
	initialize_states_SIR(Net);
	Net.x_c.clear();
	Net.y_c.clear();
	if (Net.L>0){//draw dimension of the locations from the distibution
		Net.agents_assigned_to_c.resize(Net.L), Net.densities_locations.resize(Net.L), Net.agents_in_c.resize(Net.L), Net.strengthin.resize(Net.L), Net.Rc.resize(Net.L), Net.count_agents_in_c.resize(Net.L);
		int j = 0;
		for (int i=0; i<Net.L; i++){
			if (Net.Rcmin!=Net.Rcmax){
				Net.Rc[i]= draw_radius_locations_from_power_law(Net, zeroone());
			}
			else{
				Net.Rc[i]=Net.Rcmin;
			}
		}
		
		sort (Net.Rc.begin(), Net.Rc.end());
		
		for (int i=0; i<Net.L; i++){
			Net.densities_locations[i] = double(Net.n)/(PI*Net.Rc[i]*Net.Rc[i]);
			for (int k=0; k<Net.n; k++){
				Net.agents_assigned_to_c[i].push_back(j);
				j++;
			}
		}
	}
	else{
		Net.agents_assigned_to_c.resize(Net.L+1), Net.densities_locations.resize(Net.L+1), Net.agents_in_c.resize(Net.L+1), Net.strengthin.resize(Net.L+1), Net.Rc.resize(Net.L+1), Net.count_agents_in_c.resize(Net.L+1);
		Net.Rc[0]=0.;
		Net.pjump=0.;
		Net.densities_locations[0]=0.;
		Net.noise = 0.;//if no locations are present, agents move at random
		for (int k=0; k<Net.N; k++){
			Net.agents_assigned_to_c[0].push_back(k);
		}
	}
	
	Net.x_c.push_back(Net.D/2.);//the location is placed at the center of the map for simplicity
	Net.y_c.push_back(Net.D/2.);
	
	if (Net.L>1){
		place_locations_in_space(Net, zeroone);
	}
	Net.posx.resize(Net.N), Net.posy.resize(Net.N), Net.posr.resize(Net.N), Net.theta.resize(Net.N), Net.radius.resize(Net.N), Net.v.resize(Net.N);
	
	int k = 0;
	if (Net.L>0.){
		for (int j=0; j<Net.L; j++){ 
			Net.count_agents_in_c[j] = 0;
			for (int i=0; i<Net.n; i++){
				Net.strengthin[j].push_back(0.);
			
				Net.posx[k] = zeroone()*Net.D;
				Net.posy[k] = zeroone()*Net.D;
				Net.posr[k] = sqrt(pow(Net.posx[k]-Net.x_c[j],2.)+pow(Net.posy[k]-Net.y_c[j],2.));
			
				Net.theta[k] = (zeroone()*Net.eta)*PI/180.;
				if (Net.rmin!=Net.rmax){
					Net.radius[k] = draw_radius_agents_from_power_law(Net, zeroone());
				}
				else
					Net.radius[k] = Net.rmin;
			
				while (Net.posr[k]> Net.Rc[j]){//to make sure agents start inside their assigned location
					Net.posx[k] = (zeroone()*2.*Net.Rc[j])-Net.Rc[j]+Net.x_c[j];
					Net.posy[k] = (zeroone()*2.*Net.Rc[j])-Net.Rc[j]+Net.y_c[j];
					Net.posr[k] = sqrt(pow(Net.posx[k]-Net.x_c[j],2.)+pow(Net.posy[k]-Net.y_c[j],2.));
				}
				
				Net.v[k] = Net.vexp;
			
				if (Net.posr[k]< Net.Rc[j]){
					Net.agents_in_c[j].push_back("yes");
					Net.count_agents_in_c[j]++;		 
				}
				else{
					
					Net.agents_in_c[j].push_back("no");	
				}
				k++;
			}
		}
	}
	else{
		for (int k=0; k<Net.N; k++){
			Net.posx[k] = zeroone()*Net.D;
			Net.posy[k] = zeroone()*Net.D;
			Net.v[k] = Net.vexp;
			Net.agents_in_c[0].push_back("no");
			if (Net.rmin!=Net.rmax)
				Net.radius[k] = draw_radius_agents_from_power_law(Net, zeroone());
			else
				Net.radius[k] = Net.rmin;
		}
	}
}

//reset important parameters related to the network
void reset_network(Network &Net, int &count_I, boost::uniform_01<boost::mt19937> &zeroone){
	Net.state.clear();
	Net.infected_agents.clear();
	
	for (int i=0; i<Net.N; i++) {
		Net.state.push_back(Net.possible_states["S"]);
	}

	int node_infected, sum_initial_infected=0;
	count_I=int(Net.N*Net.perc_Infected);
	for (int i=0; i<count_I; ){
		node_infected = int(zeroone()*Net.N);
		if (Net.state[node_infected]==Net.possible_states["S"]){
			Net.state[node_infected] = Net.possible_states["I"];
			Net.infected_agents.insert(node_infected);
			sum_initial_infected++;
			i++;
		}
	}
	
}

//determines the direction agents should have to go toward their location 
double Determine_angle_to_location(Network &Net, int &j, int &k){
	double angle = atan2(Net.y_c[j]-Net.posy[k],Net.x_c[j]-Net.posx[k]);
	if (angle<0.)
		angle = angle+2.*PI;
	return angle;
}

//function that regulates the agents' motion				
void update_positions_angles(Network &Net, boost::uniform_01<boost::mt19937> &zeroone){
	double suppx, suppy, suppr, suppPhi, deltar;
	double supprjump = 0.;
	int count = 0;
	int k = 0;
	if (Net.L>0.){
		for (int j=0; j<Net.L; j++){
			for (int i=0; i<Net.n; i++){
				if (Net.agents_in_c[j][i] == "yes"){
					if (zeroone()<Net.pjump){//jump outside the location
						Net.posr[k] = -log(1.-zeroone())/Net.c+Net.Rc[j];
						//cout<<k<<endl;
						//cout<<Net.posr[k]<<endl;
						supprjump += (Net.posr[k]-Net.Rc[j]);
						count++;
					
						if (Net.posr[k]>Net.D/2.){//just to make sure the agent does not jump outside the location
							Net.posr[k] = Net.D/2.;
							cout<<"a"<<endl;
						}	
						suppPhi = (zeroone()*Net.eta)*PI/180.;
			
						Net.posx[k] = Net.posr[k]*cos(suppPhi)+Net.x_c[j];
						Net.posy[k] = Net.posr[k]*sin(suppPhi)+Net.y_c[j];
					
						Net.agents_in_c[j][i] = "no";
						Net.count_agents_in_c[j]--;
						Net.theta[k] = Net.noise*(zeroone()*Net.eta)*PI/180.+Determine_angle_to_location(Net, j, k);
					}
					else{//remain inside the location
						//move at random due to fast velocity
						Net.posx[k] = (zeroone()*2.*Net.Rc[j])-Net.Rc[j]+Net.x_c[j];
						Net.posy[k] = (zeroone()*2.*Net.Rc[j])-Net.Rc[j]+Net.y_c[j];
						Net.posr[k] = sqrt(pow(Net.posx[k]-Net.x_c[j],2.)+pow(Net.posy[k]-Net.y_c[j],2.));
						while (Net.posr[k]> Net.Rc[j]){
							Net.posx[k] = (zeroone()*2.*Net.Rc[j])-Net.Rc[j]+Net.x_c[j];
							Net.posy[k] = (zeroone()*2.*Net.Rc[j])-Net.Rc[j]+Net.y_c[j];
						Net.posr[k] = sqrt(pow(Net.posx[k]-Net.x_c[j],2.)+pow(Net.posy[k]-Net.y_c[j],2.));
						}
						Net.theta[k] = (zeroone()*Net.eta)*PI/180.;
					}
				}
				else{//move outside the location
					suppx = Net.posx[k] + Net.v[k]*cos(Net.theta[k])*Net.delta_t;
					suppy = Net.posy[k] + Net.v[k]*sin(Net.theta[k])*Net.delta_t;
					suppr = sqrt(pow(suppx -Net.x_c[j],2.)+pow(suppy-Net.y_c[j],2.));

					if (suppr<= Net.Rc[j]){//ends up back in the location, in that case, it ends in a random position
						
						Net.agents_in_c[j][i] = "yes";
						Net.count_agents_in_c[j]++;
						//just a trial ...this is the problem
						Net.posx[k] = (zeroone()*2.*Net.Rc[j])-Net.Rc[j]+Net.x_c[j];
						Net.posy[k] = (zeroone()*2.*Net.Rc[j])-Net.Rc[j]+Net.y_c[j];
						Net.posr[k] = sqrt(pow(Net.posx[k]-Net.x_c[j],2.)+pow(Net.posy[k]-Net.y_c[j],2.));
						while (Net.posr[k]> Net.Rc[j]){
							Net.posx[k] = (zeroone()*2.*Net.Rc[j])-Net.Rc[j]+Net.x_c[j];
							Net.posy[k] = (zeroone()*2.*Net.Rc[j])-Net.Rc[j]+Net.y_c[j];
							Net.posr[k] = sqrt(pow(Net.posx[k]-Net.x_c[j],2.)+pow(Net.posy[k]-Net.y_c[j],2.));
						}
						Net.theta[k] = (zeroone()*Net.eta)*PI/180.;
					}
					else{
						//periodic boundary conditions
						while (suppx>Net.D)
							suppx -= Net.D;
						while (suppx<0.)
							suppx += Net.D;
				
						while (suppy<0.)
						 	suppy += Net.D;	
						while (suppy>Net.D)
						 	suppy -= Net.D;
					
						Net.posx[k] = suppx;
						Net.posy[k] = suppy;
						Net.posr[k] = sqrt(pow(Net.posx[k]-Net.x_c[j],2.)+pow(Net.posy[k]-Net.y_c[j],2.));
						Net.theta[k] = Net.noise*(zeroone()*Net.eta)*PI/180.+Determine_angle_to_location(Net, j, k);
					}
				}
				k++;
			}
		}
	}
	else{
		for (int k=0; k<Net.N; k++){
			Net.posx[k] = zeroone()*Net.D;
			Net.posy[k] = zeroone()*Net.D;
		}
		count = 1;
	}
}			


//detect all interactions based on proximity
void determine_interactions(Network &Net){
	Net.temporal_interactions.clear();
	Net.temporal_interactions.resize(Net.N);
	double max_radius;
	for (int i=0; i<Net.N; i++){
		for (int j=0; j<Net.N; j++){
			
			if (Net.radius[i]>Net.radius[j])
				max_radius = Net.radius[i];
			else
				max_radius = Net.radius[j];
			
			// interactions in square with periodic boundary conditions
			if (i<j){
				if (abs(sqrt(pow(Net.posx[i]-Net.posx[j],2.)+pow(Net.posy[i]-Net.posy[j],2.)))< max_radius){
					Net.temporal_interactions[i].push_back(j);
					Net.temporal_interactions[j].push_back(i);
				}
				else if (abs(sqrt(pow(Net.posx[i]-Net.D-Net.posx[j],2.)+pow(Net.posy[i]-Net.posy[j],2.)))< max_radius){
					Net.temporal_interactions[i].push_back(j);
					Net.temporal_interactions[j].push_back(i);
				}
				else if (abs(sqrt(pow(Net.posx[i]-Net.posx[j],2.)+pow(Net.posy[i]-Net.D-Net.posy[j],2.)))< max_radius){
					Net.temporal_interactions[i].push_back(j);
					Net.temporal_interactions[j].push_back(i);
				}
				else if (abs(sqrt(pow(Net.posx[i]-Net.D-Net.posx[j],2.)+pow(Net.posy[i]-Net.D-Net.posy[j],2.)))< max_radius){
					Net.temporal_interactions[i].push_back(j);
					Net.temporal_interactions[j].push_back(i);
				}
				else if (abs(sqrt(pow(Net.posx[i]+Net.D-Net.posx[j],2.)+pow(Net.posy[i]-Net.posy[j],2.)))< max_radius){
					Net.temporal_interactions[i].push_back(j);
					Net.temporal_interactions[j].push_back(i);
				}
				else if (abs(sqrt(pow(Net.posx[i]-Net.posx[j],2.)+pow(Net.posy[i]+Net.D-Net.posy[j],2.)))< max_radius){
					Net.temporal_interactions[i].push_back(j);
					Net.temporal_interactions[j].push_back(i);
				}
				else if (abs(sqrt(pow(Net.posx[i]+Net.D-Net.posx[j],2.)+pow(Net.posy[i]+Net.D-Net.posy[j],2.)))< max_radius){
					Net.temporal_interactions[i].push_back(j);
					Net.temporal_interactions[j].push_back(i);
				}
				else if (abs(sqrt(pow(Net.posx[i]-Net.D-Net.posx[j],2.)+pow(Net.posy[i]+Net.D-Net.posy[j],2.)))< max_radius){
					Net.temporal_interactions[i].push_back(j);
					Net.temporal_interactions[j].push_back(i);
				}
				else if (abs(sqrt(pow(Net.posx[i]+Net.D-Net.posx[j],2.)+pow(Net.posy[i]-Net.D-Net.posy[j],2.)))< max_radius){
					Net.temporal_interactions[i].push_back(j);
					Net.temporal_interactions[j].push_back(i);
				}
			}
		}
		
	}	
}	

		

//convert string to number
int StringToNumber(const string &String){
	stringstream convert(String);
	int number;
	convert>>number;
	return number;
}

//convert number to string
string NumberToString (const double &Number){
	stringstream ss;
	ss<<Number;
	return ss.str();
}

//read the file where the number of runs is contained
void read_exe(vector<int> &N_runs, const string &input_path, const string &Name_exe_file){
	string line;
	ifstream myfile;
	int supp;
	
	myfile.open(input_path+Name_exe_file);	
	if (myfile.is_open()){
		getline(myfile, line);//discard first line
		while (true){
			getline(myfile, line);
			if (myfile.eof()) 
				break;
			istringstream ss(line);
			ss >> supp ;
			N_runs.push_back(supp);
		}
		myfile.close();    
	}
	else
		cout<<"file not opened"<<endl;
}

	
//spreading of SIR disease on the urban-like environment	
void spread_disease_SIR(Network &Net, int &count_I, int &count_R, boost::uniform_01<boost::mt19937> &zeroone){
	vector<int> new_infected;
	vector<int> new_recovered;
	int supp;
	set<int>::iterator it;
	
	for (it=Net.infected_agents.begin(); it!=Net.infected_agents.end(); it++){
		//supp += Net.temporal_interactions[*it].size();
		for (int j=0; j<Net.temporal_interactions[*it].size(); j++){
			if (Net.state[Net.temporal_interactions[*it][j]] == Net.possible_states["S"] and zeroone()<Net.Lambda){
				new_infected.push_back(Net.temporal_interactions[*it][j]);
			}
		}

		if (zeroone()<Net.mu)
			new_recovered.push_back(*it);	
	}
	
	for (int i=0; i<new_infected.size(); i++){
		if (Net.state[new_infected[i]] == Net.possible_states["S"]){
			Net.state[new_infected[i]] = Net.possible_states["I"];
			count_I++;
			Net.infected_agents.insert(new_infected[i]);
		}
	}
		
		
	for (int i=0; i<new_recovered.size(); i++){
		Net.infected_agents.erase(new_recovered[i]);
		Net.state[new_recovered[i]] = Net.possible_states["R"];
		count_R++;
		count_I--;
	}	
}



//main algorithm
int main(int argc, char *argv[]){

	//seed random numbers and random numer generation
  	srand_un();
  	boost::mt19937 rng(seed);
  	static boost::uniform_01<boost::mt19937> zeroone(rng);
 		
 	string input_path = "./Data/Input/";
	string output_path = "./Data/Output/";
	string Name_exe_file = "N_runs.txt";
	string Name_file = "Parameters.txt";
	vector<int> N_runs;
	
	read_exe(N_runs, input_path, Name_exe_file);
	string open_file="on";
	
	Network Net; 
	
	ofstream myfile;
	
	
	string Name_output_file;
	vector <int> I_final;
	initialize_network(Net, input_path, Name_file, output_path, zeroone);
	for (int run_i=0; run_i<N_runs[0]; run_i++){//# of repetitions for the same value of the parameter tested
		initialize_network(Net, input_path, Name_file, output_path, zeroone);
		int count_I;
		int count_R = 0;
				
		read_parameters(Net, input_path, Name_file);
		reset_network(Net, count_I, zeroone);
				
		string Name_output_file = "Lambda_over_mu_"+NumberToString(int(100.*Net.pjump))+"_"+NumberToString(int(100.*Net.noise))+"_"+NumberToString(run_i)+".txt";
		myfile.open(output_path+Name_output_file);
		
		if (myfile.is_open()){
			for (int i=0; i<Net.T_steady_state; i++){//to reach the steady state
				update_positions_angles(Net, zeroone);
			}
			int i=0;
			while (count_I>0){
				determine_interactions(Net);
				update_positions_angles(Net, zeroone);
				spread_disease_SIR(Net, count_I, count_R, zeroone);

				myfile<<i+1<<"	"<<float(count_R)/Net.N<<endl;
				//cout<<i+1<<"	"<<float(count_R)/Net.N<<"	"<<float(count_I)/Net.N<<endl;
				i++;
			}
			myfile.close();
		}
		else
			cout<<"file not opened"<<endl;
	}		

	return 0;
}

