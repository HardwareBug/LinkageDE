#include <vector>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <math.h>
#include <time.h>
#include <random>

// Random device
/*
std::random_device RANDOM_DEVICE;
std::mt19937 MT(RANDOM_DEVICE());
*/
std::mt19937 MT((unsigned)time(NULL));

// Parameter of experiment
const unsigned int EXPERIMENT_NUM = 30;
const unsigned int MAXIMUM_EVALUATION_NUM = 1000000;
const unsigned int DEMENTION_NUM = 22;

// Parameter of evaluation function
// 0 Through, 1 Type_I, 2 Type_II, 3 Type_III
const unsigned int EVALUATION_FUNCTION = 1;
const unsigned int EVALUATION_FUNCTION_T = 2;
const unsigned int EVALUATION_FUNCTION_L = 20;
const long double EVALUATION_ALLOWABLE_ERROR = 0.001; // -1 is disable
// Useage of demention
// Type_I  : EVALUATION_FUNCTION_T + EVALUATION_FUNCTION_L
// Type_II : EVALUATION_FUNCTION_T * 2 + EVALUATION_FUNCTION_L
// Type_III: EVALUATION_FUNCTION_T * 4
const long double MINIMUM_VALUE = -2.048;
const long double MAXIMUM_VALUE = 2.047;
const long double OPTIMAL_VALUE = 1;

// DE option
// The recommended value:
// POPULATION_SIZE = 10 * dimension size
// SCALE_FACTOR = 0.5
// CROSSOVER_RATE = 0.1 (if high dependence, 0.9)
unsigned int POPULATION_SIZE = DEMENTION_NUM * 10; // Useage of evaluation number = POPULATION_SIZE
const long double SCALE_FACTOR = 0.5;
// If use default crossover
const long double CROSSOVER_RATE = 0.1;
// If use linkage crossover
const long double CROSSOVER_RATE_LOW = 0.1;
const long double CROSSOVER_RATE_HEIGH = 0.9;

// Linkage option
// 0 LINC-R, 1 LIDI-R
const unsigned int LINKAGE_TYPE = 1;
const unsigned int LINKAGE_POPULATION_SIZE = 71; // 0 is disable
// Useage of evaluation number = (3 * DIMENSION_NUM * (DIMENSION_NUM - 1) / 2 + 1) * LINKAGE_POPULATION_SIZE
// Default evaluation number is 50,000
// size 4 -> 2632
// size 6 -> 1097
// size 8 -> 589
// size 16 -> 139
// size 22 -> 71
// size 24 -> 61
// size 26 -> 52
// size 28 -> 45
// size 32 -> 34
// size 36 -> 27
const long double LINKAGE_ALLOWABLE_ERROR = 1; // LINC only

// Crossover option
// 0 Use crossoverRateTable, 1 Use linkageGroupTable, 2 Use linkageGroupTable and think linkage group size
const unsigned int CROSSOVER_TYPE = 2;

class Element {
public:
	long double value;
	long double min, max;
	
	Element(long double a, long double b) :min(a), max(b){
		setRandomValue();
	}
	void setRandomValue(){
		std::uniform_real_distribution<long double> random(min, max);
		value = random(MT);
	}
	void setValue(long double a){
		if(a < min){
			std::uniform_real_distribution<long double> random_low(min, value);
			value = random_low(MT);
		}else if(max < a){
			std::uniform_real_distribution<long double> random_high(value, max);
			value = random_high(MT);
		}else{
			value = a;
		}
	}
	void addValue(long double a){ setValue(value+a); }
	long double getValue(){ return value; }
};

class Individual {
public:
	long double fitness;
	std::vector<Element> elements;
	
	Individual();
	Individual(unsigned int individualSize) {
		for (unsigned int i = 0; i < individualSize; i++) {
			elements.push_back(Element(MINIMUM_VALUE, MAXIMUM_VALUE));
		}
	}
};

class Population {
public:
	std::vector<Individual> individuals;
	
	Population() {}
	Population(unsigned int populationSize) {
		for (unsigned int i = 0; i < populationSize; i++) {
			individuals.push_back(Individual(DEMENTION_NUM));
		}
	}
	std::string show(){
		//return "";
		//return (std::to_string(getMimFitness()) + "\n");
		
		std::string ans = "";
		for(unsigned int i = 0; i < individuals.size(); i++){
			for(unsigned int j = 0; j < individuals[i].elements.size(); j++){
				ans = ans + std::to_string(individuals[i].elements[j].value) + " ";
			}
			ans = ans + "| " + std::to_string(individuals[i].fitness) + "\n";
		}
		return ans;
		
	}
	long double getMimFitness(){
		std::sort(individuals.begin(), individuals.end(), [](const Individual &lhs, const Individual &rhs) {
			return lhs.fitness < rhs.fitness;
		});
		return individuals[0].fitness;
	}
};

class Evaluation {
public:
	Evaluation(unsigned int funcNum){
		call = { &Evaluation::through, &Evaluation::typeI, &Evaluation::typeII, &Evaluation::typeIII };
		funcPtr = (call[funcNum]);
		count = 0;
		endFlag = false;
	}
	bool getEndFlag(){ return endFlag; }
	unsigned int getCount(){ return count; }
	Population done(Population arg){
		for(unsigned int i = 0; i < arg.individuals.size(); i++){
			arg.individuals[i] = (this->*funcPtr)(arg.individuals[i]);
			check(arg.individuals[i]);
			count++;
		}
		return arg;
	}
	Individual through(Individual arg) { return arg; }
	Individual typeI(Individual arg) {
		std::vector<long double> temp1, temp2;
		
		for(int i = 0; i < arg.elements.size(); i++){
			if(i < EVALUATION_FUNCTION_T){
				temp1.push_back(arg.elements[i].value);
			}else if(i < EVALUATION_FUNCTION_T + EVALUATION_FUNCTION_L){
				temp2.push_back(arg.elements[i].value);
			}
		}
		
		arg.fitness = fRsn(temp1) + fSp(temp2);
		return arg;
	}
	Individual typeII(Individual arg) {
		std::vector<long double> temp1, temp2;
		
		long double sum = 0;
		for(int i = 0; i < arg.elements.size(); i++){
			if(i < EVALUATION_FUNCTION_T * 2 && i%2 == 0){
				temp1.clear();
				temp1.push_back(arg.elements[i].value);
			}else if(i < EVALUATION_FUNCTION_T * 2 && i%2 == 1){
				temp1.push_back(arg.elements[i].value);
				sum += fRsn(temp1);
			}
			else if(i < EVALUATION_FUNCTION_T * 2 + EVALUATION_FUNCTION_L){
				temp2.push_back(arg.elements[i].value);
			}
		}
		
		arg.fitness = sum + fSp(temp2);
		return arg;
	}
	Individual typeIII(Individual arg) {
		std::vector<long double> temp1, temp2;
		
		long double sum = 0;
		for(int i = 0; i < arg.elements.size(); i++){
			if(i < EVALUATION_FUNCTION_T * 2 && i%2 == 0){
				temp1.clear();
				temp1.push_back(arg.elements[i].value);
			}else if(i < EVALUATION_FUNCTION_T * 2 && i%2 == 1){
				temp1.push_back(arg.elements[i].value);
				sum += fRsn(temp1);
			}
			else if(i < EVALUATION_FUNCTION_T * 4 && i%2 == 0){
				temp2.clear();
				temp2.push_back(arg.elements[i].value);
			}else if(i < EVALUATION_FUNCTION_T * 4 && i%2 == 1){
				temp2.push_back(arg.elements[i].value);
				sum += fSp2(temp2);
			}
		}
		
		arg.fitness = sum;
		return arg;
	}
private:
	long double fRsn(std::vector<long double> arg){
		long double ans = 0;
		for(int i = 1; i < arg.size(); i++){
			ans += 100 * std::pow(arg[0] - std::pow(arg[i], 2), 2) + std::pow(arg[i] - OPTIMAL_VALUE, 2);
		}
		return ans;
	}
	long double fSp(std::vector<long double> arg){
		long double ans = 0;
		for(int i = 0; i < arg.size(); i++){
			ans += std::pow(arg[i] - OPTIMAL_VALUE, 2);
		}
		return ans;
	}
	long double fSp2(std::vector<long double> arg){
		return std::pow(fSp(arg), 2);
	}
	void check(Individual arg){
		if(0 <= EVALUATION_ALLOWABLE_ERROR){
			bool flag = true;
			for(int i = 0; i < arg.elements.size(); i++){
				if( !(OPTIMAL_VALUE - EVALUATION_ALLOWABLE_ERROR / 2 <= arg.elements[i].value && arg.elements[i].value <= OPTIMAL_VALUE + EVALUATION_ALLOWABLE_ERROR / 2) ) flag = false;
			}
			if(flag == true) endFlag = true;
		}
	}
	
	bool endFlag;
	unsigned int count;
	typedef Individual(Evaluation::*funcPtr_t)(Individual);
	typedef std::vector<funcPtr_t> funcPtrVec_t;
	funcPtr_t funcPtr;
	funcPtrVec_t call;
};

void mainloop() {
	std::ofstream output_result("result.txt", std::ios::app);
	
	for(unsigned int i = 0; i < EXPERIMENT_NUM; i++){
		std::cout<<"Experiment"<<i<<std::endl;
		Evaluation evaluation(EVALUATION_FUNCTION);
		
		std::vector<std::vector<long double> > crossoverRateTable;
		crossoverRateTable = std::vector<std::vector<long double> >(DEMENTION_NUM, std::vector<long double>(DEMENTION_NUM, CROSSOVER_RATE));
		for(unsigned int j = 0; j < crossoverRateTable.size(); j++){
			crossoverRateTable[j][j] = 1;
		}

		std::vector<int> dementionGroup(DEMENTION_NUM, -1);
		std::vector<std::vector<unsigned int> > linkageGroupTable;

		Population population(POPULATION_SIZE);
		population = evaluation.done(population);
		
		if(0 < LINKAGE_POPULATION_SIZE){
			std::vector<std::vector<bool> > linkage;
			linkage = std::vector<std::vector<bool> >(DEMENTION_NUM, std::vector<bool>(DEMENTION_NUM, false));
			
			for(unsigned int j = 0; j < LINKAGE_POPULATION_SIZE; j++){
				Population linkagePopulation(1);
				linkagePopulation = evaluation.done(linkagePopulation);
				
				for(unsigned int k = 0; k < DEMENTION_NUM; k++){
					for(unsigned int l = k+1; l < DEMENTION_NUM; l++){
						
						Population linkageTempPopulation(3);
						linkageTempPopulation.individuals[0] = linkagePopulation.individuals[0];
						linkageTempPopulation.individuals[1] = linkagePopulation.individuals[0];
						linkageTempPopulation.individuals[2] = linkagePopulation.individuals[0];
						
						linkageTempPopulation.individuals[0].elements[k].setRandomValue();
						linkageTempPopulation.individuals[1].elements[l].setRandomValue();
						linkageTempPopulation.individuals[2].elements[k].setValue(linkageTempPopulation.individuals[0].elements[k].getValue());
						linkageTempPopulation.individuals[2].elements[l].setValue(linkageTempPopulation.individuals[1].elements[l].getValue());
						linkageTempPopulation = evaluation.done(linkageTempPopulation);
						
						long double f00, f10, f01, f11;
						f00 = linkagePopulation.individuals[0].fitness;
						f10 = linkageTempPopulation.individuals[0].fitness;
						f01 = linkageTempPopulation.individuals[1].fitness;
						f11 = linkageTempPopulation.individuals[2].fitness;
						/*
						std::cout<<"linkageCheck("<<k<<","<<l<<")"<<std::endl;
						std::cout<<linkagePopulation.show();
						std::cout<<linkageTempPopulation.show();
						*/
						if(LINKAGE_TYPE == 0){
							// LINC-R
							if(std::abs((f10 - f00) - (f11 - f01)) > LINKAGE_ALLOWABLE_ERROR){
								//std::cout<<"Hit"<<std::endl;
								linkage[k][l] = true;
							}
						}else if(LINKAGE_TYPE == 1){
							// LIDI-R
							int sgn1, sgn2, sgn3, sgn4;
							if(f00 - f10 > 0) sgn1 = 1;
							else if(f00 - f10 == 0) sgn1 = 0;
							else sgn1 = -1;
							
							if(f01 - f11 > 0) sgn2 = 1;
							else if(f01 - f11 == 0) sgn2 = 0;
							else sgn2 = -1;
							
							if(f00 - f01 > 0) sgn3 = 1;
							else if(f00 - f01 == 0) sgn3 = 0;
							else sgn3 = -1;
							
							if(f10 - f11 > 0) sgn4 = 1;
							else if(f10 - f11 == 0) sgn4 = 0;
							else sgn4 = -1;
							
							if(!(sgn1 == sgn2 && sgn3 == sgn4)){
								//std::cout<<"Hit"<<std::endl;
								linkage[k][l] = true;
							}
						}
					}
				}
			}
			
			for(unsigned int j = 0; j < linkage.size(); j++){
				for(unsigned int k = j+1; k < DEMENTION_NUM; k++){
					linkage[k][j] = linkage[j][k];
				}
			}
			
			// Show linkage
			for(unsigned int j = 0; j < linkage.size(); j++){
				for(unsigned int k = 0; k < linkage[j].size(); k++){
					std::cout<<std::fixed<<linkage[j][k]<<" ";
				}
				std::cout<<std::endl;
			}
			std::cout<<std::endl;
			
			unsigned int groupID = -1;
			for(unsigned int j = 0; j < dementionGroup.size(); j++){
				if(dementionGroup[j] != -1) continue;
				else groupID++;
				dementionGroup[j] = groupID;
				
				std::vector<unsigned int> searchNumber;
				searchNumber.push_back(j);
				while(!searchNumber.empty()){
					int target;
					target = searchNumber.back();
					searchNumber.pop_back();
					
					for(unsigned int k = 0; k < dementionGroup.size(); k++){
						if(dementionGroup[k] == -1 && linkage[target][k] == true){
							dementionGroup[k] = groupID;
							searchNumber.push_back(k);
						}
					}
				}
			}
			/*
			// Show dementionGroup
			for(unsigned int j = 0; j < dementionGroup.size(); j++){
				std::cout<<std::fixed<<dementionGroup[j]<<" ";
			}
			std::cout<<std::endl<<std::endl;
			*/
			for(unsigned int j = 0; j < crossoverRateTable.size(); j++){
				int target;
				target = dementionGroup[j];
				for(unsigned int k = 0; k < crossoverRateTable[j].size(); k++){
					if(k == j) continue;
					else if (dementionGroup[k] == target) crossoverRateTable[j][k] = CROSSOVER_RATE_HEIGH;
					else crossoverRateTable[j][k] = CROSSOVER_RATE_LOW;
				}
			}
			
			for(unsigned int j = 0; j < groupID+1; j++){
				std::vector<unsigned int> temp;
				for(unsigned int k = 0; k < dementionGroup.size(); k++){
					if(dementionGroup[k] == j){
						temp.push_back(k);
					}
				}
				linkageGroupTable.push_back(temp);
			}
		}else{
			for(unsigned int j = 0; j < DEMENTION_NUM; j++){
				std::vector<unsigned int> temp(1, j);
				linkageGroupTable.push_back(temp);
			}
		}
		/*
		// Show crossoverRateTable
		for(unsigned int j = 0; j < crossoverRateTable.size(); j++){
			for(unsigned int k = 0; k < crossoverRateTable[j].size(); k++){
				std::cout<<std::fixed<<crossoverRateTable[j][k]<<" ";
			}
			std::cout<<std::endl;
		}
		std::cout<<std::endl;
		*/
		// Show linkageGroupTable
		for(unsigned int j = 0; j < linkageGroupTable.size(); j++){
			for(unsigned int k = 0; k < linkageGroupTable[j].size(); k++){
				std::cout<<linkageGroupTable[j][k]<<" ";
			}
			std::cout<<std::endl;
		}
		std::cout<<std::endl;
		
		//std::cout<<population.show();
		
		bool endFlag = false;
		while (evaluation.getCount() < MAXIMUM_EVALUATION_NUM) {
			Population tempPopulation = population;
			
			// Crossover
			for(unsigned int j = 0; j < POPULATION_SIZE; j++){
				std::vector<unsigned int> temp(POPULATION_SIZE);
				for(unsigned int k = 0; k < temp.size(); k++){
					temp[k] = k;
				}
				temp.erase(temp.begin() + j);
				std::shuffle(temp.begin(), temp.end(), MT);
				
				if(CROSSOVER_TYPE == 0){
					std::uniform_real_distribution<long double> randomInt(0, DEMENTION_NUM-1);
					unsigned int targetElement = randomInt(MT);
				
					for(unsigned int k = 0; k < DEMENTION_NUM; k++){
						std::uniform_real_distribution<long double> randomReal(0, 1);
						if(randomReal(MT) <= crossoverRateTable[targetElement][k]){
						
							long double base, num1, num2, setValue;
							base = population.individuals[temp[0]].elements[k].getValue();
							num1 = population.individuals[temp[1]].elements[k].getValue();
							num2 = population.individuals[temp[2]].elements[k].getValue();
						
							setValue = base + SCALE_FACTOR * (num1 - num2);
							tempPopulation.individuals[j].elements[k].setValue(setValue);
						}
					}
				}else if(CROSSOVER_TYPE == 1){
					std::uniform_real_distribution<long double> randomReal(0, 1);
					std::uniform_real_distribution<long double> randomInt(0, linkageGroupTable.size()-1);
					unsigned int targetLinkageGroup = randomInt(MT);
					
					for(unsigned int k = 0; k < linkageGroupTable.size(); k++){
						if(randomReal(MT) <= CROSSOVER_RATE_LOW || targetLinkageGroup == k){
							std::uniform_real_distribution<long double> randomInt_2(0, linkageGroupTable[k].size()-1);
							unsigned int targetErement = randomInt_2(MT);
							
							for(unsigned int l = 0; l < linkageGroupTable[k].size(); l++){
								if(randomReal(MT) <= CROSSOVER_RATE_HEIGH || targetErement == l){
									
									unsigned int element = linkageGroupTable[k][l];
									
									long double base, num1, num2, setValue;
									base = population.individuals[temp[0]].elements[element].getValue();
									num1 = population.individuals[temp[1]].elements[element].getValue();
									num2 = population.individuals[temp[2]].elements[element].getValue();
									
									setValue = base + SCALE_FACTOR * (num1 - num2);
									tempPopulation.individuals[j].elements[element].setValue(setValue);
								}
							}
						}
					}
				}else if(CROSSOVER_TYPE == 2){
					std::uniform_real_distribution<long double> randomReal(0, 1);
					std::uniform_real_distribution<long double> randomInt(0, DEMENTION_NUM-1);
					unsigned int target = randomInt(MT);
					unsigned int targetLinkageGroup = dementionGroup[target];
					
					for(unsigned int k = 0; k < linkageGroupTable.size(); k++){
						if(randomReal(MT) <= CROSSOVER_RATE_LOW || targetLinkageGroup == k){
							std::uniform_real_distribution<long double> randomInt_2(0, linkageGroupTable[k].size()-1);
							unsigned int targetErement = randomInt_2(MT);
							
							for(unsigned int l = 0; l < linkageGroupTable[k].size(); l++){
								if(randomReal(MT) <= CROSSOVER_RATE_HEIGH || targetErement == l){
									
									unsigned int element = linkageGroupTable[k][l];
									
									long double base, num1, num2, setValue;
									base = population.individuals[temp[0]].elements[element].getValue();
									num1 = population.individuals[temp[1]].elements[element].getValue();
									num2 = population.individuals[temp[2]].elements[element].getValue();
									
									setValue = base + SCALE_FACTOR * (num1 - num2);
									tempPopulation.individuals[j].elements[element].setValue(setValue);
								}
							}
						}
					}
				}
			}
			
			// Evaluation
			tempPopulation = evaluation.done(tempPopulation);
			if(evaluation.getEndFlag() == true) break;
			
			// Selection
			for(unsigned int j = 0; j < POPULATION_SIZE; j++){
				if(tempPopulation.individuals[j].fitness <= population.individuals[j].fitness){
					population.individuals[j] = tempPopulation.individuals[j];
				}
			}
			
			//std::cout<<population.show();
		}
		if(evaluation.getEndFlag() == true) output_result<<evaluation.getCount();
		else output_result<<"";
		output_result<<",";
		
		std::cout<<std::endl;;
	}
	output_result<<"\n";
	
	output_result.close();
}

int main(){
	for(unsigned int i = 0; i < 7; i++){
		if(i == 0) POPULATION_SIZE = 25;
		else if(i == 1) POPULATION_SIZE = 50;
		else if(i == 2) POPULATION_SIZE = 100;
		else if(i == 3) POPULATION_SIZE = 150;
		else if(i == 4) POPULATION_SIZE = 200;
		else if(i == 5) POPULATION_SIZE = 250;
		else if(i == 6) POPULATION_SIZE = 300;
		mainloop();
	}
}