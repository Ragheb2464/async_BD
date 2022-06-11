#ifndef CLASS_VF_EQUL_H
#define CLASS_VF_EQUL_H


ILOSTLBEGIN


class class_VF
{ 
  private:  
		int *globallyFixed;//variables that have been fixed globally so far (we dont want to fix a variable twice)
		int *toBeGlobalyFixedZero, *toBeGlobalyFixedOne;//the variables to be fixed globally in the next round
		vector <int> toBeLocallyFixedZero, toBeLocallyFixedOne;//the variables to be fixed locally in the next round
  public:        		
		class_VF(Data_S *data_S, int n_sc)		
		{
			//memory assignment to the arraies
			globallyFixed = new int[data_S->getN_arcs()]();
			toBeGlobalyFixedZero= new int[data_S->getN_arcs()]();
			toBeGlobalyFixedOne= new int[data_S->getN_arcs()]();
		}
//////////////////////////////////////////////////////////////////////////////////////////////////////////
		void GlobalVarFixation(double globalUB, double globalLB, IloNumArray &reducedCosts)
		{
			int counter =0;
		    for(int a =0; a<data_S->getN_arcs(); a++){	
				toBeGlobalyFixedZero[a] =0;//refresh the list of variables to be fixed this time
				toBeGlobalyFixedOne[a]  =0;
		        if(globallyFixed[a] == 0){	//if it is not already fixed      		      
					if(reducedCosts[a] > globalUB - globalLB){//fix to 0	
						toBeGlobalyFixedZero[a] = 0;
						globallyFixed[a] = 1;
						counter++;
					}
					if(reducedCosts[a] <  globalLB - globalUB){//fix to 1
						toBeGlobalyFixedOne[a] = 1;
						globallyFixed[a] = 1;
						counter++;
					}		
		        }
			}
			if(counter>0)
				cout<< "\n number of newly fixed variables globally are: " << counter << endl;
		}
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		void LocalVarFixation(double globalUB, double localLB, IloNumArray &reducedCosts)
		{	
			int counter=0;
			toBeLocallyFixedOne.clear();
			toBeLocallyFixedZero.clear();
		    for(int a =0; a<data_S->getN_arcs(); a++)	
				if(globallyFixed[a] == 0){//if it has not been fixed globally so-far		  			      
					if(reducedCosts[a] > globalUB - localLB){//fix to 0						
						toBeLocallyFixedZero.push_back(a);//[a] =1;
						counter++;
					}
					if(reducedCosts[a] <  localLB - globalUB){//fix to 1
						toBeLocallyFixedOne.push_back(a);//[a]  =1;
						counter++;
					}
				}			
			if(counter>0)
				cout<< "\n number of newly fixed variables locally are: " << counter << endl;
		}
//////////////////////////////////////////////////////////////////////////////////////////////////
		int gettoBeLocallyFixedOneSize()   {return toBeLocallyFixedOne.size();}
		int gettoBeLocallyFixedOne(int a)  {return toBeLocallyFixedOne[a];}
		int gettoBeLocallyFixedZeroSize()  {return toBeLocallyFixedZero.size();}
		int gettoBeLocallyFixedZero(int a) {return toBeLocallyFixedZero[a];}
		int gettoBeGlobalyFixedZero(int a) {return toBeGlobalyFixedZero[a];}
		int gettoBeGlobalyFixedOne(int a)  {return toBeGlobalyFixedOne[a];}
//////////////////////////////////////////////////////////////////////////////////////////////////
		~class_VF()
		{
		}


};//WND class



#endif
