#ifndef CLASS_GSC_H
#define CLASS_GSC_H

//ordering in ascending fashion
bool PairCompare(const std::pair<double, int>& firstElem, const std::pair<double, int>& secondElem) 
{
  return firstElem.first < secondElem.first;
}
#include <ilcplex/ilocplex.h>
#include "./class_data_S.h"
#include "./class_search_param.h"
#include "./Sub_Classes/class_poolsandmasgers.h"
// #include "./Sub_Classes/class_heurstic.h"
#include "./Sub_Classes/Sol_Choice.h"
#include "./Sub_Classes/LB_Lifting.h"
#include "./Sub_Classes/Net_Connectivity.h"
#include "./Sub_Classes/Cover_Ineq.h"
#include "./Sub_Classes/Cardinality_Ineq.h"
#include "./Sub_Classes/Pure_Cardinality.h"




ILOSTLBEGIN
#define Comm COMM_WORLD

    #include "./Sub_Classes/CallBacks.h"


class class_GSC
{ 
  private:    
		//staffs for the master cplex(master model)
		IloEnv env;		//environment 
		IloModel model;
		IloCplex cplex;
		IloObjective obj;		
 		IloNumVarArray3 x;
 		IloNumVarArray  y;
 		IloNumVarArray  theta;		
		IloRangeArray *con, *con_optcut, *con_copy, *y_heur_con, *connectivity_con, *ineq_con;//regular constraints; optimality cuts; copied optimality cuts; y fixation; network connectivity inequalities                
		//object of the required classes
		PoolsandMangers *Managers;
		MasterSolChoice *Sol_Choice;
		Data_S *data_S ; 		        //create the objects of the data class
		Search_Param	*search_param;
		class_PMCI	*raw_cardinality;
// 		MasterHeuristic *masHeuristic;
		//control and output variables 
		int  total_size_of_con_copy;					
		int  Equal_LB, iteration;
		double Best_Lower_Bound, wtime, tEndLS, phase1Time, phase1Obj, LB_best_time, interTime, interLB, interUB;
		bool transition_phase, end_transition_phase, M_D_Phase_I, continue_running;
		bool interPhaseOptimal;
		//the active/inactive counter to be used in master cleaning procedure
		vector <int >  *con_copy_binding, *con_binding, *con_copy_iter, *con_iter;		
		vector <vector <double > > ineqCoeff;//to keep the coefficient of the generated cover and cardinality cuts	  
  public:      
		class_GSC(int n_workers0, int n_sc0, string instName, string scenarioName, string insClass)		
		{
			n_sc 		=n_sc0;
			n_workers 	=n_workers0;
			//create the object of the classes   
			search_param = new Search_Param(n_sc,n_workers);
			  
			data_S = new Data_S(insClass); 		      	  //create the objects of the data class 	
			data_S->info(instName);		      			//here we fill the network information		  			  
			data_S->info_scenario(n_sc, scenarioName);   //this read the scenario file and fill the p array (probabilities)			    	          	 			 
			data_S->read_assign(n_sc);	              	//fill the arrays and vectors	  				
			if( search_param->getData_Clean_UP_Moode() ) //clean up data	    		      			   		 
				data_S->clean_data(n_sc);	       			    
			data_S->clean_memory(n_sc);                    //free some memory in the  data_S
			//reset the actual number of scenarios (sub-problems)		          
			n_sc = n_sc - search_param->getMaster_sc(); // number of scenarios are total scenarios minus global scenarios (those who are not projected)
			//print the information of the instances and the programm
			PrintInstInfo();			  
			//initialize the search paramters
			iteration	=0;	//counter of number of iterations
			M_D_Phase_I 	=false;	//indicating if the first phase of the algorithm is done		  			  			  			 
			arry_y_fixed 	=new int[data_S->getN_arcs()]();
			total_fix_Y	=0;		  
 			con_copy_binding  = new vector< int >();
			con_binding 	  = new vector< int >();
			con_copy_iter     = new vector< int >();
			con_iter 	  = new vector< int >();			  
			//initialize the global lower bound  
			GlobalLowerBound  	= new double[2];
			GlobalLowerBound[0]	= -1e75;	//lower bound
			GlobalLowerBound[1]	= -1e75;	//lower bound at the end of first phase before the transition phase
			//the information gathered from the master problem to be stored in the memory   
			Sol_send  = new double[data_S->getN_arcs() + n_sc  + 2 + 1]; // y sol and obja val, theta sol, Mp solving time ... 1 is for obj of the global scenarios
			//cluster and assign scenarios to the workers   
			Managers  = new PoolsandMangers(n_workers, n_sc, data_S, search_param);			  
			Managers->Cluster_Sce(data_S, n_workers, n_sc + search_param->getMaster_sc(), search_param->getMaster_sc());//here we cluster the scenarios
			Managers->Assign_Sce_Proc(n_workers, search_param);//here we assign each processor (including master) with its scenarios			  		  
			//let slaves start from core point generating cut	
			if(search_param->getWroker_start_from_core_sol())//if we want the workers to initially start from core point
			{
			    for(int a = 0; a < data_S->getN_arcs()- data_S->getN_od(); a++) //initialization of the M-W core point
					Sol_send[a] = search_param->getInitial_core_point();
			    for(int k = data_S->getN_arcs() - data_S->getN_od(); k < data_S->getN_arcs(); k++)		  
					Sol_send[k] = 1;				    
			    Sol_send[data_S->getN_arcs()]= 1;					    
			    for(IloInt s = 0; s < n_sc ; s++)
					Sol_send[data_S->getN_arcs() + 1 + s] = 0;					    
			    Sol_send[data_S->getN_arcs() + n_sc + 1] = 0;	
			    double Recourse_Obj = 0;
			    if(search_param->getMaster_sc() != 0)
			      Recourse_Obj = 1e10;
			    else
			      for(int a = 0; a < data_S->getN_arcs()- data_S->getN_od(); a++)//when we dont have global scenario we can easily compute the recourse cost
					Recourse_Obj += data_S->getF(a) * search_param->getInitial_core_point();
			    Sol_send[data_S->getN_arcs() + n_sc + 1 + 1] = Recourse_Obj;			
			    Managers->Sol_Pool_Manager(data_S, search_param->getMaster_sc(), n_sc, n_workers, Sol_send,  GlobalLowerBound); //add the y variable and the theta solution to the pool and create the associated memory for tracking  v_sc_time,v_w_time
			    Managers->Sol_Pool_controler(n_workers, n_sc , data_S,  search_param);	//sends solutions to the workers "asynchronously" based on the given strategy ... he wont send if there is nothing to send					
			}  		    
			//create the master model  
			MasterModel(); //this creates the master problem
			//cover inequlity class object
			if(search_param->getAddCoverIneq())
				coverCuts = new class_CI(data_S, n_sc);
			//cardinality inequlity class object
			if(search_param->getAddCardinalityIneq())
				cardinalityCuts = new class_MCI(data_S, n_sc);
			//adding intial cardinality cuts by inspecting the ...
			if(search_param->getAddRawCardinalityCuts()){
				raw_cardinality = new class_PMCI(data_S, n_sc);
				addRawCardinalityCuts();
				delete raw_cardinality;
			}				
// 			if(search_param->getRunHeuristic()){
// 				masHeuristic = new MasterHeuristic(ineqCoeff, Managers, data_S, search_param, n_sc);				
// 				//masHeuristic->cardinality_cuts(data_S, search_param, n_sc);
// 			}
// 			Sol_Choice    = new MasterSolChoice(data_S, search_param, n_sc); //this creates the master problem			  
			 
			Equal_LB = ceil(n_sc/search_param->getN_cut_wait())+1;		//to check if the solution from current call of the master problem has been stored in the memory (true) or not (false);; how many equal lower bounds we will see in first phase before ending it
			Best_Lower_Bound =0; 							//to have the best lower bound (and its time)
			stopping_ceriteria =false, transition_phase =false,  end_transition_phase =false;			  
		}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		//the master processor main loop
		void main_loop()
		{		
			//time holder
			wtime = MPI::Wtime();
			//add network connectivity cuts
			if(search_param->getNetwork_Connectivity())
			  connectivity_constraints(!M_D_Phase_I);
			//initialize the paramters which will be used 
			bool cutLoop=true;//to keep solving the master problem until new valid inequalities are added
			bool boolCI1=false, boolCI2=false, boolMCI1=false, boolMCI2=false;//is new conver of size 1,2 is found; is cardinality of size 1,2 is found
			//the main loop of the algorithm
			while( !stopping_ceriteria )				//while stopping criteria are false 
			{
				iteration++;		
				//check for the first phase stopping statisfaction 
			  	if(!M_D_Phase_I){//M-D First Phase							 
				    Managers->first_phase_ending(phase1Time, phase1Obj, wtime, GlobalLowerBound, Equal_LB, transition_phase, M_D_Phase_I, stopping_ceriteria, data_S->getN_arcs(), search_param->getTolerance(), Sol_send);									
				    if ( stopping_ceriteria )//termintation	, because the first phase might be optimal
				      break;
				    if(M_D_Phase_I){//we dont want to receive cuts for the solution discovered in the first phase during the second phase				    
				      Sol_Choice->clean_phase_I_sol( Managers, n_workers, data_S->getN_arcs(), n_sc, search_param->getTag_Useless_Sol());		    
// 				      Managers->CleanY_Sol_Pool();//clear the Y pool
    //  			      Managers->Cleanv_ub();
    //  				  Managers->Cleanv_sc_obj();
				    }
				}
				//the cut loop of the first phase if it is not done
				if(!M_D_Phase_I)//I dont want to still solve the master for one more iteration even if the first phase is done
				{
				    //if we give the master a warm start by starting from the core point
				    if(iteration==1 && search_param->getWroker_start_from_core_sol()){
					Managers->cut_pool_manager( stopping_ceriteria, n_sc, n_workers, GlobalLowerBound, data_S,  search_param); //v_sc_time,v_w_time, return number of cuts that we can receive			 				            				      	
					add_regular_cuts();
  // 				        add_copied_cuts();
				    }				   
				    //we add valid ineq as much as they exist before adding benders cut
					cutLoop=true;
					while(cutLoop){
						//solve the master
						GlobalLowerBound[0] = solve_phase_I();//, v_ac, v_ac_copy, v_inac, v_inac_copy, v_acinac_aux);										  
						//check for strong ineq
						if(search_param->getStrong_Ineq())
						  strong_inequality();	
						//cover inequalities 
						if(search_param->getAddCoverIneq())
							boolCI1 = addViloatedCoverIneq(1);				
						//cover inequalities of size 2
						if(search_param->getAddCoverIneq())
							boolCI2 = addViloatedCoverIneq(2);
						//if no cover exist then check for the cardinalities 
						if(!boolCI2 && !boolCI1){
						  //cardinality inequalities
						  if(search_param->getAddCardinalityIneq())
							  boolMCI1 = addViolatedCardinalityCuts(1); 
						  //cardinality inequalities of size 2
						  if(search_param->getAddCardinalityIneq())
							  boolMCI2 = addViolatedCardinalityCuts(2);
						}
						//if no wnew cut is derived then break the loop 
						if(!boolCI1 && !boolMCI1 && !boolCI2 &&  !boolMCI2)
						  cutLoop=false;						
				    }
				    
				    Managers->Sol_Pool_Manager(data_S, search_param->getMaster_sc(), n_sc, n_workers, Sol_send,  GlobalLowerBound); //add the y variable and the theta solution to the pool and create the associated memory for tracking  v_sc_time,v_w_time
				    					    
				    //I wrote this part to have the time when the best lower bound is found
				    if (GlobalLowerBound[0] > Best_Lower_Bound+0.5)
				    {
				      tEndLS = MPI::Wtime() - wtime;
				      Best_Lower_Bound = GlobalLowerBound[0]; 
				      cout << endl << "the best lower bound is found at the time: " <<  tEndLS << " with value of: " << Best_Lower_Bound << "  @ iteration: " << iteration-1 << endl;
				      LB_best_time = tEndLS;
				    }					    
				    
				    Managers->Sol_Pool_controler(n_workers, n_sc , data_S,  search_param);	//sends solutions to the workers "asynchronously" based on the given strategy ... he wont send if there is nothing to send					
				    Managers->cut_pool_manager( stopping_ceriteria, n_sc, n_workers, GlobalLowerBound, data_S,  search_param); //v_sc_time,v_w_time, return number of cuts that we can receive			 				            				      	
				    
// 				    Sol_Choice->regular_cuts(data_S->getN_arcs(), n_sc, Managers); //add the recived cuts to the master
// 				    Sol_Choice->copied_cuts(Sol_send, n_sc, Managers, data_S,search_param); //add the recived cuts to the master				    
// 				    Sol_Choice->solve_choice_mas(n_workers, n_sc, Managers, data_S, search_param);		    
// 				    
				    add_regular_cuts(); //add the received cuts to the master
// 				    add_copied_cuts(); //add the received cuts to the master
					/* if(search_param->getRunHeuristic())
						solve_inter_phase(); */					
				}
				else//M-D Second Phase
				{						
					if(search_param->getNetwork_Connectivity())
						connectivity_constraints(!M_D_Phase_I); 	
					//here we get the reduced costs we need, 
					if(search_param->getReducedCostFixing()){
						ReducedCosts(); //to get the reduced costs at the root node
					}
					//clean master problem
					if(search_param->getCleanMaster())
					    clean_master();//before branching we remove those cuts that have larger slack values
					 //impose the integrality requirement 
					 ConvertLPtoMIP();	//here we turn into integer those who are not	
      // 				 heuristic();	//we attempt to inject feasible solution for the second phase based on the solutions in the first phase
					interPhaseOptimal=false;
					if(search_param->getRunHeuristic()){
						(cplex).setParam(IloCplex::TiLim, min(3600.0, 36000 - phase1Time));
						solve_inter_phase();
					}  
					//look for the reduced cost based cuts
					if(search_param->getReducedCostFixing()){
						addReducedCostBasedIneq(); 
					} 
					if(!interPhaseOptimal && (36000 - (MPI::Wtime() - wtime) > 1) && (((Managers->getGlobalUpperBound(0)-phase1Obj)/Managers->getGlobalUpperBound(0)) > 0.003) ){//start second phase if any time is left
						(cplex).setParam(IloCplex::TiLim, 36000 - (MPI::Wtime() - wtime));
						solve_phase_II(); 
					}
// 				    cplex.exportModel("model.lp");			
 				    break; 
				}								
			} //end while
			
			PrintResults();  
			StopWorkers();		
		}		
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		void solve_inter_phase()
		{
			//start time of inter phase
			double interStart = MPI::Wtime() - wtime;
			double m_objval;
			//get the current pool size
			int poolSize = Managers->getY_Sol_Poolsize();
			//fix the variables
			int counter =0;
			IloExpr exprInf(env);//if the fixation yielded infeasible problem, then at least one of those fixed to zero has to be opened
			for(int a =0; a < data_S->getN_arcs(); a++){
				if(y_SOL[a] <= 0.05){
				   (*y_heur_con).add(IloRange(env, 0, y[a], 0));
				   exprInf += y[a];
				   counter++;
				}
				if(y_SOL[a] >= 0.95){
				   (*y_heur_con).add(IloRange(env, 1, y[a], 1));
				   counter++;
				}
			}
 			(model).add((*y_heur_con));
			cout << endl << "num fixed var: " << counter << " out of " << data_S->getN_arcs() << endl ;
			//solve the restricted master model	
			for(int itr =0; itr <5; itr++)
			    if(solve_fix_opt_phase()){
				  (cplex).getValues(setSolVal, setSolVar);
				  m_objval = cplex.getObjValue();
				  IloExpr expr(env);
				  cout <<endl;
				  for(int a=0; a< data_S->getN_arcs(); a++){
				    if(setSolVal[a] >1e-3)
				      cout << a << "; ";
				    
				    if(/*y_SOL[a] > 0.05 &&*/ setSolVal[a] <= 0.01 )
					      expr 	+= y[a];							
				    if(/*y_SOL[a] < 0.95  &&*/ setSolVal[a] >= 0.99 )
					      expr 	+= 1 - y[a];														
				  }
				  cout <<endl;
				  (/*(*y_heur_con)*/model).add(IloRange(env, 1, expr, IloInfinity));//y[a].setUB(0);
// 				  model.add((*y_heur_con));
				  expr.end();
				  
				  
				  
				  //get the values we need(plus those for IloIncumbentCallBack)
				  double Recourse_Obj=0;										    	
				  int t = data_S->getN_arcs() + n_sc;
				  for(IloInt a = 0; a < data_S->getN_arcs(); a++)	
					for(IloInt k = 0; k < data_S->getN_od(); k++)
					      for(IloInt s = 0; s < search_param->getMaster_sc(); s++)
						      Recourse_Obj +=  data_S->getP(Managers->getglobal_scenario_id(s)) * data_S->getC(a) * setSolVal[t++];//x_Sol[a][k][s];						  
				  //send the information which is required by workers	
					  
				 //send the information which is required by workers		
				  t=0;			
				  for(IloInt a = 0; a < data_S->getN_arcs(); a++)
					      Sol_send[a] = setSolVal[t++];//y_SOL[a];					    
				  Sol_send[data_S->getN_arcs()]= m_objval; 					    
				  for(IloInt s = 0; s < n_sc ; s++)
					      Sol_send[data_S->getN_arcs() + 1 + s] = setSolVal[t++];//theta_SOL[s];					    
				  Sol_send[data_S->getN_arcs() + n_sc + 1] = 0;	
				  Sol_send[data_S->getN_arcs() + n_sc + 1 + 1] = Recourse_Obj;
					      
				  Managers->Sol_Pool_Manager(data_S, search_param->getMaster_sc(), n_sc, n_workers, Sol_send,  GlobalLowerBound); //add the y variable and the theta solution to the pool and create the associated memory for tracking  v_sc_time,v_w_time
				  Managers->Sol_Pool_controler(n_workers, n_sc , data_S,  search_param);	//sends solutions to the workers "asynchronously" based on the given strategy ... he wont send if there is nothing to send					
				  Managers->cut_pool_manager( stopping_ceriteria, n_sc, n_workers, GlobalLowerBound, data_S,  search_param); //v_sc_time,v_w_time, return number of cuts that we can receive			 				            				      	
				  
				  if((100*(Managers->getGlobalUpperBound(0)-phase1Obj)/Managers->getGlobalUpperBound(0)) <= 0.4){
				    interPhaseOptimal=true;
				    break;
				  }
				  add_regular_cuts();				
			    }
			    else
			      break;
			 
			 
			if(!interPhaseOptimal){
			    (cplex).setParam(IloCplex::IntSolLim, 1);//we are doing this only because to find an injectable solution to the branching sol we dont need to acctually solve it
			    if(solve_fix_opt_phase())			    
				    (cplex).getValues(setSolVal, setSolVar);			  
// 			    else
// 				    abort();
			    (cplex).setParam(IloCplex::IntSolLim,  1e5);
			    //retain the trustable master problem
			    (model).remove((*y_heur_con));
			    (*y_heur_con).clear(); 
			}
			interUB 	= Managers->getGlobalUpperBound(0);//(cplex).getObjValue();			
			//the total time spend on the inter phase
			interTime	= MPI::Wtime() - wtime - interStart;
			cout << "time of the intermidiary phase: " << interTime << endl;		 
		}
///////////////////////////////////////////////////////
		bool solve_fix_opt_phase()
		{
			    bool fesProb=true;
			    if( !(cplex).solve() ){
			      fesProb=false;
			      cout << "****** Danggerrrrrr the injected solution will be rejected: " << (cplex).getStatus() << endl;
// 			      abort();
			    }
			   
			    return fesProb;
		  }
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		void PrintInstInfo()
		{
		    cout <<"\nNumber of Nodes:     			" << data_S->getN_nodes()  << "\n";
		    cout <<"\nNumber of  Arcs:     			" << data_S->getN_arcs()  << "\n";
		    cout <<"\nNumber of ODs:       			" <<  data_S->getN_od()  << "\n";
		    cout <<"\nNumber of Scenarios: 			" << n_sc + search_param->getMaster_sc()  << "\n";		    		    
		    cout <<"\n-------------------------------------------------------------------\n";
		}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		void PrintResults()
		{		  				
			//////////at the end (after main loop) we print some info*****************
			cout << "\n ----------------------------------------Outputs---------------------------------------\n";	
			cout << endl << "the time spend on the relaxation (root): " << phase1Time << endl;
			cout << endl << "the objective of the relaxation (root): " 	<< phase1Obj << endl;
			
			cout << endl << "the time spend on the intermediary phase: " << interTime << endl;
			cout << endl << "The lower bound after intermediary phase: " << interLB << endl;
			cout << endl << "The upper bound after intermediary phase: " << interUB << endl; 
			
			cout << endl <<  "The global lower bound: " << GlobalLowerBound[0] << endl;			
			cout << endl << "The global upper bound: " 	<< Managers->getGlobalUpperBound(0) << endl; 
			cout << endl << "Optimality gap in % is: " 	<< 100*((Managers->getGlobalUpperBound(0) - GlobalLowerBound[0])/Managers->getGlobalUpperBound(0)) << endl;
			
			cout << "\nThe whole run time is " << MPI::Wtime() - wtime << " seconds and numbre of Iterations: " << Managers->getY_Sol_Poolsize()  << endl;
			cout << "\n -------------------------------------------------------------------------------------\n";

			//print the constructed network
			int solID = Managers->getUbImprovSolID(Managers->getUbImprovSolIDSize()-1);
			 cout << endl << "***The opened arcs in the network: ";
			for(int a =0; a <data_S->getN_arcs(); a++)
			  if( Managers->getY_Sol_Pool(solID, a)> 1e-1)
			    cout << a  << "; ";
 		    cout << endl;  
				
		}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		void StopWorkers()
		{
		   	////we tell working processors to finish 
			continue_running = false;
			for(int w_id = 0; w_id < n_workers; w_id++)
			  MPI::Comm.Isend(&continue_running, 1, MPI::BOOL, w_id + 1, search_param->getTagconti_run()); 				      									      
		}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		void MasterModel() //constructor
		{		  
			model 				= IloModel(env);
			cplex 				= IloCplex(model);
			con   				= new IloRangeArray(env);
			con_optcut  		= new IloRangeArray(env);
			con_copy   			= new IloRangeArray(env);
			y_heur_con			= new IloRangeArray(env);
			connectivity_con 	= new IloRangeArray(env);
			ineq_con			= new IloRangeArray(env);
			interPhaseCuts 		= IloRangeArray(env);
			obj   				= IloMinimize(env);
			branchConst			= IloConstraintArray(env) ;
			
			//arraies to querry the solutions
			y_SOL = IloNumArray(env, data_S->getN_arcs());
			theta_SOL = IloNumArray (env, n_sc + search_param->getMaster_sc());
			x_Sol = IloNumArray3(env, data_S->getN_arcs());
			for(IloInt a = 0; a < data_S->getN_arcs(); a++){
			     x_Sol[a] = IloNumArray2(env, data_S->getN_od());
			     for(IloInt k = 0; k < data_S->getN_od(); k++)
					x_Sol[a][k] = IloNumArray(env, search_param->getMaster_sc());			      
			}
			
			reducedCosts = IloNumArray(env, data_S->getN_arcs());
			//***the master only includes y and theta variables****
			y = IloNumVarArray(env, data_S->getN_arcs(), 0, 1);	//ILOINT defining an variable for each arc (on dimensional variable)						
			theta = IloNumVarArray(env, n_sc + search_param->getMaster_sc(), 0, IloInfinity);	//defining variables to approximate the recourse costs
											//I generate a theta for the global scenarios too, but keep their coeficient to zero (just to ease the programming pain :P)
			(model).add(y); // we dont need to add variables explicitly to the model as they will implicitly be added to the model through the constraints
			(model).add(theta);
			//adding x variable to the master for the global scenarios 
			if(search_param->getMaster_sc() > 0){
			    x = IloNumVarArray3(env, data_S->getN_arcs());//second-stage (x) variable (three dimensional cuz each arc (i-j) is represented with a single name "a"
			    for(int a = 0; a < data_S->getN_arcs(); a++){
				    x[a] = IloNumVarArray2(env, data_S->getN_od());
				    for(int k = 0; k < data_S->getN_od(); k++)				
					    x[a][k] = IloNumVarArray(env, search_param->getMaster_sc(), 0, IloInfinity); 									
			    }
			    for(int a = 0; a < data_S->getN_arcs(); a++)
			      for(int k = 0; k < data_S->getN_od(); k++)
				(model).add(x[a][k]);
			}
			//add all variables to an array to ease the getValue and injection of solutions  
			setSolVar = IloNumVarArray( env);//, data_S->getN_arcs() + n_sc+ data_S->getN_arcs()*data_S->getN_od()*search_param->getMaster_sc() );
			setSolVal = IloNumArray(env);//, data_S->getN_arcs() + n_sc + data_S->getN_arcs()*data_S->getN_od()*search_param->getMaster_sc());
			for(IloInt a =0; a<data_S->getN_arcs(); a++)
			  setSolVar.add(y[a]);
			for(IloInt s =0; s< n_sc; s++)
			  setSolVar.add(theta[s]);
			for(IloInt a=0; a<data_S->getN_arcs(); a++)
			  for(IloInt k=0; k< data_S->getN_od(); k++)
			    for(IloInt s=0; s<search_param->getMaster_sc();s++)
			      setSolVar.add(x[a][k][s]);
			// ***** OBJECTIVE FUNCTION****			
			int ss;
			IloExpr expr(env);
			expr.clear();
			for(int a = 0; a < data_S->getN_arcs() - data_S->getN_od(); a++){
			  expr += fabs(data_S->getF(a)) * y[a];			  
			  for(int k = 0; k < data_S->getN_od(); k++)			//the second stage objective for the global scenarios 
			      for(int s = 0; s < search_param->getMaster_sc(); s++){			     
					  ss = Managers->getglobal_scenario_id(s);//global_scenario_id[s];//this gives the indicator of the scenario 
					  expr += fabs(data_S->getP(ss) * data_S->getC(a)) * x[a][k][s];//(1/master_sc) * c[a] * x[a][k][s];//
			      }
			}
			for(int s = 0; s < n_sc + search_param->getMaster_sc() ; s++)			  
			      expr += fabs(data_S->getP(s))* theta[s];//initially I need all theta not be in obj otherwise the master problem runns outbounded
			
			obj.setExpr(expr);
			(model).add(obj);	
			expr.end();
			// Flow conservation CONSTRAINT 
			for(int s = 0; s < search_param->getMaster_sc(); s++)		// flow conservation constraints for global scenarios
			{
				ss = Managers->getglobal_scenario_id(s);//global_scenario_id[s];
				for(int k = 0; k < data_S->getN_od(); k++)				  
					for(int i = 0; i < data_S->getN_nodes(); i++) 
					{ 
						  IloExpr expr(env);//expr.clear();
						  IloExpr expr1(env);
						  for(int a =0; a < data_S->getN_arcs()- data_S->getN_od(); a++)  
						  {
							  if(data_S->getArcs(a,0) == i) 
								  expr += x[a][k][s];								  
							  if(data_S->getArcs(a,1) == i) 
								  expr1 -=  x[a][k][s];
						  }
						  if (i == data_S->getOd(k,0))//if origin
						    (*con).add(IloRange(env, data_S->getD(k,ss), expr+expr1, data_S->getD(k,ss)));//for the origin only outgoing part is neccessary						
						 else
						   if (i == data_S->getOd(k,1))//if destination
							(*con).add(IloRange(env, -data_S->getD(k,ss), expr+expr1, -data_S->getD(k,ss)));						
						   else
							(*con).add(IloRange(env , 0, expr + expr1, 0));						     
						  expr.end();
						  expr1.end();
					}
			}						
			// Capcity CONSTRAINT 2
			for(int s = 0; s < search_param->getMaster_sc(); s++)		//capacity constraints for global scenarios 	
			{
				ss =  Managers->getglobal_scenario_id(s);//global_scenario_id[s];
				for(int a = 0; a < data_S->getN_arcs()- data_S->getN_od(); a++)
				{
					IloExpr expr(env);					
					for(int k = 0; k < data_S->getN_od(); k++)
						expr += x[a][k][s];					
					expr += -data_S->getU(a)*y[a];
					(*con).add(IloRange(env,-IloInfinity,expr,0));
					expr.end();
				}
			}																		
						
			//valid inequalities	
			RBF_Ineq(); //for the recourse cost bounding functions 
			if(search_param->getLooping_arcs())
				looping_arcs();
			if(search_param->getLestCapIneq())
				min_total_capacity();
			
			
			
			(model).add((*con));
				
			(cplex).setParam(IloCplex::RandomSeed, 2015);//fix cplex random seed
			(cplex).setParam(IloCplex::TiLim, 36000);
			(cplex).setParam(IloCplex::ClockType, 1);//to set it equal to cpu time since it is seq
			(cplex).setParam(IloCplex::SimDisplay, 0);			 					
			(cplex).setParam(IloCplex::EpGap, 0.003);
			(cplex).setParam(IloCplex::AdvInd, 0);
//  			(cplex).setParam(IloCplex::NumericalEmphasis, 1);

  			cplex.setParam(IloCplex::HeurFreq, -1);//to turn off the heuristic				    
			/*cplex.setParam(IloCplex::Probe, -1);
			cplex.setParam(IloCplex::Cliques,       -1);
			cplex.setParam(IloCplex::Covers,        -1);
			cplex.setParam(IloCplex::DisjCuts,      -1);
			cplex.setParam(IloCplex::FlowCovers,    -1);
			cplex.setParam(IloCplex::FlowPaths,     -1);
			cplex.setParam(IloCplex::FracCuts,      -1);
			cplex.setParam(IloCplex::GUBCovers,     -1);
			cplex.setParam(IloCplex::ImplBd,        -1);
			cplex.setParam(IloCplex::MCFCuts,       -1);
			cplex.setParam(IloCplex::MIRCuts,       -1);
			cplex.setParam(IloCplex::ZeroHalfCuts,  -1); 
      // 	cplex.setParam(IloCplex::IntParam(2132), -1); 
							  */			      
//         	(cplex).setParam(IloCplex::NodeSel, 1);//node selection strategy 
//          	(cplex).setParam(IloCplex::VarSel, 3);//to branch on variable with max infeasibility
      	    (cplex).setParam(IloCplex::Param::Threads, 1); 
			    
		}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////		
		bool addViloatedCoverIneq(int i)
		{
			bool newCI=false;//if any new CI is added, true
			//generate the cuts
			coverCuts->MainLoop(data_S, i, y_SOL, n_sc + search_param->getMaster_sc());
			//add the cuts to the model
			vector <double > aux;
			for(int q=0; q< coverCuts->getNumCoverIneq(); q++){
				aux.clear();
				aux = coverCuts->getCoverIneqCoeff(q);
				ineqCoeff.push_back(aux);//to keep the inequalities in the memory for the post optimization
				IloExpr expr(env);
				for(int a=0; a< aux.size()-1; a++)//note, the last one in the aux is the right hand side
					expr += aux[a] * y[a];
				(model).add(IloRange(env, aux[aux.size()-1], expr ,IloInfinity ));
				expr.end();
				newCI = true;
			}
			/* if(newCI){
			 (model).add(*ineq_con);
			 // (*ineq_con).clear();
			} */
			return newCI;
		}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////				
		bool addViolatedCardinalityCuts(int i)
		{
			bool newMCI =false;
			//generate the cuts
			cardinalityCuts->MainLoop(data_S, i, y_SOL, n_sc + search_param->getMaster_sc());
			//add the cuts to the model
			vector <double > aux;
			for(int q=0; q< cardinalityCuts->getNumCoverIneq(); q++){
				aux.clear();
				aux = cardinalityCuts->getCoverIneqCoeff(q);
				ineqCoeff.push_back(aux);//to keep the inequalities in the memory for the post optimization
				IloExpr expr(env);
				for(int a=0; a< aux.size()-1; a++)//note, the last one in the aux is the right hand side
					expr += aux[a] * y[a];
				(model).add(IloRange(env, aux[aux.size()-1], expr ,IloInfinity ));
				expr.end();
				newMCI =true;
			}			
			/* if(newMCI){
				 (model).add(*ineq_con);
				 (*ineq_con).clear();
			} */
			return newMCI;
		}
///////////////////////////////		
		void addRawCardinalityCuts()
		{
			//generate the cuts
			raw_cardinality->MainLoop(data_S,  n_sc + search_param->getMaster_sc());
			//add the cuts to the model
			vector <double > aux;
			for(int q=0; q< raw_cardinality->getNumCarIneq(); q++){
				aux.clear();
				aux = raw_cardinality->getCarIneqCoeff(q);
				ineqCoeff.push_back(aux);//to keep the inequalities in the memory for the post optimization
				IloExpr expr(env);
				for(int a=0; a< aux.size()-1; a++)//note, the last one in the aux is the right hand side
					expr += aux[a] * y[a];
				(model).add(IloRange(env, aux[aux.size()-1], expr ,IloInfinity ));
				expr.end();
			}			
			// (model).add(*ineq_con);
			// (*ineq_con).clear();
		}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////		
		void looping_arcs()
		{
			for(int a =0; a < data_S->getN_arcs() - data_S->getN_od() ; a++)
				if(data_S->getArcs(a,0) == data_S->getArcs(a,1))
					(*con).add(IloRange(env, 0, y[a] ,0 ));
		}//looping arcs								
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		void RBF_Ineq()
		{				
				 //LBF inequality
				  if(search_param->getLB_Lifting())//if we want to add LB lifting inequalities
				  {
					class_LB_Lifting *LBF;
					LBF = new class_LB_Lifting(n_sc, Sol_Choice, data_S, search_param);
					for(int s=0 ;s <n_sc;s++){
						IloExpr expr(env);						
						expr += LBF->getRotuingCost();
						for(int k=0; k< data_S->getN_od(); k++)
						  expr += (data_S->getD(k,s) - LBF-> getMinD(k)) * LBF-> getDual(k);
						for(int a=0; a < data_S->getN_arcs() - data_S->getN_od(); a++)
						  expr += data_S->getF(a)* (LBF->getY(a) -y[a]);
						  
						(*con).add(IloRange(env, -IloInfinity, expr - theta[s], 0 ));
						expr.end();
					  }					  
					  //lifting y variables based on LBF and recourse upper boudning 
					  if(  search_param->getY_Lifting())
					  {
						  int *coeffi;
						  coeffi = new int[data_S->getN_arcs() - data_S->getN_od()];
						  for(int a =0; a < data_S->getN_arcs() - data_S->getN_od() ; a++) 
							coeffi[a] = ceil(data_S->getC(a) * data_S->getU(a) + data_S->getF(a)) ;						  
						  int gcd = ArrayGCD(coeffi, 0, data_S->getN_arcs() - data_S->getN_od()-1);						  
						  IloExpr expr(env);
						  for(int a =0; a < data_S->getN_arcs() - data_S->getN_od() ; a++)
							expr += (coeffi[a]/gcd) *y[a];
						  (*con).add(IloRange(env, ceil(ceil(LBF->getMax_LBF())/gcd), expr , IloInfinity ));
						  expr.end();
						  delete []	coeffi;
					  }				  
					  delete LBF;					  
				  }				  
				//upper bound on theta
				if(search_param->getTeta_UB()){
				  IloExpr expr2(env);				  
				  for(int a =0; a < data_S->getN_arcs() - data_S->getN_od() ; a++) 
					  expr2 += data_S->getC(a)  * data_S->getU(a) * y[a] ;
				  for(int s = 0; s < n_sc; s++)
					(*con).add(IloRange(env, -IloInfinity, -expr2 + theta[s], 0 ));
				  expr2.end();
				}	
			}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		void connectivity_constraints(bool LP_Phase)
		{			
			//network connectivity inequalities			
			  int NodId;
			  class_Network_Connectivity *NC;
			  NC = new class_Network_Connectivity(n_sc, Sol_Choice, data_S, search_param);
			  //for the S data set and LP phase the equality constratints are used
			if(LP_Phase){
			  cout << endl << NC->getPureOriginNodeSize() << "  " << NC->getPureDestinationNodeSize() << "  " << NC->getPureTransferNodeSize() << endl;
				//single node (other magnitude will follow)
				for(int i=0; i< NC->getPureTransferNodeSize(); i++){
						NodId = NC->getTransferNode(i);
						IloExpr expr(env), expr1(env);
						for(int a =0; a < data_S->getN_arcs() - data_S->getN_od() ; a++){
						        if(data_S->getArcs(a,0) == NodId) 
								expr += data_S->getU(a)*y[a];								  
							if(data_S->getArcs(a,1) == NodId) 
								expr1 -=  data_S->getU(a)*y[a];
						}
						(*connectivity_con).add(IloRange(env, 0, expr - expr1, 0 ));
						expr.end();
						expr1.end();
				 }
				 (model).add(*connectivity_con);
			 }
			  //for the second phase, if it is R data set update the constraints
			  if(data_S->getInsClass()  == "R" && !LP_Phase){
					(model).remove(*connectivity_con);
					(*connectivity_con).clear();
					for(int i=0; i< NC->getPureTransferNodeSize(); i++){
						NodId = NC->getTransferNode(i);
						IloExpr expr(env), expr1(env);
						for(int a =0; a < data_S->getN_arcs() - data_S->getN_od() ; a++){
						    if(data_S->getArcs(a,0) == NodId) 
								expr += data_S->getU(a)*y[a];								  
							if(data_S->getArcs(a,1) == NodId) 
								expr1 -=  data_S->getU(a)*y[a];
						}						
						for(int a =0; a < data_S->getN_arcs() - data_S->getN_od() ; a++){
						    if(data_S->getArcs(a,0) == NodId)
								(*connectivity_con).add(IloRange(env, 0, expr1 -y[a], IloInfinity ));																  
							if(data_S->getArcs(a,1) == NodId) 
								(*connectivity_con).add(IloRange(env, 0, expr - y[a], 0 ));
						}						
						expr.end();
						expr1.end();
				  }
				  (model).add(*connectivity_con);
			  }
			  
			  delete NC;
		}	
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		void min_total_capacity()
		{		  
				double max_demand=0;
				double scenario_sum_demand;
				//get max demand
				for(int s=0; s <n_sc+search_param->getMaster_sc(); s++){
					scenario_sum_demand =0;
					for(int k =0; k <  data_S->getN_od() ; k++)
						scenario_sum_demand += data_S->getD(k,s);
					if(scenario_sum_demand > max_demand)
						max_demand = scenario_sum_demand;					
				}
				//find the gcd of the capacity vector 				
				int *coeffi;
				coeffi = new int[data_S->getN_arcs() - data_S->getN_od()];
				for(int a =0; a < data_S->getN_arcs() - data_S->getN_od() ; a++) 
					coeffi[a] = ceil(data_S->getU(a)) ;						  
				int gcd = ArrayGCD(coeffi, 0, data_S->getN_arcs() - data_S->getN_od()-1);
				//write the lifted inequality using MIRCuts				
				IloExpr expr(env);
				for(int a =0; a < data_S->getN_arcs() - data_S->getN_od() ; a++)
					expr += (coeffi[a]/gcd) * y[a];
				(*con).add(IloRange(env, ceil(ceil(max_demand)/gcd), expr , IloInfinity ));
				expr.end();
				delete [] coeffi;
								
			}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		//to find gcd of an array
		int ArrayGCD(int *coeffi, int first, int last)
		{
		    int a = 0, b = 0, gcd = 0;
		    if(first == last)
		    {
				gcd = coeffi[first];
				return gcd;
		    }
		    else
		    {
			a = ArrayGCD(coeffi,first,(first+last)/2);
			b = ArrayGCD(coeffi,(first+last)/2+1,last);
			if(a < 0)
			{
			    a = -a;
			}
			if(b < 0)
			{
			    b = -b;
			}
			while(a != b)
			{
			    if(a > b)
			    {
					a = a-b;
			    }
			    else
			    {
					b = b-a;
			    }
			}
				gcd = a;
				return gcd;
		    }
		}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		//this function convert the first-stage variable to integer 
		void ConvertLPtoMIP()
		{		     
		    (model).add(IloConversion(env, y, ILOINT));       		      
		    cout << endl <<endl<< "***All design variables are converted into {0,1} integer" << endl<<endl;	   
		}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////		
		bool strong_inequality()
		{
			bool newSI=false;
			int coutn_Ineq = 0;
			double Viol_Str_Ineq;//amoun of the vilation of stron inequality 
			int ss;
			//disaggregated constraints for the set of global scenarios 
			for(int s = 0; s < search_param->getMaster_sc(); s++)	 	
			{
				ss =  Managers->getglobal_scenario_id(s);//global_scenario_id[s];
				for(int k = 0; k < data_S->getN_od(); k++)		
				    for(int a = 0; a < data_S->getN_arcs()- data_S->getN_od(); a++)
						if(x_Sol[a][k][s] > 1.0e-2){
							//add only violated ones; (note, if a constraint is added previously it wont be violated now)
							if (x_Sol[a][k][s] - min(fabs(data_S->getU(a)),fabs(data_S->getD(k,ss))) * y_SOL[a] > 1.0e-2){
								IloExpr expr(env);
								expr += x[a][k][s];						  
								expr -= min(fabs(data_S->getU(a)), fabs(data_S->getD(k,ss))) * y[a];								
								(*con).add(IloRange(env, -IloInfinity, expr, 0));
								coutn_Ineq++;
								expr.end();							
							}
						}				      
			}
			if(coutn_Ineq >0){
				(model).add((*con));
				newSI= true;
				cout << endl << "Number of violated strong inequalities are: " << coutn_Ineq  << endl;
			}
 			
			return newSI;
		}		
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////		
		double solve_phase_I()//, vector<int> *v_ac, vector<int> *v_ac_copy, vector<int> *v_inac, vector<int> *v_inac_copy, vector<int> *v_acinac_aux ) 
		{                             			
			double m_objval, tStart, tEnd;//, *slackvals;
			double Recourse_Obj = 0;			//the second stage obj associated to the global scenarios
			//---------------------solve ------------------------------
			// tStart= (cplex).getTime();	//this on linux gives cpu time not the wall clock
						
			if( !(cplex).solve() ) 
			{			  			
			  env.error()  << " ** Failed to optimize Master **> " <<  (cplex).getCplexStatus() << endl;			 
			  abort();
			  return 1;
			}
			//to clean up the master
			if(search_param->getCleanMaster())
			  clean_master();
  			// cout << endl << "Time spend on master " << (cplex).getTime() - tStart<<endl;	
			//retrive the required values
			m_objval = (cplex).getObjValue();
			(cplex).getValues(setSolVal, setSolVar);
			// (cplex).getValues(y_SOL, y);			
			int ss, t=0;				    
			for(int a = 0; a < data_S->getN_arcs(); a++)  {	
				y_SOL[a] = setSolVal[a];
			    for(int k = 0; k < data_S->getN_od(); k++){
					// x_Sol[a][k] = ;//(cplex).getValues( x_Sol[a][k], x[a][k]);//the second stage objective for the global scenarios 
					for(int s = 0; s < search_param->getMaster_sc(); s++){	
						x_Sol[a][k][s] = setSolVal[data_S->getN_arcs() + n_sc + t++];
						ss = Managers->getglobal_scenario_id(s);//ss = global_scenario_id[s];
						Recourse_Obj +=  data_S->getP(ss) * data_S->getC(a) * x_Sol[a][k][s];
					}				  
				}
			}    
			// (cplex).getValues(theta_SOL , theta);
			for(int s =0; s <n_sc; s++)
			   theta_SOL[s] = setSolVal[data_S->getN_arcs()+s];
			    
			//----------put the information we obtain from solving the master problem into an array
			//to be send to the processor-0 in the main			
			for(int a = 0; a < data_S->getN_arcs(); a++)
			    Sol_send[a] = y_SOL[a];			
			Sol_send[data_S->getN_arcs()]= m_objval;			
			for(int s = 0; s < n_sc ; s++)
			    Sol_send[data_S->getN_arcs() + 1 + s] = theta_SOL[s];			
			Sol_send[data_S->getN_arcs() + n_sc + 1] = tEnd;	
			Sol_send[data_S->getN_arcs() + n_sc + 1 + 1] = Recourse_Obj;
			//--------------------------------------
// 						
			return m_objval; 
		}
/////////////////////////////////////////////////Master Clean up/////////////////////////////////////////////////////////////////////////////////
		void clean_master()		   
		{
			    
			    //regular cuts			    
			     if((*con_optcut).getSize() >0 && !M_D_Phase_I){			       
				IloNumArray  Reg_slack(env);
				(cplex).getSlacks(Reg_slack, *con_optcut);
				for(int i = 0; i <(*con_optcut).getSize(); i++) {//for every cut in the MP					   				     
					if ((*con_binding).size() > i){										      
					      if(Reg_slack[i] <= 1e5)//active					  
						  (*con_binding)[i]++;//increase binding counter
					}
					else{ 					
					      if(Reg_slack[i] <= 1e5)					      
						(*con_binding).push_back(1);					      
					      else					      				      
						(*con_binding).push_back(0);					      
					}
				  }
				Reg_slack.end();
			     }
				
			if(M_D_Phase_I){
			      //define the array 
			      IloRangeArray temp_con(env);
			      //keep the cuts which have been active 
			      for(int i = 0; i < (*con_binding).size(); i++) 
				  if( (*con_binding)[i] >= 1)//which constraint to keep						  					    
					(temp_con).add((*con_optcut)[i]);
			      //safe gaurd, keep those cuts that have been added late 
			      for(int i= (*con_binding).size(); i < (*con_optcut).getSize(); i++ )
				  (temp_con).add((*con_optcut)[i]);
			      			      					
			      //remove the useless cuts
			      if((*con_optcut).getSize() - (temp_con).getSize()>10){
				  cout<<"\n I am cleaning the copied cuts by removing: "<<  (*con_optcut).getSize() - (temp_con).getSize()  << " cuts out of:  " << (*con_optcut).getSize()<< endl;
				  (model).remove((*con_optcut));
				  (*con_optcut).clear();
				  (*con_optcut).add((temp_con));
				  (model).add((*con_optcut));
			      }
			      (temp_con).end();  
			}
				
				
				
				
				
				
				
			/*	
// 				(cplex).getSlacks(Copy_slack, *con_copy );
				env.out() << endl << "Reg slacks " << Reg_slack << endl;
			    			    
				IloRangeArray *temp_con;
				temp_con  = new IloRangeArray(env);
				for(int i = 0; i <(*con_optcut).getSize(); i++) //for every cut in the MP						
				  if(  Reg_slack[i] < search_param->getReg_Slack ()  )//which constraint to keep					  					    
				      (*temp_con).add((*con_optcut)[i]);					  					
				    															
  				cout<<"\n I am cleaning the copied cuts by removing: "<<  (*con_optcut).getSize() - (*temp_con).getSize()  << " cuts out of:  " << (*con_optcut).getSize()<< endl;					
 				(model).remove((*con_optcut));
 				(*con_optcut).clear();
 				(*con_optcut).add((*temp_con));
 				(model).add((*con_optcut));
				delete  temp_con;
				(*temp_con).end();
				Reg_slack.end();
			     }
			     			     
			     //copied cuts
			    if((*con_copy).getSize() >0){			
				env.out() << endl << "Copyslacks " << Copy_slack << endl;			    
			    
				IloRangeArray *temp_con;
				temp_con  = new IloRangeArray(env);
				for(int i = 0; i <(*con_copy).getSize(); i++) //for every cut in the MP						
				  if(  Copy_slack[i] < search_param->getCopy_Slack()  )//which constraint to keep					  					    
				      (*temp_con).add((*con_copy)[i]);					  					

				cout<<"\n I am cleaning the copied cuts by removing: "<<  (*con_copy).getSize() - (*temp_con).getSize()  << " cuts out of:  " << (*con_copy).getSize()<< endl;					
 				(model).remove((*con_copy));
 				(*con_copy).clear();
 				(*con_copy).add((*temp_con));
 				(model).add((*con_copy));
				delete  temp_con;
				(*temp_con).end();
				Copy_slack.end();
			    }
			    */
			    
			    
			   
				  //****** the con_copy***				  
				/*if ((*con_copy).getSize() > 0)
				  for(int i = 0; i <(*con_copy).getSize(); i++) //for every cut in the MP	
				   {
				     
					if ((*con_copy_binding).size() > i)
					{
					      (*con_copy_iter)[i]++;//increase counter of its existence in the master
					      if((cplex).getSlack((*con_copy)[i]) <= 10)//active					  
						  (*con_copy_binding)[i]++;//increase binding counter
// 					      cout << endl << i << "  --> " << (cplex).getSlack((*con_copy)[i]);//(*con_copy_iter)[i] << "   " << (*con_copy_binding)[i] << endl;
					}
					else 
					{
					      (*con_copy_iter).push_back(1);
					      if((cplex).getSlack((*con_copy)[i]) <= 10)					      
						(*con_copy_binding).push_back(1);					      
					      else					      				      
						(*con_copy_binding).push_back(0);
// 					      cout << endl << i << "  --> " << (cplex).getSlack((*con_copy)[i]);
					}
				  }								  				  				  				  

				  //****** the regular constraints***				  
				if ((*con_optcut).getSize() > 0)
				  for(int i = 0; i <(*con_optcut).getSize(); i++) //for every cut in the MP	
				   {
				     
					if ((*con_binding).size() > i)
					{
					      
					      (*con_iter)[i]++;					      
					      if((cplex).getSlack((*con_optcut)[i]) <= 10)//active					  
						  (*con_binding)[i]++;//increase binding counter
// 					      cout << endl << i << "  --> " << (*con_iter)[i] << "   " << (*con_binding)[i] << endl;
					}
					else 
					{
					      (*con_iter).push_back(1);
					      if((cplex).getSlack((*con_optcut)[i]) <= 10)					      
						(*con_binding).push_back(1);					      
					      else					      				      
						(*con_binding).push_back(0);					      
					}
				  }								  				  				  				  


				  //clean up for the copy-con
			     if (   (*con_copy).getSize() >= 1000   )//we remove cuts if they have been at least for two iterations in the master and never being binding 
			     {
					IloRangeArray *temp_con;
					temp_con  = new IloRangeArray(env);
					for(int i = 0; i <(*con_copy).getSize(); i++) //for every cut in the MP						
					  if( ((*con_copy_binding)[i]  >= 1) || ((*con_copy_iter)[i] <= 10 )   )//which constraint to keep					  					    
					      (*temp_con).add((*con_copy)[i]);					  					

				      if ( (*con_copy).getSize() - (*temp_con).getSize()  >= 10*n_sc )//if we have at least 30*n_sc cuts to remove
				      {															
  					cout<<"\n I am cleaning the copied cuts by removing: "<<  (*con_copy).getSize() - (*temp_con).getSize()  << " cuts out of:  " << (*con_copy).getSize()<< endl;
					
 					(model).remove((*con_copy));
 					(*con_copy).clear();
 					(*con_copy).add((*temp_con));
 					(model).add((*con_copy));
					total_size_of_con_copy = (*con_copy).getSize();				      
				      
					(*con_copy_binding).clear();
					(*con_copy_iter).clear();
				      }
				      
				      delete  temp_con;
				      (*temp_con).end();
		
			     }
			     
			       //clean up for the copy-con
			     if ((*con_optcut).getSize() >= 50*n_sc)//assume each time we add n_sc cuts of the regular ones, so we check for removal if 30 iter are passed
			     {
					IloRangeArray *temp_con1;
					temp_con1   = new IloRangeArray(env);
					for(int i = 0; i <(*con_optcut).getSize(); i++) //for every cut in the MP						
					  if( ((*con_binding)[i]  >= 1) || ((*con_iter)[i] <= 30 )  )//which constraint to keep					  					    
					      (*temp_con1).add((*con_optcut)[i]);					  					

				      if ( (*con_optcut).getSize() - (*temp_con1).getSize() >= 2*n_sc )//if we have at least 3*n_sc cuts to remove
				      {															
  					cout<<"\n I am cleaning the regular cuts by removing: "<<  (*con_optcut).getSize() - (*temp_con1).getSize()  << " cuts out of:  " << (*con_optcut).getSize() << endl;
					
 					(model).remove((*con_optcut));
 					(*con_optcut).clear();
 					(*con_optcut).add((*temp_con1));
 					(model).add((*con_optcut));
					
					(*con_iter).clear();
					(*con_binding).clear();
				      }
				      
				      delete  temp_con1;
				      (*temp_con1).end();
		
			     }
			     
			     
			     
			     
			*/     
			
		  }
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////	
		int add_regular_cuts() //"NOTE: we only have optimality cuts
		{		  			  
			  double  g;
			  int sc_id, sol_id, vio_check, Num_added_cuts=0, solId=1e-2;//the last one is to check we are not adding the same CBC cut several times

			for(int i = 0; i < Managers->getv_cpSize(); i++)
			{
				sol_id = Managers->getv_cp(i, 0);//(*v_cp)[i][0];	//solution ID
				sc_id  = Managers->getv_cp(i, 1);//(*v_cp)[i][1];	//this indicates which scenario this cut belongs to
				g      = Managers->getv_cp(i, 2);//(*v_cp)[i][2];
				
				vio_check =0;				
				for(int a =0; a < data_S->getN_arcs(); a++)
					vio_check += Managers->getv_cp(i, a+3)/*(*v_cp)[i][a+3]*/ * Sol_send[a];//3 + n_arcs + n_sc is to get to the initial location in the array in which the G[a] values are stored																		 //CopyCut*n_sc*n_arcs to get to the initial location of the G[a] for the considered scenario
				vio_check += g - Sol_send[data_S->getN_arcs() + 1 + sc_id]; //Sol_send[n_arcs + 1 + CopyCut]give the theta
								  
			    if( vio_check > 1e-2 )//add the egular cut if it violate the current solution
			    {
					Num_added_cuts++;				    
				    IloExpr expr(env);				    
				    for(int a = 0; a < data_S->getN_arcs(); a++)
						expr += Managers->getv_cp(i, a+3)/*(*v_cp)[i][a+3]*/*y[a];
				    
// 				    cout << endl << g << "+";
// 				    for(int a = 0; a < data_S->getN_arcs(); a++)
// 				      if(Managers->getv_cp(i, a+3) != 0)
// 						cout <<Managers->getv_cp(i, a+3) << "*y[" << a << "]  ";
// 				    cout << " <= theta" << endl;
				    
				    expr += g;
				    expr += -theta[sc_id];			       
				    (*con_optcut).add(IloRange( env, -IloInfinity, expr, 0));//, con_name));					      
				    expr.end();
			    }				  
			    //add combinatorial cuts
				if(search_param->getAddCombinatorialCuts() && Managers->getv_cp(i, 3+data_S->getN_arcs())<1 && solId != Managers->getv_cp(i, 0)){//if we want to add CBC, and the solution has been infeasible, and it has not been added previousely					
					solId = Managers->getv_cp(i, 0);					
					if(IntegerSol(solId)){
						cout << endl << "adding CBC" << endl;
						IloExpr expr(env);
						for(IloInt a = 0; a < data_S->getN_arcs(); a++)
							if(Managers->getY_Sol_Pool(solId, a) <= 0.001)
								expr += y[a];
						(model).add(IloRange( env, 1, expr, IloInfinity));							 
						expr.end();
					}
				}
				
			}
			  
			if(Num_added_cuts >0)
			    (model).add((*con_optcut));
			//clean up the cut pool		 			 			
			Managers->Cleanv_cp();
			return Num_added_cuts;
		}			 
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////			  
		bool IntegerSol(int solId)	//to check if a specific solution is integer
		{
			bool solIsInt=true;
			for(IloInt a = 0; a < data_S->getN_arcs(); a++)
				if(Managers->getY_Sol_Pool(solId, a) <= 0.99 && Managers->getY_Sol_Pool(solId, a) >= 0.01){
					solIsInt =false;
					break;
				}
			
			return solIsInt;
		}		
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////			  
		void add_copied_cuts( ) 
		{			
			  double  g,  Vio_Check;
			  int sc_id, sol_id;
			  //********************************here we copy the cut for other scenarios *******************
			  bool add_cut;			  
			  int cut_to_add_id[search_param->getNum_copy_cut_add()];
			  double cut_to_add_violation[search_param->getNum_copy_cut_add()];
			  int id;
			  			      
			  if(search_param->getNum_copy_cut_add() > 0)
				  for(int CopyCut = 0; CopyCut < n_sc; CopyCut++ )
				  {	
					add_cut = true;
					for (int ss = 0; ss < search_param->getMaster_sc(); ss++)	//if CopyCut is not for global scenarios 
					    if (Managers->getglobal_scenario_id(ss)/*global_scenario_id[ss]*/ == CopyCut )
						add_cut = false;

					    if ( add_cut) //i dont want to add cut for the global scenarios and the current scenairo				  
					    {
						for(int j =0; j <search_param->getNum_copy_cut_add(); j++)
						  cut_to_add_violation[j] = -1e5;
						
						for (int i =0; i < Managers->getv_cp_copySize()/*(*v_cp_copy).size()*/; i++)//here we determine which cuts to add
						  {		
						        sol_id = Managers->getv_cp_copy(i,0);//(*v_cp_copy)[i][0];
						        sc_id  = Managers->getv_cp_copy(i,1);//(*v_cp_copy)[i][1];							
							if ( sc_id == CopyCut )
							  break;
							
							//check if it is violated by the current solution 
							Vio_Check = 0;							
							for(int a =0; a < data_S->getN_arcs(); a++)
							      Vio_Check += Managers->getv_cp_copy(i,2 + n_sc + CopyCut*data_S->getN_arcs() + a)/*(*v_cp_copy)[i][ 2 + n_sc + CopyCut*n_arcs + a]*/ * Sol_send[a];//3 + n_arcs + n_sc is to get to the initial location in the array in which the G[a] values are stored
																		//CopyCut*n_sc*n_arcs to get to the initial location of the G[a] for the considered scenario							
							Vio_Check += Managers->getv_cp_copy(i,2 + CopyCut)/*(*v_cp_copy)[i][2 + CopyCut]*/ - Sol_send[data_S->getN_arcs() + 1 + CopyCut]; //Sol_send[n_arcs + 1 + CopyCut]give the theta
														
							if(Vio_Check >= 10)
							{
	  // 						  sort(cut_to_add_violation, cut_to_add_violation + num_violated_cuts);
	  // 						  sort(cut_to_add_id, cut_to_add_id + num_violated_cuts);
							    for(int j =0; j <search_param->getNum_copy_cut_add(); j++)
								if(  cut_to_add_violation[j] <= 0)
								{		
								  cut_to_add_violation[j] =Vio_Check;
								  cut_to_add_id[j] =i;
								  Vio_Check = -1e6;
								  break;
								}
							    for(int j=0; j<search_param->getNum_copy_cut_add(); j++)
								  if(  Vio_Check >= cut_to_add_violation[j] )
								  {
								      cut_to_add_violation[j] =Vio_Check;
								      cut_to_add_id[j] =i;
								      break;
								  }								  								   
							}
							
							/*
 								  if(CopyCut==1 && sc_id==0)
 								  {
								    cout << endl << "master cut: ";
								    cout << endl << (*v_cp_copy)[i][2+ CopyCut] << "  ";
								    for(int a = 0; a < n_arcs; a++)
								      if((*v_cp_copy)[i][2 + n_sc + CopyCut*n_arcs + a] != 0)
									cout << (*v_cp_copy)[i][ 2 + n_sc + CopyCut*n_arcs + a] << "y[" << a << "] " ;
								      cout << "<= theta" << endl;
 								  }*/
 								  
						  }
						  
						  //*************************************
						  total_size_of_con_copy=0;
						  for(int j=0; j<search_param->getNum_copy_cut_add(); j++)
						      if( cut_to_add_violation[j] >= 10)						 
							  {		
							        total_size_of_con_copy++;
								id = cut_to_add_id[j];
	  // 						      cout << endl << cut_to_add_violation[j];
								IloExpr expr(env);
								for(int a = 0; a < data_S->getN_arcs(); a++)
								    expr += Managers->getv_cp_copy(id, 2 + n_sc + CopyCut*data_S->getN_arcs() + a)/*(*v_cp_copy)[id][2 + n_sc + CopyCut*n_arcs + a]*/ * y[a];							      
								expr += Managers->getv_cp_copy(id, 2 + CopyCut);//(*v_cp_copy)[id][2 + CopyCut];								
								expr += -theta[CopyCut];
								
 								(*con_copy).add(IloRange( env, -IloInfinity, expr, 0));//, con_name));
								expr.end();																
							  }						  						  
					  }//if
					//******************************** end of copy the cut for other scenarios *******************
				  }//for
				  
			//(model).remove((*con_copy));
			if(  total_size_of_con_copy > 0)
			{
  			  (model).add((*con_copy));
			  total_size_of_con_copy = (*con_copy).getSize();
 			  cout << endl << "The size of copied Bender cuts are: " <<  (*con_copy).getSize()  << endl;
			}
						 			 			
			Managers->Cleanv_cp_copy();
			
		}
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		void solve_phase_II( )
		{		    
				// Set up the cut callback to be used for separating Benders' cuts
//  			     (cplex).setParam(IloCplex::Param::Preprocessing::Presolve, IloFalse); 
//  			    (cplex).setParam(IloCplex::Param::Threads, 1); 
//  			     (cplex).setParam(IloCplex::Param::MIP::Strategy::Search, IloCplex::Traditional);
       
   			    (cplex).use(BendersLazyCallback(env,x,y,theta, Managers, Sol_Choice, data_S, search_param));//
 			    (cplex).use(MyHeuristic(env,  x, y, theta, Managers, data_S, search_param));
 			    (cplex).use(ControlIncumSol(env, Managers));
			    //(cplex).use(MyBranch(env, y, Managers, data_S));
			    
			    
					   // if(search_param->getComm_strategy() != 0)
 			    // (cplex).use(BendersUserCallback(env,x,y,theta, Managers, Sol_Choice, data_S, search_param));
			    				
//  		   (cplex).use(MySelect(env));
//  		   (cplex).use(mySolve(env,y)); 			   
// 			   (cplex).use(myContinous(env,y)); 			    
//  		   (cplex).use(loggingCallback(env));
			    if( !(cplex).solve() )			    
			      cout << endl << "the problem is infeasible or unbounded ... take care, kiss kiss :D " << endl;			    
				else
					GlobalLowerBound[0] = (cplex).getBestObjValue();
		  }
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		void ReducedCosts()
		{
 		    if((cplex).solve()){
			   (cplex).getReducedCosts(reducedCosts, y);
			   env.out() << "Reduced Costs = " << reducedCosts << endl;	
			}		   
		} 
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		void addReducedCostBasedIneq()
		{
			//keep the list of variables fixed to 0 or 1
			vector <int> fixedZero, fixedOne;
			delete [] arry_y_fixed;
			arry_y_fixed 	=new int[data_S->getN_arcs()]();
			//we first search for singleton set
			int zeroCount=0, oneCount=0;
			for(int a =0; a < data_S->getN_arcs(); a++){
				fixedZero.push_back(0);
				fixedOne.push_back(0);
				if(reducedCosts[a] >= Managers->getGlobalUpperBound(0) - GlobalLowerBound[0]){
					(*y_heur_con).add(IloRange(env, 0, y[a], 0));
					fixedZero[a] = 1;
					arry_y_fixed[a] = 1;
					zeroCount++;
				}
				if(reducedCosts[a] <= -Managers->getGlobalUpperBound(0) + GlobalLowerBound[0]){
					(*y_heur_con).add(IloRange(env, 1, y[a], 1));
					fixedOne[a] = 1;
					arry_y_fixed[a] = 1;
					oneCount++;
				}
			}
			cout << "\n ******Num of fixed vars to zero and one: (" << zeroCount << ", " << oneCount << ")\n" << endl;
			//cardinality of two
			zeroCount=0;
			for(int a =0; a < data_S->getN_arcs(); a++)
				if(fixedZero[a] == 0 && reducedCosts[a]>0)
					for(int aa =0; aa < data_S->getN_arcs(); aa++)
						if(aa != a && fixedZero[aa] == 0 && reducedCosts[aa]>0)
							if(reducedCosts[a] + reducedCosts[aa] > Managers->getGlobalUpperBound(0) - GlobalLowerBound[0]){
								(*y_heur_con).add(IloRange(env, 0, y[a] + y[aa], 1));
								zeroCount++;
							}
			oneCount=0;
			for(int a =0; a < data_S->getN_arcs(); a++)
				if(fixedOne[a] == 0 && reducedCosts[a] <0)
					for(int aa =0; aa < data_S->getN_arcs(); aa++)
						if(aa != a && fixedOne[aa] == 0 && reducedCosts[aa] <0)
							if(reducedCosts[a] + reducedCosts[aa] < -Managers->getGlobalUpperBound(0) + GlobalLowerBound[0]){
								(*y_heur_con).add(IloRange(env, 1, y[a] + y[aa], 2));
								oneCount++;
							}
			cout << "\n ******Num of added Q and P constraints: (" << zeroCount << ", " << oneCount << ")\n" << endl;
			//add the constraints to the model
			model.add(*y_heur_con);
		}
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////		
 		void heuristic()
 		{
//  		     MasterHeuristic *masHeuristic = new MasterHeuristic(ineqCoeff, Managers, data_S, search_param, n_sc);
// 		     masHeuristic->cardinality_cuts(data_S, search_param, n_sc);
// 		     cout << endl << "size of getBest_UB_sol_id is: " << Managers->getUbImprovSolIDSize() << endl;
// 		     
// 		     for(int i=0; i < Managers->getUbImprovSolIDSize(); i++)
// 		     {
// 		      masHeuristic->y_fixing(Managers, data_S, i);
// 		      Sol_send =  masHeuristic->heuristic_solve(Managers, data_S, search_param, n_sc);    
// 		     
// 		      Managers->Sol_Pool_Manager(data_S, search_param->getMaster_sc(), n_sc, n_workers, Sol_send,  GlobalLowerBound); //add the y variable and the theta solution to the pool and create the associated memory for tracking  v_sc_time,v_w_time
// 		      Managers->Sol_Pool_controler(n_workers, n_sc , data_S,  search_param);	//sends solutions to the workers "asynchronously" based on the given strategy ... he wont send if there is nothing to send					
// 		     }
// 		     
// 		      delete masHeuristic;
// 		     
// 		     
// 		     
// 		     
// 		     
// 		      double ave;
// 		       for(IloInt a=0; a<data_S->getN_arcs()-data_S->getN_od(); a++){
// 			 ave=0;
// 			 for(int i=0; i <Managers->getY_Sol_Poolsize(); i++)
// 			   ave += Managers->getY_Sol_Pool(i, a);
// 			 ave /= Managers->getY_Sol_Poolsize();
// 			 if(ave<=0.3)
// 			   (*y_heur_con).add(IloRange(env, 0, y[a], 0));
// 			 if(ave >= 0.7)
// 			   (*y_heur_con).add(IloRange(env, 1, y[a], 1));
// 		      }
// 	      
// 		     /* int sol_id = Managers->getY_Sol_Poolsize()-1;//Managers->getBest_UB_sol_id();//0;
// 		      
// 		      for(IloInt a=0; a<data_S->getN_arcs(); a++){
// 			if(Managers->getY_Sol_Pool(sol_id, a) <= 0.1 )
// 			  (*y_heur_con).add(IloRange(env, 0, y[a], 0));
// 			if(Managers->getY_Sol_Pool(sol_id, a) >= 0.8)
// 			  (*y_heur_con).add(IloRange(env, 1, y[a], 1));
// 		      }*/
// 		      (model).add((*y_heur_con));
// 		      (cplex).solve();
// 		      cout << endl << endl<< "??????????????????:  " << (cplex).getObjValue() << endl;		      
// 		      
// 		         double m_objval = (cplex).getObjValue();
// 			    int ss;				    
// 			    (cplex).getValues(y_SOL, y);
// 			    double Recourse_Obj = 0;
// 			    for(int a = 0; a < data_S->getN_arcs(); a++)  			      
// 			      for(int k = 0; k < data_S->getN_od(); k++){
// 				 (cplex).getValues( x_Sol[a][k], x[a][k]);//the second stage objective for the global scenarios 
// 				for(int s = 0; s < search_param->getMaster_sc(); s++){			     
// 				      ss = Managers->getglobal_scenario_id(s);//ss = global_scenario_id[s];
// 				      Recourse_Obj +=  data_S->getP(ss) * data_S->getC(a) * x_Sol[a][k][s];
// 				  }				  
// 				}			    
// 			    (cplex).getValues(theta_SOL , theta);
// 			    
// 		      (model).remove((*y_heur_con));
// 		      delete y_heur_con;
// 			    
// 			//----------put the information we obtain from solving the master problem into an array
// 			//to be send to the processor-0 in the main			
// 			  for(int a = 0; a < data_S->getN_arcs(); a++)
// 			      Sol_send[a] = y_SOL[a];			
// 			  Sol_send[data_S->getN_arcs()]= m_objval;			
// 			  for(int s = 0; s < n_sc ; s++)
// 			      Sol_send[data_S->getN_arcs() + 1 + s] = theta_SOL[s];			
// 			  Sol_send[data_S->getN_arcs() + n_sc + 1] = 1;	
// 			  Sol_send[data_S->getN_arcs() + n_sc + 1 + 1] = Recourse_Obj;		      
// 		      		      
// 		     //initialize the injected incumbent 
// 		     for(IloInt a=0; a<data_S->getN_arcs(); a++)
// 			      best_y[a] = y_SOL[a];
// 		     for(IloInt s=0; s<n_sc;s++)
// 			       best_theta[s] = theta_SOL[s]; /*getValue(theta[s])*/;			      
// 		     for(IloInt a=0; a<data_S->getN_arcs(); a++)
// 		         for(IloInt k=0; k< data_S->getN_od(); k++)
// 			    for(IloInt s=0; s<search_param->getMaster_sc();s++)
// 			       best_x[a][k][s] = x_Sol[a][k][s];
// 						     		     
// 		     
// 		   Managers->Sol_Pool_Manager(data_S, search_param->getMaster_sc(), n_sc, n_workers, Sol_send,  GlobalLowerBound); //add the y variable and the theta solution to the pool and create the associated memory for tracking  v_sc_time,v_w_time
// 		   Managers->Sol_Pool_controler(n_workers, n_sc , data_S,  search_param);	//sends solutions to the workers "asynchronously" based on the given strategy ... he wont send if there is nothing to send					
// 	      
// // 		  
 		}
		
		~class_GSC()
		{
		}


};//WND class



#endif
