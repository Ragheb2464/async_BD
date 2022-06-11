#ifndef CLASS_MP_H
#define CLASS_MP_H

//ordering in ascending fashion
bool PairCompare(const std::pair<double, int>& firstElem, const std::pair<double, int>& secondElem) 
{
  return firstElem.first < secondElem.first;
}
#include <ilcplex/ilocplex.h>
//#include "./class_data_S.h"
#include "./LB_Lifting.h"
#include "./Net_Connectivity.h"
#include "./Cover_Ineq.h"
#include "./Cardinality_Ineq.h"
#include "./Pure_Cardinality.h"
#include "./OptRound.h"



ILOSTLBEGIN
#define Comm COMM_WORLD

    #include "./CallBacks.h"


class class_MP
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
		IloRangeArray *con, *con_optcut, *con_copy, *y_heur_con, *connectivity_con, *ineq_con, *artificialCon, *heurYPoolCon;//regular constraints; optimality cuts; copied optimality cuts; y fixation; network connectivity inequalities                
		//object of the required classes
		
		//Data_S *data_S ; 		        //create the objects of the data class		
		class_PMCI	*raw_cardinality;
// 		MasterHeuristic *masHeuristic;
		//control and output variables 
		int  total_size_of_con_copy, Sol_send_size;
		bool  continue_running;
		bool interPhaseOptimal;
		//the active/inactive counter to be used in master cleaning procedure
		vector <int >  *con_copy_binding, *con_binding, *con_copy_iter, *con_iter;		
		vector <vector <double > > ineqCoeff;//to keep the coefficient of the generated cover and cardinality cuts	  
		int numCopyCutUsed;
		vector <int> cutIDsToRemove;
  public:      
		class_MP(Data_S *data_S, MasterSolChoice *Sol_Choice, Search_Param	*search_param, PoolsandMangers *Managers, int n_workers0, int n_sc0, string instName, string scenarioName, string insClass)		
		{
			n_sc 		=n_sc0;
			n_workers 	=n_workers0;
			numCopyCutUsed=0;
			/*//create the object of the classes   			  
			data_S = new Data_S(insClass); 		      	  //create the objects of the data class 	
			data_S->info(instName);		      			//here we fill the network information		  			  
			data_S->info_scenario(n_sc, scenarioName);   //this read the scenario file and fill the p array (probabilities)			    	          	 			 
			data_S->read_assign(n_sc);	              	//fill the arrays and vectors	  				
			if( search_param->getData_Clean_UP_Moode() ) //clean up data	    		      			   		 
				data_S->clean_data(n_sc);	       			    
			data_S->clean_memory(n_sc);                    //free some memory in the  data_S
			//reset the actual number of scenarios (sub-problems)		          
			n_sc = n_sc - search_param->getMaster_sc(); // number of scenarios are total scenarios minus global scenarios (those who are not projected)			  
			*///initialize the search parameters
			arry_y_fixed 	=new int[data_S->getN_arcs()]();
			total_fix_Y	=0;		  
 			con_copy_binding  = new vector< int >();
			con_binding 	  = new vector< int >();
			con_copy_iter     = new vector< int >();
			con_iter 	  	  = new vector< int >();			  
			//initialize the global lower bound  
			GlobalLowerBound  	= new double[2];
			GlobalLowerBound[0]	= -1e75;	//lower bound
			GlobalLowerBound[1]	= -1e75;	//lower bound at the end of first phase before the transition phase
			//the information gathered from the master problem to be stored in the memory   
			Sol_send_size = data_S->getN_arcs() + n_sc  + 2 + 1;
			Sol_send  = new double[Sol_send_size]; // y sol and obja val, theta sol, Mp solving time ... 1 is for obj of the global scenarios					 
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
			    Managers->Sol_Pool_Manager(false, search_param, data_S, search_param->getMaster_sc(), n_sc, n_workers, Sol_send,  GlobalLowerBound[0], GlobalLowerBound[1]); //add the y variable and the theta solution to the pool and create the associated memory for tracking  v_sc_time,v_w_time
			    Managers->Sol_Pool_controler(false, n_workers, n_sc , data_S,  search_param);	//sends solutions to the workers "asynchronously" based on the given strategy ... he wont send if there is nothing to send					
			}  		    
			//create the master model  
			MasterModel(data_S, Sol_Choice,search_param,Managers); //this creates the master problem
			//cover inequlity class object
			if(search_param->getAddCoverIneq())
				coverCuts = new class_CI(data_S, n_sc);
			//cardinality inequlity class object
			if(search_param->getAddCardinalityIneq())
				cardinalityCuts = new class_MCI(data_S, n_sc);
			//adding intial cardinality cuts by inspecting the ...
			if(search_param->getAddRawCardinalityCuts()){
				raw_cardinality = new class_PMCI(data_S, n_sc);
				addRawCardinalityCuts(data_S, search_param,2);
				delete raw_cardinality;
			}
			optRound = new OptRound(Managers,  data_S, search_param,  n_sc) ;
		}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		bool solve_inter_phase(Data_S *data_S, bool stopping_ceriteria, double& interUB, double& interTime, double phase1Obj, double phase1Time, PoolsandMangers *Managers, Search_Param *search_param)
		{
			(cplex).setParam(IloCplex::TiLim, min(0.05*search_param->getTimeLim(), search_param->getTimeLim() - phase1Time));
			(cplex).setParam(IloCplex::EpGap, 0.003);
			// cplex.setOut(env.getNullStream());
			//start time of inter phase
			double interStart = MPI::Wtime();
			double m_objval;
			//get the current pool size
			int poolSize = Managers->getY_Sol_Poolsize();
			//fix the variables
			int counter =0;
			IloExpr exprInf(env);//if the fixation yielded infeasible problem, then at least one of those fixed to zero has to be opened
			if(search_param->getRunHardHeur()){
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
			}
			if(search_param->getRunSoftHeur()){
				for(int a =0; a < data_S->getN_arcs(); a++){
					if(y_SOL[a] <= 1e-2){
					   (*y_heur_con).add(IloRange(env, 0, y[a], 0));				   
					   counter++;
					}
					else if(y_SOL[a] <= 0.05)
						exprInf += y[a];
					if(y_SOL[a] >= 1-1e-2){
					   (*y_heur_con).add(IloRange(env, 1, y[a], 1));
					   counter++;
					}
					else if(y_SOL[a] >= 0.95)
						exprInf += 1- y[a];
				}
				(*y_heur_con).add(IloRange(env, 0, exprInf, 1));
			}
			exprInf.end();
 			(model).add((*y_heur_con));
			cout << endl << "num fixed var: " << counter << " out of " << data_S->getN_arcs() << endl ;
			//solve the restricted master model	
			for(int itr =0; itr < search_param->getNumHeurIter(); itr++)
			    if(solve_fix_opt_phase()){
				
				
					//extracting all the solutions from the pool
					// vector <double> aux;
					cout << endl << "************  "  << (cplex).getSolnPoolNsolns() << endl;
					for(int solPoolID=0; solPoolID < max((cplex).getSolnPoolNsolns()-1.0, 1.0); solPoolID++){
							// aux.clear();
							// (cplex).getValues(setSolVal, setSolVar, solPoolID);
							// IloExpr expr(env);
							// for(int a=0; a< data_S->getN_arcs(); a++){
								// aux.push_back(setSolVal[a]);							
								// if( setSolVal[a] <= 0.01 )
								  // expr 	+= y[a];							
								// if( setSolVal[a] >= 0.99 )
								  // expr 	+= 1 - y[a];
							// }
							// (*heurYPoolCon).add(IloRange(env, 1, expr, IloInfinity));
							// expr.end();
							// aux.push_back(cplex.getObjValue(solPoolID));
							// heurSolPool.push_back(aux);
					
							
							//
							   m_objval = cplex.getObjValue(solPoolID);
							   // if(cplex.getObjValue() == m_objval)
								// cout << endl << solPoolID << endl;
							  (cplex).getValues(setSolVal, setSolVar, solPoolID);							  
							   IloExpr expr(env);
							  cout <<endl;
							  for(int a=0; a< data_S->getN_arcs(); a++){
								if(setSolVal[a] >1e-3)
								  cout << a << "; ";						
								if( setSolVal[a] <= 0.01 )
									  expr 	+= y[a];							
								if( setSolVal[a] >= 0.99 )
									  expr 	+= 1 - y[a];														
							  }
							cout <<endl;
						  (*heurYPoolCon).add(IloRange(env, 1, expr, IloInfinity));//y[a].setUB(0);
		// 				  model.add((*y_heur_con));
						  expr.end(); 
						  // model.add(*heurYPoolCon);
						  
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
								  
						  Managers->Sol_Pool_Manager(true, search_param, data_S, search_param->getMaster_sc(), n_sc, n_workers, Sol_send,  GlobalLowerBound[0], GlobalLowerBound[1]); //add the y variable and the theta solution to the pool and create the associated memory for tracking  v_sc_time,v_w_time
						  GlobalLowerBound[1]=GlobalLowerBound[0];
						  Managers->Sol_Pool_controler(true, n_workers, n_sc , data_S,  search_param);	//sends solutions to the workers "asynchronously" based on the given strategy ... he wont send if there is nothing to send					
						  if(search_param->getCreate_art_SPs())
							Managers->art_cut_pool_manager(n_sc, n_workers, GlobalLowerBound[0], data_S, search_param);			 				            				      							  
						  Managers->cut_pool_manager(false, stopping_ceriteria, n_sc, n_workers, GlobalLowerBound[0], data_S,  search_param); //v_sc_time,v_w_time, return number of cuts that we can receive			 				            				      							  
						  interPhaseOptimal=false;
						  if((100*(Managers->getGlobalUpperBound(0)-phase1Obj)/Managers->getGlobalUpperBound(0)) <= 1){
							interPhaseOptimal=true;
							break;
						  }
					}
					add_regular_cuts(data_S, Managers, search_param);
					model.add(*heurYPoolCon);					
									  
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
			interTime	= MPI::Wtime() - interStart;
			cout << "time of the intermidiary phase: " << interTime << endl;	
			
			return interPhaseOptimal;
		}
///////////////////////////////////////////////////////
		bool solve_fix_opt_phase()
		{
			    bool fesProb=true;
			    if( !(cplex).solve() ){
			      fesProb=false;
			      cout << "****** Danggerrrrrr the injected solution will be rejected: " << (cplex).getStatus() << endl;
 			      // abort();
			    }
			   
			    return fesProb;
		  }
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		void StopWorkers(Search_Param *search_param)
		{
		   	////we tell working processors to finish 
			continue_running = false;
			for(int w_id = 0; w_id < n_workers; w_id++)
			  MPI::Comm.Isend(&continue_running, 1, MPI::BOOL, w_id + 1, search_param->getTagconti_run()); 				      									      
		}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		void MasterModel(Data_S *data_S, MasterSolChoice *Sol_Choice, Search_Param *search_param, PoolsandMangers *Managers) //constructor
		{		  
			model 				= IloModel(env);
			cplex 				= IloCplex(model);
			con   				= new IloRangeArray(env);
			con_optcut  		= new IloRangeArray(env);
			con_copy   			= new IloRangeArray(env);
			artificialCon			= new IloRangeArray(env);
			y_heur_con			= new IloRangeArray(env);
			heurYPoolCon		= new IloRangeArray(env);
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
					x_Sol[a][k] = IloNumArray(env, search_param->getMaster_sc()+search_param->getNumArtificialScenarios());			      
			}
			
			reducedCosts = IloNumArray(env, data_S->getN_arcs());
			//***the master only includes y and theta variables****
			y = IloNumVarArray(env, data_S->getN_arcs(), 0, 1);	//ILOINT defining an variable for each arc (on dimensional variable)						
			theta = IloNumVarArray(env, n_sc + search_param->getMaster_sc(), 0, IloInfinity);	//defining variables to approximate the recourse costs
											//I generate a theta for the global scenarios too, but keep their coeficient to zero (just to ease the programming pain :P)
			(model).add(y); // we dont need to add variables explicitly to the model as they will implicitly be added to the model through the constraints
			(model).add(theta);
			//adding x variable to the master for the global scenarios 
			if(search_param->getMaster_sc() + search_param->getNumArtificialScenarios()> 0){
			    x = IloNumVarArray3(env, data_S->getN_arcs());//second-stage (x) variable (three dimensional cuz each arc (i-j) is represented with a single name "a"
			    for(int a = 0; a < data_S->getN_arcs(); a++){
				    x[a] = IloNumVarArray2(env, data_S->getN_od());
				    for(int k = 0; k < data_S->getN_od(); k++)				
					    x[a][k] = IloNumVarArray(env, search_param->getMaster_sc()+search_param->getNumArtificialScenarios(), 0, IloInfinity); 									
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
			if(search_param->getNumArtificialScenarios()>0)
			  for(IloInt a=0; a<data_S->getN_arcs(); a++)
			      for(IloInt k=0; k< data_S->getN_od(); k++)				  
				    setSolVar.add(x[a][k][search_param->getMaster_sc()]);					
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
			for(int s = 0; s < search_param->getMaster_sc(); s++){		// flow conservation constraints for global scenarios			
				ss = Managers->getglobal_scenario_id(s);//global_scenario_id[s];
				for(int k = 0; k < data_S->getN_od(); k++)				  
					for(int i = 0; i < data_S->getN_nodes(); i++){ 
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
			//dummies are always open
			// for(int a = data_S->getN_arcs()- data_S->getN_od(); a < data_S->getN_arcs(); a++)			
				// model.add(IloRange(env,1,y[a],1));
			//valid inequalities	
			RBF_Ineq(data_S, Sol_Choice, search_param); //for the recourse cost bounding functions 
			if(search_param->getLooping_arcs())
				looping_arcs(data_S);
			if(search_param->getLestCapIneq())
				min_total_capacity(data_S, search_param);
			
			//constraints related to artifical scenario			
			(model).add((*con));
			//
			if(search_param->getNumArtificialScenarios()>0){
				for(int k = 0; k < data_S->getN_od(); k++)				  
					for(int i = 0; i < data_S->getN_nodes(); i++){ 
						  IloExpr expr(env);//expr.clear();
						  IloExpr expr1(env);
						  for(int a =0; a < data_S->getN_arcs()/*- data_S->getN_od()*/; a++){
							  if(data_S->getArcs(a,0) == i) 
								  expr += x[a][k][search_param->getMaster_sc()];								  
							  if(data_S->getArcs(a,1) == i) 
								  expr1 -=  x[a][k][search_param->getMaster_sc()];
						  }
						  if (i == data_S->getOd(k,0))//if origin
						    (*artificialCon).add(IloRange(env, data_S->getArtificialScenario(k), expr+expr1, data_S->getArtificialScenario(k)));//for the origin only outgoing part is neccessary						
						 else
						   if (i == data_S->getOd(k,1))//if destination
							(*artificialCon).add(IloRange(env, -data_S->getArtificialScenario(k), expr+expr1, -data_S->getArtificialScenario(k)));						
						   else
							(*artificialCon).add(IloRange(env , 0, expr + expr1, 0));						     
						  expr.end();
						  expr1.end();
					}
				for(int a = 0; a < data_S->getN_arcs()/*- data_S->getN_od()*/; a++){
					IloExpr expr(env);					
					for(int k = 0; k < data_S->getN_od(); k++)
						expr += x[a][k][search_param->getMaster_sc()];					
					expr += -data_S->getU(a)*y[a];
					(*artificialCon).add(IloRange(env,-IloInfinity,expr,0));
					expr.end();
				}
				IloExpr expr(env);
				expr.clear();
				for(int a = 0; a < data_S->getN_arcs() /*- data_S->getN_od()*/; a++)			  
					for(int k = 0; k < data_S->getN_od(); k++)			//the second stage objective for the global scenarios 				   			     
						  expr +=  data_S->getC(a) * x[a][k][search_param->getMaster_sc()];//(1/master_sc) * c[a] * x[a][k][s];//
				for(int s=0; s< n_sc; s++)
				  expr -= data_S->getAlphaWeights(s)*theta[s];
				(*artificialCon).add(IloRange(env,0,expr,0));
				expr.end();
				(model).add(*artificialCon);
			}
			//
			(cplex).setParam(IloCplex::RandomSeed, 2015);//fix cplex random seed
			(cplex).setParam(IloCplex::TiLim, 36000);
			(cplex).setParam(IloCplex::ClockType, 1);//to set it equal to cpu time since it is seq
			(cplex).setParam(IloCplex::SimDisplay, 0);			 					
			(cplex).setParam(IloCplex::EpGap, 0.01);
			(cplex).setParam(IloCplex::AdvInd, 0);
			// (cplex).setParam(IloCplex::RootAlg, IloCplex::Primal);
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
		bool addViloatedCoverIneq(Data_S *data_S, int i, Search_Param *search_param)
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
		bool addViolatedCardinalityCuts(Data_S *data_S, int i,Search_Param	*search_param)
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
///////////////////////////////////////////////////////////////////////////////////////		
		void addRawCardinalityCuts(Data_S *data_S, Search_Param	*search_param, int cutSetCardinality)
		{
			//generate the cuts
			raw_cardinality->MainLoop(data_S,  n_sc + search_param->getMaster_sc(), cutSetCardinality);
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
		void looping_arcs(Data_S *data_S)
		{
			for(int a =0; a < data_S->getN_arcs() - data_S->getN_od() ; a++)
				if(data_S->getArcs(a,0) == data_S->getArcs(a,1))
					(*con).add(IloRange(env, 0, y[a] ,0 ));
		}//looping arcs								
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		void RBF_Ineq(Data_S *data_S, MasterSolChoice *Sol_Choice, Search_Param	*search_param)
		{				
				 //LBF inequality
				  if(search_param->getLB_Lifting())//if we want to add LB lifting inequalities
				  {
					class_LB_Lifting *LBF;
					
					for(int s=0 ;s <n_sc; s++){
						LBF = new class_LB_Lifting(n_sc, Sol_Choice, data_S, search_param, s);
						IloExpr expr(env);						
						expr += LBF->getRotuingCost();
						for(int k=0; k< data_S->getN_od(); k++)
						  expr += (data_S->getD(k,s) - LBF-> getMinD(k)) * LBF-> getDualVal(k);
						for(int a=0; a < data_S->getN_arcs() - data_S->getN_od(); a++)
						  expr += data_S->getF(a)* (LBF->getY(a) -y[a]);
						  
						(*con).add(IloRange(env, -IloInfinity, expr - theta[s], 0 ));
						expr.end();
						delete LBF;
					  }					  
					  //lifting y variables based on LBF and recourse upper boudning 
					 /* if(  search_param->getY_Lifting())
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
				}*/	
			}
			}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		void connectivity_constraints(Data_S *data_S, bool LP_Phase, MasterSolChoice *Sol_Choice, Search_Param	*search_param)
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
				 cout << endl << "added connectivity cuts " << (*connectivity_con).getSize() << endl;
				 ODconnectivityCuts(data_S, search_param);
			 }
			  //for the second phase, if it is R data set update the constraints
			  if(/* data_S->getInsClass()  == "R" && */ !LP_Phase){
					// (model).remove(*connectivity_con);
					(*connectivity_con).endElements();//clear();
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
				  if((*connectivity_con).getSize() > 0)
					(model).add(*connectivity_con);
				 cout << endl << "added connectivity cuts " << (*connectivity_con).getSize() << endl;
			  }
			  
			  delete NC;
		}	
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		void ODconnectivityCuts(Data_S *data_S, Search_Param *search_param)
		{			
			//adding initial cardinality cuts by inspecting the ...
			if(!search_param->getAddRawCardinalityCuts()){
				raw_cardinality = new class_PMCI(data_S, n_sc);
				addRawCardinalityCuts(data_S, search_param,1);
				delete raw_cardinality;
			}
		}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		void min_total_capacity(Data_S *data_S, Search_Param	*search_param)
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
		bool strong_inequality(Data_S *data_S, PoolsandMangers *Managers, Search_Param	*search_param)
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
			//for the created artificial scenario 
			if( search_param->getNumArtificialScenarios() > 0){
				for(int k = 0; k < data_S->getN_od(); k++)		
				    for(int a = 0; a < data_S->getN_arcs()- data_S->getN_od(); a++)
						if(x_Sol[a][k][search_param->getMaster_sc()] > 1.0e-2){
							if (x_Sol[a][k][search_param->getMaster_sc()] - min(data_S->getU(a),data_S->getArtificialScenario(k)) * y_SOL[a] > 1e-2){
								IloExpr expr(env);
								expr += x[a][k][search_param->getMaster_sc()];						  
								expr -= min(data_S->getU(a),data_S->getArtificialScenario(k)) * y[a];								
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
		double solve_phase_I(Data_S *data_S, bool M_D_Phase_I, PoolsandMangers *Managers, Search_Param	*search_param)//, vector<int> *v_ac, vector<int> *v_ac_copy, vector<int> *v_inac, vector<int> *v_inac_copy, vector<int> *v_acinac_aux ) 
		{                             			
			double m_objval, tStart, tEnd;//, *slackvals;
			double Recourse_Obj = 0;			//the second stage obj associated to the global scenarios
			//---------------------solve ------------------------------
			tStart= (cplex).getTime();	//this on linux gives cpu time not the wall clock
						
			if( !(cplex).solve() ) 
			{			  			
			  env.error()  << " ** Failed to optimize Master **> " <<  (cplex).getCplexStatus() << endl;			 
			  abort();
			  return 1;
			}
			//to clean up the master
			if(search_param->getCleanMaster())
			  clean_master(data_S, M_D_Phase_I);
  			cout << endl << "Time spend on master " << (cplex).getTime() - tStart<<endl;	
			//retrive the required values			
			/* (cplex).getReducedCosts(reducedCosts, y);
			env.out() << "Reduced Costs = " << reducedCosts << endl; */
			
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
// 					if(search_param->getNumArtificialScenarios()>0)
// 					  t += search_param->getNumArtificialScenarios();
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
		void clean_master(Data_S *data_S, bool M_D_Phase_I)		   
		{			    
			    //regular cuts			    
			     if((*con_optcut).getSize() >0 && !M_D_Phase_I){			       
					IloNumArray  Reg_slack(env);
					(cplex).getSlacks(Reg_slack, *con_optcut);
					for(int i = 0; i <(*con_optcut).getSize(); i++) {//for every cut in the MP					   				     
						if ((*con_binding).size() > i){										      
							  if(Reg_slack[i] <= 1e2)//active					  
							  (*con_binding)[i]++;//increase binding counter
						}
						else{ 					
							  if(Reg_slack[i] <= 1e2)					      
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
			      cout << (*con_optcut).getSize() - (temp_con).getSize() << endl;			      					
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
		 }
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////	
		int add_regular_cuts(Data_S *data_S, PoolsandMangers *Managers, Search_Param	*search_param) //"NOTE: we only have optimality cuts
		{		  			  
			  double  g, g1;
			  int sc_id, sol_id, vio_check, Num_added_cuts=0, Num_added_cuts1=0, solId=1e-2;//the last one is to check we are not adding the same CBC cut several times
			  // double viol=0.0;

			
			for(int j = 0; j < search_param->getNumAggrCluster(); j++){
				IloExpr expr(env), expr1(env);
				for(int k = 0; k < Managers->getClusterSize(j); k++)
					for(int i = 0; i < Managers->getv_cpSize(); i++)
						if(Managers->getv_cp(i, 0)== Managers->getY_Sol_Poolsize()-1 && Managers->getv_cp(i, 1) == Managers->getClusterIDs(j, k)  ){
							// sol_id = Managers->getv_cp(i, 0);//(*v_cp)[i][0];	//solution ID
							// cout <<  j << " , " << Managers->getv_cp(i, 1) << " - " << Managers->getClusterIDs(j, k) << endl;
							// cout << Managers->getv_cp(i, 0) << " - " << Managers->getv_cp(i, 1)<< endl;
							vio_check=0;
							sol_id = Managers->getY_Sol_Poolsize()-1;
							sc_id  = Managers->getv_cp(i, 1);//(*v_cp)[i][1];	//this indicates which scenario this cut belongs to
							for(int a = 0; a < data_S->getN_arcs(); a++)
								vio_check += Managers->getv_cp(i, a+3)* Managers->getY_Sol_Pool(sol_id, a) ;
							vio_check += Managers->getv_cp(i, 2);							
							vio_check -= Managers->getV_theta_sol(sol_id, sc_id);//theta[sc_id];
							if(vio_check >1){
								Num_added_cuts++;
								// sc_id  = Managers->getv_cp(i, 1);//(*v_cp)[i][1];	//this indicates which scenario this cut belongs to
								for(int a = 0; a < data_S->getN_arcs(); a++)
										expr += Managers->getv_cp(i, a+3)*y[a];
								expr += Managers->getv_cp(i, 2);							
								expr -= theta[sc_id];
							}
						}
						else if(Managers->getv_cp(i, 0) != Managers->getY_Sol_Poolsize()-1 && Managers->getv_cp(i, 1) == Managers->getClusterIDs(j, k) ){
							vio_check=0;
							sol_id = Managers->getY_Sol_Poolsize()-1;
							sc_id  = Managers->getv_cp(i, 1);//(*v_cp)[i][1];	//this indicates which scenario this cut belongs to
							for(int a = 0; a < data_S->getN_arcs(); a++)
								vio_check += Managers->getv_cp(i, a+3)* Managers->getY_Sol_Pool(sol_id, a) ;
							vio_check += Managers->getv_cp(i, 2);							
							vio_check -= Managers->getV_theta_sol(sol_id, sc_id);//theta[sc_id];							
							if(vio_check>1){
								Num_added_cuts1++;
								sc_id  = Managers->getv_cp(i, 1);//(*v_cp)[i][1];	//this indicates which scenario this cut belongs to
								for(int a = 0; a < data_S->getN_arcs(); a++)
										expr1 += Managers->getv_cp(i, a+3)*y[a];
								expr1 += Managers->getv_cp(i, 2);							
								expr1 -= theta[sc_id];
							}
						}
				if(Num_added_cuts>0)
					(*con_optcut).add(IloRange( env, -IloInfinity, expr, 0));//, con_name));					      
				if(Num_added_cuts1>0)
					(*con_optcut).add(IloRange( env, -IloInfinity, expr1, 0));
				expr.end(); 
				expr1.end(); 
			}
			  
			if(Num_added_cuts >0)
			    (model).add((*con_optcut));
			// cout << endl << "nym of added cuts: " <<Managers->getv_cpSize()- Num_added_cuts << endl;
			//clean up the cut pool		 			 			
			Managers->Cleanv_cp();
			// if(search_param->getPropagatCut())
				// add_copied_cuts2(data_S, Managers, search_param);
			if(search_param->getCreate_art_SPs())
				add_art_cuts(data_S, Managers, search_param);
			return Num_added_cuts;
		}
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////	
		int add_art_cuts(Data_S *data_S, PoolsandMangers *Managers, Search_Param	*search_param) //"NOTE: we only have optimality cuts
		{		  			  
			  double  g, g1;
			  int sc_id, sol_id, vio_check, Num_added_cuts=0, solId=1e-2;//the last one is to check we are not adding the same CBC cut several times
			  // double viol=0.0;

			for(int i = 0; i < Managers->getv_cp_artSize(); i++)
			{
				sol_id = Managers->getv_cp_art(i, 0);//(*v_cp)[i][0];	//solution ID
				sc_id  = Managers->getv_cp_art(i, 1);//(*v_cp)[i][1];	//this indicates which scenario this cut belongs to
				g      = Managers->getv_cp_art(i, 2);//(*v_cp)[i][2];
		  
			   
					Num_added_cuts++;				    
				    IloExpr expr(env);				    
				    for(int a = 0; a < data_S->getN_arcs(); a++){
						expr += Managers->getv_cp_art(i, a+3)/*(*v_cp)[i][a+3]*/*y[a];
						// viol += Managers->getv_cp(i, a+3) *  Managers->getY_Sol_Pool(sol_id, a) ;
				    }
				    expr += g;
				    // viol += g;
					int w_id = sc_id-n_sc;
					double prob_sum=0.0;
					for(int s=0; s< Managers->get_aux_worker_sc_size(w_id); s++)
						prob_sum += data_S->getP(Managers->get_aux_worker_sc_ids(w_id, s));
					for(int s=0; s< Managers->get_aux_worker_sc_size(w_id); s++)
							expr -= (data_S->getP(Managers->get_aux_worker_sc_ids(w_id, s))/prob_sum) * theta[Managers->get_aux_worker_sc_ids(w_id, s)];
				    // viol -= Managers->getV_theta_sol(sol_id, sc_id);
				    (*con_optcut).add(IloRange( env, -IloInfinity, expr, 0));//, con_name));					      
				    expr.end();
			}
			  
			if(Num_added_cuts >0)
			    (model).add((*con_optcut));
			//clean up the cut pool		 			 			
			Managers->Cleanv_cp_art();
			// Managers->Cleanv_cp_copy();
			return Num_added_cuts;
		}
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////			  
		bool IntegerSol(Data_S *data_S, int solId, PoolsandMangers *Managers)	//to check if a specific solution is integer
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
		int add_copied_cuts2(Data_S *data_S, PoolsandMangers *Managers, Search_Param	*search_param) 
		{
			  double  g, g1;
			  int sc_id, sol_id, vio_check, Num_added_cuts=0, solId=1e-2;//the last one is to check we are not adding the same CBC cut several times
			  double viol=0.0;
			  cout << endl << "cop pool size: " << Managers->getv_cp_copySize() << endl;
			for(int i = 0; i < Managers->getv_cp_copySize(); i++)
			/* if(Managers->getv_cp_copy(i, 1) >= 64) */{
				viol=0.0;
				sol_id = Managers->getv_cp_copy(i, 0);//(*v_cp)[i][0];	//solution ID
				sc_id  = Managers->getv_cp_copy(i, 1);//(*v_cp)[i][1];	//this indicates which scenario this cut belongs to
				g      = Managers->getv_cp_copy(i, 2);//(*v_cp)[i][2];
				IloExpr expr(env);				    
				for(int a = 0; a < data_S->getN_arcs(); a++){
					expr += Managers->getv_cp_copy(i, a+3)*y[a];
					viol += Managers->getv_cp_copy(i, a+3)*setSolVal[a] ;
				}
				expr += g;
				viol += g;
				expr -= theta[sc_id];
				viol -= setSolVal[data_S->getN_arcs() + sc_id];//Managers->getV_theta_sol(sol_id, sc_id);
				if(/* true &&  */viol >= 10){
					cutIDsToRemove.push_back(i);
					(*con_optcut).add(IloRange(env, -IloInfinity, expr, 0));//, con_name));
					Num_added_cuts++;						
				}
				expr.end();			
			}			  
			if(Num_added_cuts >0){
			    (model).add((*con_optcut));
				numCopyCutUsed += Num_added_cuts;
				cout << endl << "^^^^^^^^^^^^^Num of CopyCut used: " <<  numCopyCutUsed  << endl;			
			}
			//clean up the cut pool	
			// Managers->Cleanv_cp_copy();
			// for(int i=0; i< cutIDsToRemove.size(); i++)
				// Managers->Cleanv_cp_copyID(i);
			return Num_added_cuts;
		}
///////////////////////////////////////////////////////////////////////////////////////////				
		void add_copied_cuts(Data_S *data_S, PoolsandMangers *Managers, Search_Param	*search_param) 
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
 			  cout << endl << "^^^^^^^^^^^^^The size of copied Bender cuts are: " <<  (*con_copy).getSize()  << endl;
			}
						 			 			
			Managers->Cleanv_cp_copy();
			
		}
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		void solve_phase_II(Data_S *data_S, double wtime, MasterSolChoice *Sol_Choice, PoolsandMangers *Managers, 	Search_Param	*search_param)
		{		    
				(cplex).setParam(IloCplex::EpGap, 0.01);
				(cplex).setParam(IloCplex::Param::MIP::Strategy::LBHeur, 1);
				(cplex).setParam(IloCplex::TiLim, search_param->getTimeLim() -  wtime);
				// Set up the cut callback to be used for separating Benders' cuts
//  			     (cplex).setParam(IloCplex::Param::Preprocessing::Presolve, IloFalse); 
//  			    (cplex).setParam(IloCplex::Param::Threads, 1); 
//  			     (cplex).setParam(IloCplex::Param::MIP::Strategy::Search, IloCplex::Traditional);
		
   			    (cplex).use(BendersLazyCallback(env,x,y,theta, Managers, Sol_Choice, data_S, search_param));//
 			    (cplex).use(MyHeuristic(env, x, y, theta, Managers, data_S, search_param));
 			    (cplex).use(ControlIncumSol(env, Managers));
			    (cplex).use(MyBranch(env, y, Managers, data_S));
			    
			    
				if(search_param->getComm_strategy() != 0)
					(cplex).use(BendersUserCallback(env,x,y,theta, Managers, Sol_Choice, data_S, search_param));
			    				
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
			   // env.out() << "Reduced Costs = " << reducedCosts << endl;	
			}		   
		} 
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		void addReducedCostBasedIneq(Data_S *data_S, PoolsandMangers *Managers)
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
		double getGlobalLowerBound(int i){return GlobalLowerBound[i];}
		void setGlobalLowerBound(int i, double GlobalLowerBound0){ GlobalLowerBound[i]=GlobalLowerBound0;}
		double* getSol_send()  {return Sol_send;}
		int getSol_send_size() {return Sol_send_size;}
		int getreducedCostsSize(){return reducedCosts.getSize();}
		double getreducedCosts(int a){return reducedCosts[a];}
		double gety_SOL(int a){return y_SOL[a];}
		~class_MP()
		{
		}


};//WND class



#endif
