#ifndef CLASS_MAIN_H
#define CLASS_MAIN_H


#include <ilcplex/ilocplex.h>
#include "./class_search_param.h"
#include "./Sub_Classes/class_data_S.h"

#include "./Sub_Classes/class_poolsandmasgers.h"
#include "./Sub_Classes/Sol_Choice.h"
#include "./Sub_Classes/class_MP.h"
#include "./Sub_Classes/class_DEF.h"


ILOSTLBEGIN
#define Comm COMM_WORLD



class class_main
{ 
  private:    
		//object of the required classes
		class_MP *MP;
		class_DEF *DEF;//deterministic equivalent
		PoolsandMangers *Managers;
		MasterSolChoice *Sol_Choice;		
		Search_Param *search_param;
		Data_S *data_S; 		        //create the objects of the data class
		//control and output variables 
		int n_sc, n_workers;
		double *auxSol_send;
		int  Equal_LB, iteration;
		double Best_Lower_Bound, wtime, tEndLS, phase1Time, phase1Obj, LB_best_time, interTime, interLB, interUB;
		bool transition_phase, M_D_Phase_I;
		bool interPhaseOptimal, stopping_ceriteria, addProgCuts;		  
  public:      
		class_main(int n_workers0, int n_sc0, string instName, string scenarioName, string insClass)		
		{
			//time holder
			wtime = MPI::Wtime();
			n_sc 		=n_sc0;
			n_workers 	=n_workers0;   
			//create the object of the classes   
			search_param = new Search_Param(n_sc,n_workers);
			DEF = new class_DEF();
			
			data_S = new Data_S(insClass); 		      	  //create the objects of the data class 	
			data_S->info(instName);		      			//here we fill the network information		  			  
			data_S->info_scenario(search_param, n_sc, scenarioName);   //this read the scenario file and fill the p array (probabilities)			    	          	 			 
			// data_S->SceCreRetention(search_param, n_sc);
			data_S->read_assign(n_sc);	              	//fill the arrays and vectors	  				
			if( search_param->getData_Clean_UP_Moode() ) //clean up data	    		      			   		 
				data_S->clean_data(n_sc);	       			    
			data_S->clean_memory(n_sc);                    //free some memory in the  data_S
			cout << endl << "Time spend on creation-retention problem is: " << MPI::Wtime() - wtime << endl ;
			//reset the actual number of scenarios (sub-problems)		          
			n_sc = n_sc - search_param->getMaster_sc(); // number of scenarios are total scenarios minus global scenarios (those who are not projected)
			//print the information of the instances and the programm
			PrintInstInfo();			  
			//initialize the search paramters
			iteration	=0;	//counter of number of iterations
			M_D_Phase_I =false;	//indicating if the first phase of the algorithm is done		  			  			  			 				  					  
			//cluster and assign scenarios to the workers  
			Managers  = new PoolsandMangers(n_workers, n_sc, data_S, search_param);				
			Managers->Cluster_Sce(search_param, data_S, n_workers, n_sc + search_param->getMaster_sc(), search_param->getMaster_sc());//here we cluster the scenarios
			Managers->Assign_Sce_Proc(n_workers, search_param);//here we assign each processor (including master) with its scenarios			  		  
			//let slaves start from core point generating cut	
			Sol_Choice    = new MasterSolChoice(data_S, search_param, n_sc); //this creates the master problem			  			
			MP  = new class_MP(data_S, Sol_Choice, search_param, Managers, n_workers, n_sc/*+search_param->getMaster_sc()*/, instName, scenarioName, insClass);
			auxSol_send  = new double[MP->getSol_send_size()];
			Equal_LB = ceil(n_sc/search_param->getN_cut_wait())+1;		//to check if the solution from current call of the master problem has been stored in the memory (true) or not (false);; how many equal lower bounds we will see in first phase before ending it
			Best_Lower_Bound =0; 							//to have the best lower bound (and its time)
			stopping_ceriteria =false, transition_phase =false;	
			addProgCuts=true;
		}
////////////////////////////////the master processor main loop////////////////////////////////////////////////////////////////////////////////////////////////////////////////////		
		void main_loop()
		{		
			//add network connectivity cuts
			if(search_param->getNetwork_Connectivity())
			  MP->connectivity_constraints(data_S, !M_D_Phase_I, Sol_Choice, search_param);
			//initialize the paramters which will be used 
			bool cutLoop=true;//to keep solving the master problem until new valid inequalities are added
			bool boolCI1=false, boolCI2=false, boolMCI1=false, boolMCI2=false;//is new conver of size 1,2 is found; is cardinality of size 1,2 is found
			//the main loop of the algorithm
			wtime = MPI::Wtime();
			while( !stopping_ceriteria )				//while stopping criteria are false 
			{
				iteration++;		
				//check for the first phase stopping statisfaction 
			  	if(!M_D_Phase_I){//M-D First Phase	
					auxSol_send = MP->getSol_send();
				    Managers->first_phase_ending(search_param, phase1Time, phase1Obj, wtime, MP->getGlobalLowerBound(0), Equal_LB, transition_phase, M_D_Phase_I, stopping_ceriteria, data_S->getN_arcs(), search_param->getTolerance(), auxSol_send);									
				    if ( stopping_ceriteria )//termintation	, because the first phase might be optimal
				      break;
				    if(M_D_Phase_I){//we dont want to receive cuts for the solution discovered in the first phase during the second phase				    
				      Sol_Choice->clean_phase_I_sol( Managers, n_workers, data_S->getN_arcs(), n_sc, search_param->getTag_Useless_Sol());		    
				    }
				}
				//the cut loop of the first phase if it is not done
				if(!M_D_Phase_I)//I dont want to still solve the master for one more iteration even if the first phase is done
				{
				    //if we give the master a warm start by starting from the core point
				    if(iteration==1 && search_param->getWroker_start_from_core_sol()){
						Managers->cut_pool_manager( true, stopping_ceriteria, n_sc, n_workers, MP->getGlobalLowerBound(0), data_S,  search_param); //v_sc_time,v_w_time, return number of cuts that we can receive			 				            				      	
						MP->add_regular_cuts(data_S, Managers, search_param);
  // 				        add_copied_cuts();
				    }	



					











					
				    //we add valid ineq as much as they exist before adding benders cut
					cutLoop=true;
					while(cutLoop){
						//solve the master
						double GlobalLowerBound0 = MP->solve_phase_I(data_S, M_D_Phase_I, Managers, search_param);//, v_ac, v_ac_copy, v_inac, v_inac_copy, v_acinac_aux);										  
						MP->setGlobalLowerBound(0, GlobalLowerBound0);
						//check for strong ineq
						if(true)
						  MP->strong_inequality(data_S, Managers, search_param);	
						//cover inequalities 
						if(search_param->getAddCoverIneq())
							boolCI1 = MP->addViloatedCoverIneq(data_S, 1, search_param);				
						//cover inequalities of size 2
						if(search_param->getAddCoverIneq())
							boolCI2 = MP->addViloatedCoverIneq(data_S, 2, search_param);
						//if no cover exist then check for the cardinalities 
						if(!boolCI2 && !boolCI1){
						  //cardinality inequalities
						  if(search_param->getAddCardinalityIneq())
							  boolMCI1 = MP->addViolatedCardinalityCuts(data_S,1,search_param); 
						  //cardinality inequalities of size 2
						  if(search_param->getAddCardinalityIneq())
							  boolMCI2 = MP->addViolatedCardinalityCuts(data_S,2,search_param);
						}
						//if no wnew cut is derived then break the loop 
						if(!boolCI1 && !boolMCI1 && !boolCI2 &&  !boolMCI2)
						  cutLoop=false;						
				    }
					
					
					
					
				    //going through the copy cuts to find and add violated ones. 
					if(search_param->getPropagatCut() && addProgCuts && iteration > 1){
						// while (true){
							MP->add_copied_cuts2(data_S, Managers, search_param);
							double GlobalLowerBound0 = MP->solve_phase_I(data_S, M_D_Phase_I, Managers, search_param);//, v_ac, v_ac_copy, v_inac, v_inac_copy, v_acinac_aux);										  
							cout << endl << "lower bound before / after adding the propagated cuts: " << MP->getGlobalLowerBound(0) << " / " << GlobalLowerBound0 << endl ;
							if(GlobalLowerBound0 < MP->getGlobalLowerBound(0) + 1 )
								addProgCuts=false;
							MP->setGlobalLowerBound(0, GlobalLowerBound0);	
							// break;  
						// }
						// //we clean memory heritage from phase 1						
					}
					Managers->Cleanv_cp_copy();











					
					auxSol_send = MP->getSol_send();					
					Managers->Sol_Pool_Manager(false, search_param, data_S, search_param->getMaster_sc(), n_sc, n_workers, auxSol_send, MP->getGlobalLowerBound(0), MP->getGlobalLowerBound(1)); //add the y variable and the theta solution to the pool and create the associated memory for tracking  v_sc_time,v_w_time
					MP->setGlobalLowerBound(1, MP->getGlobalLowerBound(0));// GlobalLowerBound[1] = GlobalLowerBound[0];
				    					    
				    //I wrote this part to have the time when the best lower bound is found
				    if (MP->getGlobalLowerBound(0) > Best_Lower_Bound+0.5)
				    {
				      tEndLS = MPI::Wtime() - wtime;
				      Best_Lower_Bound = MP->getGlobalLowerBound(0); 
				      cout << endl << "the best lower bound is found at the time: " <<  tEndLS << " with value of: " << Best_Lower_Bound << "  @ iteration: " << iteration-1 << endl;
				      LB_best_time = tEndLS;
				    }					    
				    
				    Managers->Sol_Pool_controler(M_D_Phase_I, n_workers, n_sc , data_S,  search_param);	//sends solutions to the workers "asynchronously" based on the given strategy ... he wont send if there is nothing to send					
				    if(search_param->getCreate_art_SPs())
						Managers->art_cut_pool_manager(n_sc, n_workers, GlobalLowerBound[0], data_S, search_param);			 				            				      							  
					Managers->cut_pool_manager(true, stopping_ceriteria, n_sc, n_workers, MP->getGlobalLowerBound(0), data_S,  search_param); //v_sc_time,v_w_time, return number of cuts that we can receive			 				            				      					    
// 				    Sol_Choice->regular_cuts(data_S->getN_arcs(), n_sc, Managers); //add the recived cuts to the master
// 				    Sol_Choice->copied_cuts(Sol_send, n_sc, Managers, data_S,search_param); //add the recived cuts to the master				    
// 				    Sol_Choice->solve_choice_mas(n_workers, n_sc, Managers, data_S, search_param);		    
// 				    
				    MP->add_regular_cuts(data_S, Managers, search_param); //add the received cuts to the master
					
					
					
					
					
					
					
					
					
					
					
					
					
					
					
					
					
					
	// 				    add_copied_cuts(); //add the received cuts to the master
						/* if(search_param->getRunHeuristic())
							solve_inter_phase(); */					
				}
				else//M-D Second Phase
				{
					if(true){// break;
					Managers->Cleanv_cp_copy();//we clean memory heritage from phase 1
					if(search_param->getNetwork_Connectivity())
						MP->connectivity_constraints(data_S, !M_D_Phase_I, Sol_Choice, search_param); 	
					//here we get the reduced costs we need, 
					if(search_param->getReducedCostFixing()){
						MP->ReducedCosts(); //to get the reduced costs at the root node
					}
					if(search_param->getRunWholeAlgAsHeuristic()){
						double heursticObj = DEF->main_loop(n_sc, search_param, data_S, MP);
						Managers->setGlobalUpperBound(0,heursticObj);
					}
					else{
						//clean master problem
						if(search_param->getCleanMaster())
							MP->clean_master(data_S, M_D_Phase_I);//before branching we remove those cuts that have larger slack values
						 //impose the integrality requirement 
						 MP->ConvertLPtoMIP();	//here we turn into integer those who are not	
		  // 				 heuristic();	//we attempt to inject feasible solution for the second phase based on the solutions in the first phase					
						if(search_param->getRunHeuristic()){						
							interPhaseOptimal = MP->solve_inter_phase(data_S, stopping_ceriteria, interUB, interTime, phase1Obj, phase1Time, Managers, search_param);
						}  
						//look for the reduced cost based cuts
						if(search_param->getReducedCostFixing()){
							MP->addReducedCostBasedIneq(data_S, Managers); 
						}
						if(!interPhaseOptimal && (search_param->getTimeLim() - (MPI::Wtime() - wtime) > 1) && (((Managers->getGlobalUpperBound(0)-phase1Obj)/Managers->getGlobalUpperBound(0)) > 0.003) ){//start second phase if any time is left						
							MP->solve_phase_II(data_S, MPI::Wtime() - wtime, Sol_Choice, Managers, search_param); 
						}
	// 				    cplex.exportModel("model.lp");			
						
					}
					}
					break; 					
				}								
			} //end while
			
			PrintResults();  
			MP->StopWorkers(search_param);		
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
			
			cout << endl <<  "The global lower bound: " << MP->getGlobalLowerBound(0) << endl;			
			cout << endl << "The global upper bound: " 	<< Managers->getGlobalUpperBound(0) << endl; 
			cout << endl << "Optimality gap in % is: " 	<< 100*((Managers->getGlobalUpperBound(0) - MP->getGlobalLowerBound(0))/Managers->getGlobalUpperBound(0)) << endl;			
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
		~class_main()
		{
		}


};//WND class



#endif
