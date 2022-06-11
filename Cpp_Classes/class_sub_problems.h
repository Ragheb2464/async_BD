#ifndef CLASS_SP_H
#define CLASS_SP_H



#include "./Sub_Classes/class_opt.h"
#include "./Sub_Classes/class_sub_problem_manager.h"

#include "./Sub_Classes/class_data_S.h"
#include "./class_search_param.h"
//  ILOSTLBEGIN
#define Comm COMM_WORLD



class class_sub_problems
{ 
  private:
	  int iteration;
	  int *sc_worker, n_sc_worker;
	  int *aux_sc_worker, aux_n_sc_worker;
	  double tEndworker;
	  bool continue_running;
	  Data_S *data_S ; 		        //create the objects of the data class
	  Search_Param	*search_param;
	  
	  int aloc_memo, sc_id, sol_id, t, tt, ttt; 
	  double *G,  **core_y_SOL, *g_Values, ***Send_cut, *y_SOL, **X_Value;//,*Recv_data,
	  double Sub_Problem_Obj;
	  double wtime;
	  bool solve_sub_problems;//this for the cplex to say do not have warm start for the very beging sub-problem you solve
	  bool infSubProblem, infPapaSubProblem;		  
	  bool Phase_I;

	
  public:   
		class_sub_problems(int id, int n_sc, int n_workers, string instName, string scenarioName, string insClass)		
		{
				  search_param = new Search_Param(n_sc,n_workers);
				  
				  data_S = new Data_S(insClass); 		        //create the objects of the data class 								  
				  data_S->info(instName);		      //here we fill the network information	  	  				  
				  data_S->info_scenario(search_param, n_sc, scenarioName);   //this read the scenario file and fill the p array (probabilities)
				  data_S->read_assign(n_sc);	              //fill the arrays and vectors	  					
				  if( search_param->getData_Clean_UP_Moode() ) //clean up data	    		      			   		 
					data_S->clean_data(n_sc);	       				    
				  data_S->clean_memory(n_sc);                    //free some memory in the  data_S								  
				  // data_S->SceCreRetention(search_param, n_sc);
				  n_sc = n_sc - search_param->getMaster_sc(); // number of scenarios are total scenarios minus global scenarios (those who are not projected)
	 
				  iteration=0;
				  
				  MPI::Comm.Recv(&n_sc_worker, 1, MPI::INT, 0, search_param->getTagscnumber());				  
				  sc_worker = new int[n_sc_worker];
				  MPI::Comm.Recv(sc_worker, n_sc_worker, MPI::INT, 0, search_param->getTaglistsc());//receive the local scenarios	 
				  MPI::Comm.Recv(&aux_n_sc_worker, 1, MPI::INT, 0, search_param->getTagscnumber()+2);				  
				  aux_sc_worker = new int[aux_n_sc_worker];
				  MPI::Comm.Recv(aux_sc_worker, aux_n_sc_worker, MPI::INT, 0, search_param->getTaglistsc()+2);//receive the local scenarios	 
			    
				///set the demand values for artifical scenarios
				double prob_sum=0.0;
				for(int s=0; s < aux_n_sc_worker; s++)
					prob_sum += data_S->getP(aux_sc_worker[s]);
				// cout <<" artifical demands: " << endl;
				for(int k=0; k< data_S->getN_od(); k++){
					double sum=0.0;
					for(int s=0; s < aux_n_sc_worker; s++)
						sum += (data_S->getP(aux_sc_worker[s])/prob_sum) * fabs(data_S->getD(k,aux_sc_worker[s]));
					data_S->setD(k, n_sc + search_param->getMaster_sc(), sum);
					// cout << sum << "; ";
				}
			// cout <<endl;
				//reorder scenarios based on their demand sum
				// rearrangeScenarios(); 
				
 				    
				  MPI::Request Send_Cut_Req[n_sc_worker];
			            
			    
				//-------------------------------------------------------
				
				
				 if(search_param->getNum_cut_send() >= n_sc_worker)//maybe numbr of local scenarios are much smaller than the number of cuts we want to send
				    search_param->setNum_cut_send(n_sc_worker);
				  
				 
				 X_Value = new double*[data_S->getN_arcs()];
				 for(int a=0; a< data_S->getN_arcs(); a++)
				   X_Value[a] = new double[data_S->getN_od()]();
				 
				 
				  y_SOL = new double[data_S->getN_arcs()];
				  
				  G 	   = new double[data_S->getN_arcs() + data_S->getN_arcs()*data_S->getN_od()];// one from Papadakos SP + one from regular SP 
				  core_y_SOL = new double*[n_sc_worker];
				  for(int sce=0; sce <n_sc_worker; sce++ )
				     core_y_SOL[sce] = new double[data_S->getN_arcs()];
				  
				  g_Values = new double[n_sc];		//based on the generated cut for a scenario we generate valid cuts for other scenarios as well, so to do so, I need the g information calculated for each scenario based on that dual values; this array keeps that information 
				  //double aux_g_Values[search_param->getNum_cut_send() * n_sc];
				  //double aux_G_Values[search_param->getNum_cut_send() * n_sc * data_S->getN_arcs()];

				 
// 			         aloc_memo = 1 + n_sc_worker*(n_arcs + 3) + n_sc_worker * n_sc + n_sc_worker * n_arcs * n_sc; // 1 refers for sol_id; 3 for g , time and obj; 
												  //n_sc_worker * n_sc is for the g calculated for other scenarios becasue we use that information in adding cuts to the master problem in order to generate additional valid cuts
												  // n_sc_worker *n_arcs*n_sc for the  G[a] information whihc we will use to generate valid cuts for other scenarios (for each scenario sub-problem we G[a] for each scenarioÂ¼
// 				
				 // aloc_memo = 6 + 1 + data_S->getN_arcs() + data_S->getN_arcs() + 1;//(data_S->getN_arcs() + 5) +  n_sc +  data_S->getN_arcs() * n_sc + n_sc + 1;
				  aloc_memo = 1 + 6 + n_sc + data_S->getN_arcs() + data_S->getN_arcs()*data_S->getN_od();//1 because of art scenario
				  
// 				  num_cut_ready_to_send*(5+n_arcs) + num_cut_send * n_sc
// 				  cout << endl << "the aloc_memo of the worker " << id << " is: " << aloc_memo << endl;
				    
				  Send_cut = new double**[n_sc_worker+1];  //     
				  for(int s=0; s <n_sc_worker+1; s++)
				  {
				    Send_cut[s] = new double*[search_param->getN_interations()];
				    for(int sol_id=0; sol_id < search_param->getN_interations(); sol_id++)
				      Send_cut[s][sol_id] = new double[aloc_memo]();
				  }
				  //-------------------------------------------------------
				  vector <OptimalcutModel*> optModels;//creating the container of the object from the class of solving sub-problems
				  for(int s=0; s <n_sc_worker+1; s++){//+1 for artificial scenario
				    int sceID = sc_worker[s];
 				    optModels.push_back( new OptimalcutModel(data_S, search_param) );				  
				  }
				  
// 				  OptimalcutModel *optModels = new OptimalcutModel(data_S, search_param);//creating the object from the class of solving sub-problems
				  sub_manager *Manger = new sub_manager( data_S->getN_arcs());//creating the object from the class of solving sub-problems				  
				  
				  //*********************************************************************
				  for(int sce=0; sce <n_sc_worker; sce++ )
				  {
				      for(int a = 0; a < data_S->getN_arcs()- data_S->getN_od(); a++) //initialization of the M-W core point
				          core_y_SOL[sce][a] = search_param->getInitial_core_point();
				      //int dummy_arcs = n_arcs - n_od;//for the dummy arcs the y_sol and core_y_SOL are equal to 1 always, updating the core point should not change these values
				      for(int k = data_S->getN_arcs() - data_S->getN_od(); k < data_S->getN_arcs(); k++)		  
				          core_y_SOL[sce][k] = 1;		
				  }
// 				  Manger->setInitialY(data_S->getN_arcs(), n_sc_worker, core_y_SOL);
				  //*********************************************************************
				  wtime = MPI::Wtime();
				  //-------------------------------------------------------
				  continue_running = true;
				  Phase_I =true;//true if we are still at the first phase				  
				  
				  
				  MPI::Request  req = MPI::Comm.Irecv( &continue_running, 1, MPI::BOOL, 0, search_param->getTagconti_run() );	
				  		  				  
				  while(continue_running)
				  {	 
				    if (req.Test())
					break;
				      //******************************************************************
				      for(int s = 0; s < n_sc_worker; s++)
				      {
						sc_id = sc_worker[s];//this gives the id of scenario; i.e., actual_sc
						//******************************************************************
						
			 //receive info
						solve_sub_problems=false;
						if(Manger->Y_pool_manager(Phase_I, n_sc_worker, s, id, solve_sub_problems, y_SOL, sol_id, data_S->getN_arcs(), search_param))
							if(search_param->getCreate_art_SPs()){
								// cout << endl << "evaluate artificial scenario " << id << endl;
								infSubProblem	=false;
								infPapaSubProblem =false;
								// sc_id = n_sc+search_param->getMaster_sc();
								Sub_Problem_Obj = optModels[n_sc_worker]->solve(Phase_I, sc_worker, infSubProblem, infPapaSubProblem, search_param, data_S, X_Value, y_SOL, core_y_SOL, s, n_sc+search_param->getMaster_sc(), n_sc, G, g_Values);// givies the objective of the scenario multiple to its probability
								Manger->fill_send_vector(data_S->getN_od(), infPapaSubProblem, infSubProblem, Send_cut, sol_id, n_sc+id-1, g_Values, G, Sub_Problem_Obj, tEndworker, data_S->getN_arcs(), n_sc, n_sc_worker);//fill the Send_cut// 							
								// cout << endl << Sub_Problem_Obj <<  " evaluate artificial scenario " << id << endl;
								//sending to the master		  
								Send_Cut_Req[s] = MPI::Comm.Isend(Send_cut[n_sc_worker][sol_id], aloc_memo, MPI::DOUBLE, 0, n_sc*sol_id + n_sc + id);
								Send_Cut_Req[s].Free();
							}
						
						
						if(solve_sub_problems)//if we have new solution to use for this scenario
						{
							iteration++; // counter of workers iterations (how many messages have been sent)
							//updating the M-W core point for each sub-problem
							if(search_param->getUpdateCorePoint())
							  for(int a = 0; a < data_S->getN_arcs() - data_S->getN_od(); a++) //becuase I dont want to change the value of dummy arcs which is equal to 1
							    core_y_SOL[s][a] = (y_SOL[a] + core_y_SOL[s][a])/2;							
							//solving the sub-problem  
							infSubProblem	=false;
							infPapaSubProblem =false;
			   				wtime = MPI::Wtime();
							Sub_Problem_Obj = optModels[s]->solve(Phase_I, sc_worker, infSubProblem, infPapaSubProblem, search_param, data_S, X_Value, y_SOL, core_y_SOL, s, sc_id, n_sc, G, g_Values);// givies the objective of the scenario multiple to its probability
							tEndworker = MPI::Wtime() - wtime;
							// cout << endl << "time spend on SP: " << tEndworker << endl;
								// abort();  								  
							//******************************************************************
							Manger->fill_send_vector(data_S->getN_od(), infPapaSubProblem, infSubProblem, Send_cut, sol_id, sc_id, g_Values, G, Sub_Problem_Obj, tEndworker, data_S->getN_arcs(), n_sc, s);//fill the Send_cut
// 							if(id ==1)
// 							Manger->UB_Other_Sce(s, sol_id, Send_cut, y_SOL, X_Value, sc_id, Sub_Problem_Obj, n_sc, data_S);	
 							// abort();
							//sending to the master		  
							Send_Cut_Req[s] = MPI::Comm.Isend(Send_cut[s][sol_id], aloc_memo , MPI::DOUBLE, 0,  n_sc*sol_id+ sc_id);
 							Send_Cut_Req[s].Free();
							
	      //  		 			MPI::Comm.Send(Send_cut, aloc_memo , MPI::DOUBLE, 0, sol_id);	
//   							cout << endl << "send cut for solution id: " << sol_id << " for scenario: " << sc_id ;
// 							delete [] Send_cut[s][sol_id];
						}
				    
				      }
				    
				    
				  } //end while
				  
				  
				  delete [] Send_cut, G, g_Values, core_y_SOL; //Recv_data,
				  
				  
		}
////////////////////////////////////////////////////////////////////////////////////////////
		void rearrangeScenarios()
		{
			double sum;
			vector <pair <double,int>  > aux ;
			for(int s=0; s< n_sc_worker; s++){
				sum=0;
				for(int k=0; k < data_S->getN_od(); k++)
					sum += data_S->getD(k,sc_worker[s]);
				aux.push_back(make_pair(sum,sc_worker[s]));
			}
			sort(aux.begin(), aux.end(), PairCompare);//sort based on the double value in a non-decreasing order
			int t=0;
			for(int s=n_sc_worker-1; s>=0 ;s--)
				sc_worker[t++] = aux[s].second;// cout << aux[s].second << "  " << aux[s].first << endl; 
			
		}
////////////////////////////////////////////////////////////////////////////////////////////		
		~class_sub_problems()
		{

		}


};
#endif
