#ifndef CLASS_POOLS_MANAGER_H
#define CLASS_POOLS_MANAGER_H

#include "./CorePoint.h"
#include "./scenCluster.h"
ILOSTLBEGIN
#define Comm COMM_WORLD
class PoolsandMangers
{ 

  private: 
	
			  
		 //search variables and params
		 vector < vector <double> > *Y_Sol_Pool,  *v_cp, *v_cp_copy, *v_cp_art, *v_cp_copy_art,*v_sc_obj, *v_theta_sol; 		
		 vector <vector <int> > *v_workers_sc, *v_csw, clusterIDs;
		 vector <vector <int> > aux_v_workers_sc;
		 vector <vector <bool> > *v_css, /*//*v_ss, */*v_sw; 						
		 vector <double> *v_lb, *v_ub, v_aux;					
		 vector <int   >  *v_ciw, *v_nscs, *v_sid;
		 bool *v_iw;
		 int    sol_id, worker_id, Num_of_Equal_LB, aloc_memo;	
		 double ***Cuts_recv, GlobalUpperBound[2],  LB_best_time, UB_best_time, wtime, tEndLS; 	//I defined this to have the optimality gap as caluclated in the convergance_control function			  
		 double *Sol_send;
		 int *global_scenario_id;
		 int Potential_Cuts , identical_sol, identical_sol_counter, Num_newly_added_cuts;		//to check if the solution from current call of the master problem has been stored in the memory (true) or not (false);; how many equal lower bounds we will see in first phase before ending it
		 bool **Req_is_Made;
		 MPI::Request **req_rec_cut;
		 vector <int> ubImprovSolID;//id of solutions who improve the upper bound
		 double currentUB;// to get the upper bound on the recent solution used to generate cuts
		 double* yCore;
		 int numAwaitingCuts;//number of cuts expected to be delivered 		 
  public:
	PoolsandMangers(int n_workers, int n_sc, Data_S* data_S, Search_Param* search_param) 
	{	  
	  
			  identical_sol =0;
			  Potential_Cuts = 1e2;
			  numAwaitingCuts=0;
	  			  
			  Y_Sol_Pool 	= new vector< vector <double> >;
 			  v_theta_sol 	= new vector< vector <double> >;
			  v_cp 			= new vector< vector <double> >;
			  v_cp_copy		= new vector< vector <double> >;
			  v_cp_art 		= new vector< vector <double> >;
			  v_cp_copy_art	= new vector< vector <double> >;
			  v_sc_obj		= new vector< vector <double> >;			  
			  v_css 		= new vector< vector < bool > >;
			  v_sw 			= new vector< vector < bool > >;		  
			  v_csw 		= new vector< vector <int> >; 
			  v_workers_sc 	= new vector< vector <int> >;			  
			  v_ciw			= new vector<int>;
			  v_nscs      	= new vector<int>;
			  v_sid			= new vector<int>;			  
			  v_lb      	= new vector<double>;
			  v_ub 	 		= new vector<double>;				  		  
			  v_iw 			= new bool[n_workers];
			  			  
			  GlobalUpperBound[0] 	= 1e75;	//current upper bound
			  GlobalUpperBound[1] 	= 1e75;	//previous upper bound			  
			  Num_of_Equal_LB  		= 0.0;		//number of consequitive iterations that the lower bound has not changed 			  			  			  			  

			  Sol_send  = new double[data_S->getN_arcs() + n_sc  + 2 + 1]; // y sol and obja val, theta sol, Mp solving time ... 1 is for obj of the global scenarios
			  
			  yCore		= new double[data_S->getN_arcs()];
			  class_CorePoint *CorePoint = new class_CorePoint();
			CorePoint->Model(n_sc, search_param, data_S);
			CorePoint->Solve();
			// for(int sce=0; sce <n_sc_worker; sce++ )
			for(int a = 0; a < data_S->getN_arcs(); a++) //initialization of the M-W core point
			    yCore[a] = CorePoint->getySol(a);
			delete CorePoint;
			 /* for(int a=0; a<data_S->getN_arcs(); a++)
				yCore[a] = search_param->getInitial_core_point();*/
			  
			  // aloc_memo = 6 + 1 + data_S->getN_arcs() + data_S->getN_arcs() + 1;//(data_S->getN_arcs() + 5) + n_sc +  data_S->getN_arcs() * n_sc + n_sc + 1;//the last 1 is to indicate the solution was feasible (1) or infeasible (0)	
			  aloc_memo = 1 + 6 + n_sc + data_S->getN_arcs() + data_S->getN_arcs()*data_S->getN_od();//(data_S->getN_arcs() + 5) +  n_sc +  data_S->getN_arcs() * n_sc + n_sc + 1;
				  			  			  
 			  Req_is_Made =new bool*[n_sc+n_workers];
 			  for(int w_id=0; w_id < n_sc+n_workers; w_id++)
 			      Req_is_Made[w_id] =new bool[search_param->getN_interations()]();			  
			  Cuts_recv = new double**[n_sc+n_workers];			  
			  for(int s_id = 0; s_id < n_sc+n_workers; s_id++){
			     Cuts_recv[s_id] = new double*[search_param->getN_interations()];
			     for(int sol_id =0; sol_id < search_param->getN_interations(); sol_id++)
				  Cuts_recv[s_id][sol_id] = new double[aloc_memo];
			  }
			  req_rec_cut      = new MPI::Request*[n_sc+n_workers];
 			  for(int s_id=0; s_id < n_sc+n_workers; s_id++)
 			    req_rec_cut[s_id]  = new MPI::Request[search_param->getN_interations()];			  			  
	}
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	  void Cluster_Sce(Search_Param* search_param, Data_S *data_S, int n_workers, int n_sc, int master_sc)
	  {		    
	      vector <int> v_aux, v_aux1, gScID;
	      for(int s= n_sc -1; s >= n_sc - master_sc; s--)//we have put the global scenarios at the end of scenario vector
				gScID.push_back(s);
	      int numSc4Worker 	 = ( (n_sc - master_sc) - ((n_sc - master_sc)% n_workers) )/n_workers;	//number of scenarios to be send to the workers      
	      
	      for(int s=0; s<n_sc; s++)
			if( !(find(gScID.begin(), gScID.end(), s) != gScID.end()) ){		  
				  v_aux.push_back(s);
			  if(v_aux.size() >= numSc4Worker){
				(*v_workers_sc).push_back(v_aux);
				v_aux.clear();
			  }
		}
	      for(int s=0; s<v_aux.size(); s++)
		 (*v_workers_sc)[s].push_back(v_aux[s]);
	      //global scenarios
		  //if(gScID.size() > 0)
	      (*v_workers_sc).push_back(gScID);	


		ClusMainFunc *scenClus = new ClusMainFunc();
		scenClus->MainLoop(false, data_S, search_param,  n_sc- master_sc, n_workers);
		(*v_workers_sc).clear();
		for(int i=0; i< n_workers; i++){
			v_aux.clear();
			v_aux1.clear();
			for(int j=0; j < scenClus->GetClusterSize(i); j++)
				v_aux.push_back(scenClus->GetCluster(i, j));
			for(int j=0; j < scenClus->GetOriginalClusterSize(i); j++)
				v_aux1.push_back(scenClus->GetOrginalCluster(i, j));
			(*v_workers_sc).push_back(v_aux);
			(aux_v_workers_sc).push_back(v_aux1);
		}
		(*v_workers_sc).push_back(gScID);
		delete scenClus;
		
		
		// cout << endl << "EEEEEEE" << endl;
		ClusMainFunc *scenClus1 = new ClusMainFunc();
		scenClus1->Cleaning();
		// cout <<endl << "I am hereeeeeee\n";
		scenClus1->MainLoop(search_param->getCutAggr(), data_S, search_param,  n_sc- master_sc, search_param->getNumAggrCluster());
		// cout <<endl << "I am hereeeeeee\n";
		for(int i=0; i< search_param->getNumAggrCluster(); i++){
			v_aux.clear();
			for(int j=0; j < scenClus1->GetOriginalClusterSize(i); j++)
				v_aux.push_back(scenClus1->GetOrginalCluster(i, j));
			clusterIDs.push_back(v_aux);
		} 
		delete scenClus1;
		  
	  }
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////	  
	  void Assign_Sce_Proc(int n_workers, Search_Param *search_param)
	  {	    
	      global_scenario_id = new int[search_param->getMaster_sc()];
	      // sending group of scenarios for each worker  
	      for(int w_id = 0; w_id < n_workers ; w_id++) {	  
			  int n_sc_worker = (*v_workers_sc)[w_id].size();			//how many scenario
			  int aux_n_sc_worker = (aux_v_workers_sc)[w_id].size();
			  MPI::Comm.Isend(&n_sc_worker, 1, MPI::INT, w_id + 1, search_param->getTagscnumber()); //send it to the processor w_id + 1 (+1 because we dont want to send anything to the GSC which is processor number 0)
			  MPI::Comm.Isend(&aux_n_sc_worker, 1, MPI::INT, w_id + 1, search_param->getTagscnumber()+2); //send it to the processor w_id + 1 (+1 because we dont want to send anything to the GSC which is processor number 0)
			  int *sc_worker;
			  int *aux_sc_worker;
			  sc_worker = new int[n_sc_worker];
			  aux_sc_worker = new int[aux_n_sc_worker];
			  //int sc_worker [n_sc_worker];
	//    		  cout << endl << "I am worker-" << w_id+1 << " receiving scenario " ;
			  for(int i=0; i < n_sc_worker; i++){
				sc_worker[i] = (*v_workers_sc)[w_id][i];
				cout  << sc_worker[i] << " " ;
			  }
			  cout << endl;
			  for(int i=0; i < aux_n_sc_worker; i++)
				aux_sc_worker[i] = (aux_v_workers_sc)[w_id][i];
			  
			  MPI::Comm.Isend(sc_worker, n_sc_worker, MPI::INT, w_id + 1, search_param->getTaglistsc());
			  MPI::Comm.Isend(aux_sc_worker, aux_n_sc_worker, MPI::INT, w_id + 1, search_param->getTaglistsc()+2);
			  //cout << endl << "I am here" << endl;
			  delete [] sc_worker, aux_sc_worker;
	      }	
	       
	      //scenarios for the master problem			  
	      for(int i=0; i < (*v_workers_sc)[n_workers].size(); i++){		 
			  global_scenario_id[i] = (*v_workers_sc)[n_workers][i];
			  cout  << endl <<"I am the master problem and recieving global scenarios: " << global_scenario_id[i] << endl ;
	      }
	  }
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////	  
		bool Sol_Pool_Manager(bool M_D_Phase_I, Search_Param *search_param, Data_S *data_S, int master_sc, int n_sc, int n_workers, double *Sol_send, double GlobalLowerBound0, double GlobalLowerBound1)//vector <vector <double> > *v_sc_time, vector <vector <double> > *v_w_time, vector <double> *v_mp_time,
		{
			vector <bool> v_aux1; // checks that a worker used a solution or no	
			vector <double> v_aux2;
			int t;
			
			double fsobjval;
			bool record_sol = true;
			int pool_size = (*Y_Sol_Pool).size();
			
			//------ Save the opt solution in the solution pool
			 if ( GlobalLowerBound0 <= GlobalLowerBound1 + 0.5 ) // here I check if the lower bound is not improving, inccrease the indicator otherwise restart it from zero
				Num_of_Equal_LB++;
			 else		 
				Num_of_Equal_LB = 0; 
			 
			if( pool_size > 11  ){//here I check if we need to record the solution in the memory in accordance with the y variables	
				  bool sol_1=false, sol_2=false, sol_3=false, sol_4=false, sol_5=false;
				  record_sol = false;		      
				  for(int a = 0; a <data_S->getN_arcs() - data_S->getN_od(); a++){//checking with the last solution		      
					  if ( fabs(Sol_send[a] - (*Y_Sol_Pool)[pool_size-1][a]) > 0.01 )
						  sol_1 = true;			      
					  if ( fabs(Sol_send[a] - (*Y_Sol_Pool)[pool_size-2][a]) > 0.01 )
						  sol_2 = true;	
					if ( fabs(Sol_send[a] - (*Y_Sol_Pool)[pool_size-3][a]) > 0.01 )
						  sol_3 = true;	
					if ( fabs(Sol_send[a] - (*Y_Sol_Pool)[pool_size-4][a]) > 0.01 )
						  sol_4 = true;
					if ( fabs(Sol_send[a] - (*Y_Sol_Pool)[pool_size-5][a]) > 0.01 )
						  sol_5 = true;
					  if(sol_1 && sol_2  && sol_3  && sol_4  && sol_5){
						record_sol = true;
						break;
					  }
				 }	
			}
						
			
			if (record_sol){
				  identical_sol =0;
				  
				  v_aux2.clear();					
				  for(int a = 0; a< data_S->getN_arcs(); a++){	//Y SOL
					if(!M_D_Phase_I && pool_size < search_param->getWarmStartIterations() && search_param->getWarmStart()){//if we have warm start and we are still in the first phase
						Sol_send[a] = (Sol_send[a]+yCore[a])/2;
						yCore[a]=Sol_send[a];
					}
					v_aux2.push_back(Sol_send[a]);
				}
				  v_aux2.push_back(0);//I use this in the procedure of identifying useless solution(i.e., solve_choice_mas) .. 0 indicates the solution is still among useful solutions
				  (*Y_Sol_Pool).push_back(v_aux2);
				  //for the solution we have added to the memory since it is new, we dont have recvied the cuts and recourse obj of the sub problems associated to this solution
				  v_aux1.clear();
				  for(int s =0; s <n_sc; s++)
					v_aux1.push_back(false);
				  (*v_css).push_back(v_aux1); 		// cut indicator for each sol_id and Scenario (for this solution which scenarios has send back the cut)
							  
				  v_aux1.clear();
				  for(int w_id =0; w_id <n_workers; w_id++)
				v_aux1.push_back(false);
				  (*v_sw).push_back(v_aux1); 		// sol_id - workers indicator 		
							  
				  v_aux2.clear();
				  for(int s=0; s<n_sc;s++)
					v_aux2.push_back(1e10);
				  v_aux2.push_back(Sol_send[data_S->getN_arcs() + n_sc  + 1 + 1]);	//for the expected recourse cost of all global scenarios, to be used in the sol choice func
				  (*v_sc_obj).push_back(v_aux2);		// scenario obj for each sol_id		      
				  //--------------------------------------------------------------
				  //recourd theta values
				  v_aux2.clear();
				  for(IloInt s = 0; s < n_sc ; s++)
					v_aux2.push_back(Sol_send[data_S->getN_arcs() + 1 + s]);
				  (*v_theta_sol).push_back(v_aux2);				      		      
				  //------------------------------------------------------------
				  (*v_lb).push_back(Sol_send[data_S->getN_arcs()]);     // lower bound : master problem obj val
				  
				  fsobjval = 0.0; 				// first stage obj val
				  for(int a = 0; a < data_S->getN_arcs(); a++)
					fsobjval += Sol_send[a]*data_S->getF(a);
				  //cout << endl << "the first stage is " << fsobjval << endl;
				  fsobjval += Sol_send[data_S->getN_arcs() + n_sc  + 1 + 1];
				  //cout << endl << "the cost of global sc is " << Sol_send[n_arcs + n_sc  + 1 + 1] << endl;
				  
				  (*v_ub).push_back(fsobjval);		// the initial amount for ub is the fSobjval ()	
				  //cout << endl << "the number of solutions are: " << (*v_ub).size() << endl;
				  //cout << endl << " the UoOoOoer boOound is: " << fsobjval << endl;
				  
				  //--------------------------------------------------------------
				  (*v_nscs).push_back(0);  		//number of sc that returec cuts for each solution.	   
			}
			else{
			  identical_sol++;// =true;
			  cout << endl << "I did not record the first stage solution since it is already in there --> " << identical_sol << endl;		  
			}	
			return record_sol;
		}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	//this function picks the solutions based on the given strategy and send it to the available workers
	void Sol_Pool_controler(bool M_D_Phase_I, int n_workers, int n_sc, Data_S *data_S, Search_Param *search_param)
	{			
		int w_id, sol_id, n_sc_worker, sc_id;
		double   *Data_To_worker;
		vector<int>  v_aux;		
		
		MPI::Request Send_Sol;
		
		(*v_csw).clear(); 				// v_csw: current sol that are going to be sent to workers 
				
		 	// cout << endl << (*Y_Sol_Pool).size() << endl;
		for(int iw_id = 0; iw_id < n_workers; iw_id++){
			  v_aux.clear();	
			  for(int i = (*Y_Sol_Pool).size()-1; i >= 0; i--)			  
			      if((*v_sw)[i][iw_id] == false){//if this solution has not been send to the worker w_id			      
					  v_aux.push_back(i);//which solution
					  (*v_sw)[i][iw_id] = true;	//update the status of that solution ---> it has been send to the worker
					  // cout << endl << "picked solution is  " << i << endl;
					  break;					  
			      }
// 							  
			  if( v_aux.size() == 0 )
			    (*v_ciw).push_back(iw_id);//if no solution is found to be send add it to the completely idel list 
			    
			  //now send the solutions to the selected worker and update the memories 		
			  if ( v_aux.size() > 0 ) {
			     Data_To_worker = new double[v_aux.size()*(data_S->getN_arcs() + 1 + 1)];
			     w_id 	= iw_id;//(*v_csw)[i][0];	//this is for which worker (note it is recorder from 0, so 1 for example means the processor number 2)				   
			     sol_id    = v_aux[0];
			     (*v_sid).push_back(sol_id);
			     //------------------------------------------------------------------------------------
			     Data_To_worker[0] = sol_id;	//transfer the solution into an one dimensional vector 
		  // 	     cout << endl <<"I send message"<<endl;
			     for(int a = 0; a < data_S->getN_arcs(); a++)
				Data_To_worker[ a + 1] = (*Y_Sol_Pool)[sol_id][a];					
			     //if we have something for an idle worker ... first tell him to continue then give him the work piece
		// 	     cout << endl << "I am worker " << w_id << " recieiving solution id " << sol_id << endl;					
			     v_iw[w_id] = false;//since we are sendin it a solution, it will be occupid, (is it free, no --> = false)
										
			     n_sc_worker = (*v_workers_sc)[w_id].size();
			     //for(int s=0; s< n_sc_worker; s++)			//the scenarios in this worker will recieve the solution so, update the status to indicate that solution has been sent to those scenarios 
			     //{
			        //sc_id=(*v_workers_sc)[w_id][s];
// 				(*v_ss)[sol_id][sc_id] = true;  
			     //}
			     //indicate if the first phase is finished or not
			    if(!M_D_Phase_I)//is first phase still on
			       Data_To_worker[data_S->getN_arcs() + 1] = 1;
			    else
			       Data_To_worker[data_S->getN_arcs() + 1] = 0;
			    ///////
			    Send_Sol = MPI::Comm.Isend(Data_To_worker, data_S->getN_arcs() + 1 + 1, MPI::DOUBLE, w_id + 1, search_param->getTagysol());
			    Send_Sol.Free();//to free the handle
				
			    for(int sc_id =0; sc_id <n_sc ; sc_id++)
			      for(int s =0; s < n_sc_worker; s++)
					if( (*v_workers_sc)[w_id][s] == sc_id ){
						req_rec_cut[sc_id][sol_id] = MPI::Comm.Irecv(Cuts_recv[sc_id][sol_id], aloc_memo, MPI::DOUBLE, w_id + 1, n_sc*sol_id + sc_id);
						Req_is_Made[sc_id][sol_id] = true;// One is allowed to call MPI_TEST with a null or inactive request argument. In such a case the operation returns with flag = true and empty status. 
					}
				//cuts from artificial scenarios
				if(search_param->getCreate_art_SPs()){
					req_rec_cut[n_sc + w_id ][sol_id] = MPI::Comm.Irecv(Cuts_recv[n_sc + w_id][sol_id], aloc_memo, MPI::DOUBLE, w_id + 1, n_sc*sol_id + n_sc + w_id + 1);
					Req_is_Made[n_sc + w_id ][sol_id] = true;// One is allowed to call MPI_TEST with a null or inactive request argument. In such a case the operation returns with flag = true and empty status. 
				}
				delete[] Data_To_worker;
		    }				
	    }	  
	}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//this function waits for n_cut_wait workers to return their cut ... and add them to the cut pool
	int cut_pool_manager(bool Phase1, bool &stopping_ceriteria, int n_sc, int n_workers,  double GlobalLowerBound0, Data_S *data_S, Search_Param *search_param)//, vector< vector <double> > *v_sc_time  vector< vector <double> > *v_w_time,
	{	                                                                                             
		int t, sol_id, sc_id, n_sc_worker;
// 		vector<int> v_aux;//I use this to find out which worker have send back the cuts
		vector<double> v_aux, v_aux1;
		bool UpperBoundCheck;
		//---------------
		int  N_cut_wait =  search_param->getN_cut_wait();//how many cuts to wiat for

		//counting number of potential cuts to be received		    
		Potential_Cuts=0;
		for(int i=0; i< (*v_css).size(); i++)
		  for(int s=0; s < n_sc; s++)
		    if( Req_is_Made[s][i] )
				Potential_Cuts++;
		if(Potential_Cuts > 0 ){
			cout << "Num waiting cuts " << Potential_Cuts <<endl;
			numAwaitingCuts = Potential_Cuts;		
	//      	if(Num_of_Equal_LB[0] >=1 && Potential_Cuts >= 2*n_sc)
	//      		  n_cut_wait = n_sc  ;
			if(Potential_Cuts < search_param->getN_cut_wait())//if we have less cuts update N_cut_wait
				  N_cut_wait = Potential_Cuts;
					
					
			int test_count = 0;
			MPI::Status  status;
			/////First we check for N_cut_wait asscoiated to the last solution
			//if the MP has generated the same solution, we check for cuts associated to any solution 
			if(true){
			/* int */ sol_id = (*v_css).size() - 1;
				if(identical_sol >= 1){//if we are in the second phase and the algorithm fails to change the sol, wait for all the cuts associated to the previous solution
					Potential_Cuts =0;
					for(int s=0; s < n_sc; s++)
						if( Req_is_Made[s][sol_id] )
							Potential_Cuts++;
					N_cut_wait = Potential_Cuts;
					// cout << endl << "$$$$$$$$$$$$$$ we wait for cuts of the last solution: " << N_cut_wait << endl;
				}
				do
				{
					 // for(int sol_id=0; sol_id< (*v_css).size(); sol_id++)
						for(int s_id = 0; s_id < n_sc; s_id++)
							if(Req_is_Made[s_id][sol_id])
							   if( req_rec_cut[s_id][sol_id].Test() ) {				      					
					// 				    sol_id   = Cuts_recv[sc_id][sol_id][0];
									sc_id    = Cuts_recv[s_id][sol_id][1];
									if(s_id != sc_id)
										cout << endl << s_id << " - " << sc_id << endl;
									(*v_css)[sol_id][sc_id] = true;		//say the cut associate to solution sol_id is generate by scenario sc_id
									Req_is_Made[sc_id][sol_id] =false;
									//cut from Papadakos SP							
									v_aux.clear();
									v_aux.push_back(sol_id);  			//v_cp is sol_id, sc_id, g and G
									v_aux.push_back(sc_id);									
									v_aux.push_back(Cuts_recv[sc_id][sol_id][4 + sc_id]); 		// assign g																					
									// cout << "MP:: g for scen " << sc_id <<  " -> " << Cuts_recv[sc_id][sol_id][4 +sc_id] << endl;
									t=0;
									for(int a=0; a< data_S->getN_arcs(); a++){									      
										double aux=0.0;
										aux += Cuts_recv[sc_id][sol_id][4 + n_sc + 1 + t++];										
										for(int k = 0; k < data_S->getN_od(); k++)
											aux += data_S->getD(k,sc_id) * Cuts_recv[sc_id][sol_id][4 + n_sc + 1 + t++];
										v_aux.push_back(aux);	// assign G									
									}									
									v_aux.push_back(Cuts_recv[sc_id][sol_id][4 + n_sc + 1 + data_S->getN_arcs() +  data_S->getN_arcs()*data_S->getN_od()]);//the indicator of feasibility or infeasiblity of the used solution for regular SP
									v_aux.push_back(Cuts_recv[sc_id][sol_id][4 + n_sc + 1 + data_S->getN_arcs() +  data_S->getN_arcs()*data_S->getN_od()+1]);//the indicator of feasibility or infeasiblity of the used solution for the papdakos SP
									(*v_cp).push_back(v_aux);	//add it to the vector of constraints (it is a single cut plus the additional information to generate cut for other sencarios)																			
									//****************************************************************************																		
									if(Phase1 && search_param->getPropagatCut() && sc_id <= 60 )
										for(int ss=0; ss< n_sc; ss++)
											if(ss != sc_id) {
												v_aux.clear();
												v_aux.push_back(sol_id);  			//v_cp is sol_id, sc_id, g and G
												v_aux.push_back(ss);									
												v_aux.push_back(Cuts_recv[sc_id][sol_id][4 + ss]); 		// assign g																					
												// cout << "MP:: g for scen " << sc_id <<  " -> " << Cuts_recv[sc_id][sol_id][4 +sc_id] << endl;
												t=0;
												for(int a=0; a< data_S->getN_arcs(); a++){									      
													double aux=0.0;
													aux += Cuts_recv[sc_id][sol_id][4 + n_sc + 1 + t++];										
													for(int k = 0; k < data_S->getN_od(); k++)
														aux += data_S->getD(k,ss) * Cuts_recv[sc_id][sol_id][4 + n_sc + 1 + t++];
													v_aux.push_back(aux);	// assign G									
												}
												v_aux.push_back(0);//to keep number of iterations the cut has been in the pool
												(*v_cp_copy).push_back(v_aux);
											}									
									//****************************************************************************
									currentUB = (*v_ub)[sol_id];
									(*v_sc_obj)[sol_id][sc_id]  = Cuts_recv[sc_id][sol_id][2];
									for(int scen =0; scen <n_sc; scen++)							
										currentUB +=  (*v_sc_obj)[sol_id][scen];				
									//************check for updating the upper bound
									if ( currentUB < GlobalUpperBound[0] ){
										GlobalUpperBound[1] = GlobalUpperBound[0];//previous upper bound
										GlobalUpperBound[0] = currentUB;
										UB_best_time = MPI::Wtime() - wtime;
										// cout << endl << "This is solution id: " << sol_id << " and the its obj is: " << GlobalUpperBound[0] << endl;										
										if((fabs((GlobalUpperBound[0] -  GlobalLowerBound0)/(1 +  GlobalLowerBound0)) <= search_param->getTolerance()))
											  test_count=1e5;//if we converged do not wait for any further cuts
				// 					    if(ubImprovSolID.size() == 0)
				// 					      ubImprovSolID.push_back(sol_id);
										bool record_best_UB_id=true;
										for(int ii = ubImprovSolID.size()-1; ii>0; ii--)					    
										  if( ubImprovSolID[ii] == sol_id ){
											record_best_UB_id = false;
											break;
										  }					    
									   if(record_best_UB_id)
										 ubImprovSolID.push_back(sol_id);
									}																																			
									test_count++;								  
							   }				
				}while(test_count < N_cut_wait); // wait for n_wait_workers
				cout << endl << "num cuts received for sol id " << sol_id << " are " << test_count << endl;
			}
			///////check for the cuts from previous solutions
			if(test_count < Potential_Cuts ){
				 for(int sol_id=0; sol_id< (*v_css).size(); sol_id++)
					for(int s_id = 0; s_id < n_sc; s_id++)
						if(Req_is_Made[s_id][sol_id])
						   if( req_rec_cut[s_id][sol_id].Test() ){
									sc_id    = Cuts_recv[s_id][sol_id][1];
									if(s_id != sc_id)
										cout << endl << s_id << " - " << sc_id << endl;
									(*v_css)[sol_id][sc_id] = true;		//say the cut associate to solution sol_id is generate by scenario sc_id
									Req_is_Made[sc_id][sol_id] =false;
									//cut from Papadakos SP							
									v_aux.clear();
									v_aux.push_back(sol_id);  			//v_cp is sol_id, sc_id, g and G
									v_aux.push_back(sc_id);									
									v_aux.push_back(Cuts_recv[sc_id][sol_id][4+sc_id]); 		// assign g												
									t=0;
									for(int a=0; a< data_S->getN_arcs(); a++){									      
										double aux=0.0;
										aux += Cuts_recv[sc_id][sol_id][4 + n_sc + 1 + t++];										
										for(int k = 0; k < data_S->getN_od(); k++)
											aux += data_S->getD(k,sc_id) * Cuts_recv[sc_id][sol_id][4 + n_sc + 1 + t++];
										v_aux.push_back(aux);	// assign G									
									}									
									v_aux.push_back(Cuts_recv[sc_id][sol_id][4 + n_sc + 1 + data_S->getN_arcs() +  data_S->getN_arcs()*data_S->getN_od()]);//the indicator of feasibility or infeasiblity of the used solution for regular SP
									v_aux.push_back(Cuts_recv[sc_id][sol_id][4 + n_sc + 1 + data_S->getN_arcs() +  data_S->getN_arcs()*data_S->getN_od()+1]);//the indicator of feasibility or infeasiblity of the used solution for the papdakos SP
									(*v_cp).push_back(v_aux);	//add it to the vector of constraints (it is a single cut plus the additional information to generate cut for other sencarios)																			
									//****************************************************************************
									
									if(Phase1 && search_param->getPropagatCut() && sc_id <= 64)
										for(int ss=0; ss< n_sc; ss++)
											if(ss != sc_id) {
												v_aux.clear();
												v_aux.push_back(sol_id);  			//v_cp is sol_id, sc_id, g and G
												v_aux.push_back(ss);									
												v_aux.push_back(Cuts_recv[sc_id][sol_id][4 + ss]); 		// assign g																					
												// cout << "MP:: g for scen " << sc_id <<  " -> " << Cuts_recv[sc_id][sol_id][4 +sc_id] << endl;
												t=0;
												for(int a=0; a< data_S->getN_arcs(); a++){									      
													double aux=0.0;
													aux += Cuts_recv[sc_id][sol_id][4 + n_sc + 1 + t++];										
													for(int k = 0; k < data_S->getN_od(); k++)
														aux += data_S->getD(k,ss) * Cuts_recv[sc_id][sol_id][4 + n_sc + 1 + t++];
													v_aux.push_back(aux);	// assign G									
												}
												v_aux.push_back(0);//to keep number of iterations the cut has been in the pool											
												(*v_cp_copy).push_back(v_aux);
											}									
								//****************************************************************************
								currentUB = (*v_ub)[sol_id];
								(*v_sc_obj)[sol_id][sc_id]  = Cuts_recv[sc_id][sol_id][2];
								for(int scen =0; scen <n_sc; scen++)							
									currentUB +=  (*v_sc_obj)[sol_id][scen];				
								//************check for updating the upper bound
								if ( currentUB < GlobalUpperBound[0] ){
									GlobalUpperBound[1] = GlobalUpperBound[0];//previous upper bound
									GlobalUpperBound[0] = currentUB;
									UB_best_time = MPI::Wtime() - wtime;
									// cout << endl << "This is solution id: " << sol_id << " and the its obj is: " << GlobalUpperBound[0] << endl;										
									// if((fabs((GlobalUpperBound[0] -  GlobalLowerBound0)/(1 +  GlobalLowerBound0)) <= search_param->getTolerance()))
										  // test_count=1e5;//if we converged do not wait for any further cuts
									bool record_best_UB_id=true;
									for(int ii = ubImprovSolID.size()-1; ii>0; ii--)					    
									  if( ubImprovSolID[ii] == sol_id ){
										record_best_UB_id = false;
										break;
									  }					    
								   if(record_best_UB_id)
									 ubImprovSolID.push_back(sol_id);
								}																																		
								test_count++;								  
						   }				
			}				
			cout << endl << "num cuts received: " << test_count << endl;
			numAwaitingCuts -= test_count;
		}		
		return Potential_Cuts;
	}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	int auxcut_pool_manager(bool Phase1, int n_sc, int n_workers, double GlobalLowerBound0, Data_S *data_S, Search_Param *search_param)
	{	
		int test_count=0;
		if(numAwaitingCuts > 0){
			int Potential_Cuts=0;
			int t, sol_id, sc_id, n_sc_worker;
			vector<double> v_aux, v_aux1;
			bool UpperBoundCheck;
				for(int sol_id=0; sol_id< (*v_css).size(); sol_id++)
					for(int s_id = 0; s_id < n_sc; s_id++)
						if(Req_is_Made[s_id][sol_id])
						   if( req_rec_cut[s_id][sol_id].Test() ){	
									sc_id    = Cuts_recv[s_id][sol_id][1];
									if(s_id != sc_id)
										cout << endl << s_id << " - " << sc_id << endl;
									(*v_css)[sol_id][sc_id] = true;		//say the cut associate to solution sol_id is generate by scenario sc_id
									Req_is_Made[sc_id][sol_id] =false;
									//cut from Papadakos SP							
									v_aux.clear();
									v_aux.push_back(sol_id);  			//v_cp is sol_id, sc_id, g and G
									v_aux.push_back(sc_id);									
									v_aux.push_back(Cuts_recv[sc_id][sol_id][4+sc_id]); 		// assign g												
									t=0;
									for(int a=0; a< data_S->getN_arcs(); a++){									      
										double aux=0.0;
										aux += Cuts_recv[sc_id][sol_id][4 + n_sc + 1 + t++];										
										for(int k = 0; k < data_S->getN_od(); k++)
											aux += data_S->getD(k,sc_id) * Cuts_recv[sc_id][sol_id][4 + n_sc + 1 + t++];
										v_aux.push_back(aux);	// assign G									
									}							
									v_aux.push_back(Cuts_recv[sc_id][sol_id][4 + n_sc + 1 + data_S->getN_arcs() +  data_S->getN_arcs()*data_S->getN_od()]);//the indicator of feasibility or infeasiblity of the used solution for regular SP
									v_aux.push_back(Cuts_recv[sc_id][sol_id][4 + n_sc + 1 + data_S->getN_arcs() +  data_S->getN_arcs()*data_S->getN_od()+1]);//the indicator of feasibility or infeasiblity of the used solution for the papdakos SP
									(*v_cp).push_back(v_aux);	//add it to the vector of constraints (it is a single cut plus the additional information to generate cut for other sencarios)																			
									//****************************************************************************
									
									if(Phase1 && search_param->getPropagatCut() && sc_id <= 64)
										for(int ss=0; ss< n_sc; ss++)
											if(ss != sc_id) {
												v_aux.clear();
												v_aux.push_back(sol_id);  			//v_cp is sol_id, sc_id, g and G
												v_aux.push_back(ss);									
												v_aux.push_back(Cuts_recv[sc_id][sol_id][4 + ss]); 		// assign g																					
												// cout << "MP:: g for scen " << sc_id <<  " -> " << Cuts_recv[sc_id][sol_id][4 +sc_id] << endl;
												t=0;
												for(int a=0; a< data_S->getN_arcs(); a++){									      
													double aux=0.0;
													aux += Cuts_recv[sc_id][sol_id][4 + n_sc + 1 + t++];										
													for(int k = 0; k < data_S->getN_od(); k++)
														aux += data_S->getD(k,ss) * Cuts_recv[sc_id][sol_id][4 + n_sc + 1 + t++];
													v_aux.push_back(aux);	// assign G									
												}
												v_aux.push_back(0);//to keep number of iterations the cut has been in the pool																						
												(*v_cp_copy).push_back(v_aux);
											}
									
								//****************************************************************************
								currentUB = (*v_ub)[sol_id];
								(*v_sc_obj)[sol_id][sc_id]  = Cuts_recv[sc_id][sol_id][2];
								for(int scen =0; scen <n_sc; scen++)							
									currentUB +=  (*v_sc_obj)[sol_id][scen];				
								//************check for updating the upper bound
								if ( currentUB < GlobalUpperBound[0] ){
									GlobalUpperBound[1] = GlobalUpperBound[0];//previous upper bound
									GlobalUpperBound[0] = currentUB;
									UB_best_time = MPI::Wtime() - wtime;
									// cout << endl << "This is solution id: " << sol_id << " and the its obj is: " << GlobalUpperBound[0] << endl;										
									if((fabs((GlobalUpperBound[0] -  GlobalLowerBound0)/(1 +  GlobalLowerBound0)) <= search_param->getTolerance()))
										  test_count=1e5;//if we converged do not wait for any further cuts
									bool record_best_UB_id=true;
									for(int ii = ubImprovSolID.size()-1; ii>0; ii--)					    
									  if( ubImprovSolID[ii] == sol_id ){
										record_best_UB_id = false;
										break;
									  }					    
								   if(record_best_UB_id)
									 ubImprovSolID.push_back(sol_id);
								}																																		
								test_count++;								  
						   }
				cout << endl << ">>>>>>>num cuts received in the fractional node: " << test_count << endl;					   							
				numAwaitingCuts -= test_count;
		}
		return test_count;
	}
///////////////////////////////////////////////////////////////////////////////////////////////
	int art_cut_pool_manager(int n_sc, int n_workers, double GlobalLowerBound0, Data_S *data_S, Search_Param *search_param)
	{	
		int test_count=0;
		int Potential_Cuts=0;
		int t, sol_id, sc_id, n_sc_worker;
		vector<double> v_aux, v_aux1;
		bool UpperBoundCheck;
		sol_id=getY_Sol_Poolsize()-1;
		if(identical_sol == 0)
			while(test_count < n_workers){
				for(int w_id=0; w_id < n_workers; w_id++){
						int s_id =n_sc + w_id;//from this we can realize from which worker and scenario set the cut is associated to
						if(Req_is_Made[s_id][sol_id])
							 if( req_rec_cut[s_id][sol_id].Test() ){									
									sc_id    = Cuts_recv[s_id][sol_id][1];
									Req_is_Made[sc_id][sol_id] =false;
									//cut from Papadakos SP									
									v_aux.clear();
									v_aux.push_back(sol_id);  			//v_cp is sol_id, sc_id, g and G
									v_aux.push_back(sc_id);									
									v_aux.push_back(Cuts_recv[sc_id][sol_id][4+n_sc]); 		// assign g												
									t=0;
									for(int a=0; a< data_S->getN_arcs(); a++){									      
										double aux=0.0;
										aux += Cuts_recv[sc_id][sol_id][4 + n_sc + 1 + t++];										
										for(int k = 0; k < data_S->getN_od(); k++)
											aux += Cuts_recv[sc_id][sol_id][4 + n_sc + 1 + t++];
										v_aux.push_back(aux);	// assign G									
									}							
									v_aux.push_back(Cuts_recv[sc_id][sol_id][4 + n_sc + 1 + data_S->getN_arcs() +  data_S->getN_arcs()*data_S->getN_od()]);//the indicator of feasibility or infeasiblity of the used solution for regular SP
									v_aux.push_back(Cuts_recv[sc_id][sol_id][4 + n_sc + 1 + data_S->getN_arcs() +  data_S->getN_arcs()*data_S->getN_od()+1]);//the indicator of feasibility or infeasiblity of the used solution for the papdakos SP
									(*v_cp_art).push_back(v_aux);	//add it to the vector of constraints (it is a single cut plus the additional information to generate cut for other sencarios)																			
									//****************************************************************************
/*									
									
									//cut from regular SP
									v_aux1.clear();
									v_aux1.push_back(sol_id);
									v_aux1.push_back(sc_id);
									//for the g values
									v_aux1.push_back(Cuts_recv[sc_id][sol_id][5 + data_S->getN_arcs()]);
									//for G[a] values
									for (int a = 0; a <  data_S->getN_arcs() ; a++)			  
										v_aux1.push_back(Cuts_recv[sc_id][sol_id][6 + data_S->getN_arcs() + a]);											
									// (*v_cp_copy_art).push_back(v_aux1);	
*/
									test_count++;								  						 
							}
				}
			}
			cout << endl << "********num cuts received for artificial scenarios: " << test_count << endl;					   							
		
		return test_count;
	}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	bool intSol(Data_S *data_S, int sol_id)
	{
		bool solWasInt=true;		
		for(int a = 0; a < data_S->getN_arcs() - data_S->getN_od(); a++)
			if(getY_Sol_Pool(sol_id, a) > 1e-5 && getY_Sol_Pool(sol_id, a) < 1- 1e-5){
				solWasInt=false;
				break;
			}
		return solWasInt;
	}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	  //this is to control the convergence (stopping) of the whole algorithm    , 
	  void first_phase_ending(Search_Param *search_param, double &phase1Time, double &phase1Obj, double wtime, double GlobalLowerBound0, int Equal_LB, bool &transition_phase, bool &M_D_Phase_I, bool &stopping_ceriteria, int n_arcs, double tolerance, double *Sol_send)
	  {			    
		    //we want to finish first phase at the first time the lower bound we equal when we are coming bak from the transition phase
		    if (  ((Potential_Cuts == 0) || ( Num_of_Equal_LB  >= 10) || (fabs((GlobalUpperBound[0] - GlobalLowerBound0)/(1 + GlobalLowerBound0)) <= search_param->getTolerance()) || (MPI::Wtime() - wtime >= search_param->getTimeLimI()) )   )//;|| (fabs((GlobalUpperBound[0] - GlobalLowerBound[0])/(1 + GlobalLowerBound[0])) <= 0.001) || (fabs(GlobalUpperBound[0] - GlobalLowerBound[0]) <= 5)) ) //if first phase is done, move to the transition phase
		    {
				cout << endl << "________________________First-Phase is done_____________________" <<endl;
				 /*cout << endl <<  "The lower bound of: " << GlobalLowerBound[0] << endl;
				cout << endl << "The upper bound of: " << GlobalUpperBound[0] <<  endl;
				cout << endl << "The time consumption in seconds: " << MPI::Wtime() - wtime<< endl; 
				cout <<  "_______________________________________________________________"<<endl;*/							
				// cout << endl <<search_param->getTolerance() <<  " I am stopping phase 1 because: " << fabs((GlobalUpperBound[0] - GlobalLowerBound0)/(1 + GlobalLowerBound0))  << endl;
				
				phase1Time = MPI::Wtime() - wtime;
				phase1Obj  = GlobalLowerBound0;
				//if first phase is optimal and integer, we are done
				bool first_phase_global_optimal = true;
				for(int a =0; a <n_arcs; a++)
				  if( Sol_send[a] > 0.01 && Sol_send[a] < 0.99 )//if on arcs has taken non intger value the first phase is not optimal				    
				    first_phase_global_optimal = false;					  
				 				
				M_D_Phase_I = true;//ok let say the first phase is done (true)
					
				if(!first_phase_global_optimal || Equal_LB == 1)					    
				    Num_of_Equal_LB  = 0;			//as these two (i.e, number of iterations that the lower bound has not increased and the global upper bound) parameters will be used again in the second phase to check for the stopping criteria we will reset them		      				
				else
				  if((fabs((GlobalUpperBound[0] - GlobalLowerBound0)/(1 + GlobalLowerBound0)) <= 0.005))
					{
					  cout << endl << "the problem ending because the first phase is global optimal" << endl;
					  stopping_ceriteria =true;
					}	
				//reset the upper bounds if we keep going with algorithm 
				if(!stopping_ceriteria){
					GlobalUpperBound[0] = 1e75;
					GlobalUpperBound[1] = 1e75;
				}					
		  }		
	  }
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	      int getglobal_scenario_id(int s) 					{ return global_scenario_id[s];  }
	      int getY_Sol_Poolsize() 							{return (*Y_Sol_Pool).size();    }	      
	      double getY_Sol_Pool(int sol_id, int a) 			{return (*Y_Sol_Pool)[sol_id][a];  }
	      double getV_theta_sol(int sol_id, int s)			{return (*v_theta_sol)[sol_id][s];}		
	      void setY_Sol_Pool(int i, int j, int k)  			{(*Y_Sol_Pool)[i][j] = k; }
	      double getv_cp(int i, int j) 						{return (*v_cp)[i][j]; }	     	     
	      int getv_cpSize()  								{return (*v_cp).size();  }	     
	      double getv_cp_copy(int i, int j) 				{return (*v_cp_copy)[i][j];    }	      
	      int getv_cp_copySize()  							{return (*v_cp_copy).size(); }
		  double getv_cp_art(int i, int j) 					{return (*v_cp_art)[i][j]; }	     	     
	      int getv_cp_artSize()  							{return (*v_cp_art).size();  }	     
	      double getv_cp_copy_art(int i, int j) 			{return (*v_cp_copy_art)[i][j];    }	      
	      int getv_cp_copy_artSize()  						{return (*v_cp_copy_art).size(); }
		  int get_worker_sc_ids(int w_id, int scen_id) 		{return (*v_workers_sc)[w_id][scen_id];}
		  int get_worker_sc_size(int w_id) 					{return (*v_workers_sc)[w_id].size();}
		  int get_aux_worker_sc_ids(int w_id, int scen_id) 	{return (aux_v_workers_sc)[w_id][scen_id];}
		  int get_aux_worker_sc_size(int w_id) 				{return (aux_v_workers_sc)[w_id].size();}
		  bool getv_css(int i, int s) 						{return (*v_css)[i][s];  }
	      void setV_css(int s, int i, bool state) 			{(*v_css)[i][s] = state;  }	      
		  void setReq_is_Made(int s,int i,bool state)		{Req_is_Made[s][i] =state;}	      
	      double getv_sc_obj(int i, int s) 					{return (*v_sc_obj)[i][s]; }
	      double getGlobalUpperBound(int i)					{return GlobalUpperBound[i];}	     
	      void setGlobalUpperBound(int i,double value)		{GlobalUpperBound[i] = value;}
	      int getPotential_Cuts()							{return Potential_Cuts;	}	     
	      int getUbImprovSolIDSize() 						{return ubImprovSolID.size();}	     
	      int getUbImprovSolID(int i)						{return ubImprovSolID[i];}	     
		  double getCurrentUB() 							{return currentUB;}
		
		
		
		int getClusterSize(int clusterID) {return clusterIDs[clusterID].size();}
		int getClusterIDs(int clusterID, int element) {return clusterIDs[clusterID][element];}
	     
		
		
		void Cleanv_cp()  			{delete v_cp; v_cp= new vector< vector <double> >;}	     
	    void Cleanv_cp_copy()		{delete v_cp_copy;  v_cp_copy = new vector< vector <double> >;}
	    void Cleanv_cp_copyID(int i){v_cp_copy->erase(v_cp_copy->begin()+i);}
		void Cleanv_cp_art()  		{delete v_cp_art;	v_cp_art = new vector< vector <double> >; }	     
	    void Cleanv_cp_copy_art()	{delete v_cp_copy_art;  v_cp_copy_art = new vector< vector <double> >;    }	    
	    void CleanY_Sol_Pool()		{delete Y_Sol_Pool;   Y_Sol_Pool = new vector< vector <double> >;	    }	    
	    void Cleanv_ub()			{delete  v_ub;  v_ub = new vector<double>();	    }
	    void Cleanv_sc_obj()		{delete v_sc_obj;v_sc_obj = new vector< vector <double> >;}

/////////////////////////////////////////////////////////////////////////////////////////////////


	~PoolsandMangers()
	{

	}


};
#endif
