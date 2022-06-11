 #ifndef CLASS_PARA_H
#define CLASS_PARA_H


// ILOSTLBEGIN
#define Comm COMM_WORLD  
class class_GSC_Param
{
	  private:
// 			  PoolsandMangers *Managers  = new PoolsandMangers();
// 			  MasterSolChoice *Sol_Choice;
			  
			  //i need to define these parameters as global variables to be used in callback
			  
			  //all the input variables to the GSC class
			  
			  int Tag_Useless_Sol,  TagNumYsol,  num_cut_send, **arcs, **od,  n_od, n_workers, iteration, n_interations, n_p, master_sc, n_nodes, n_arcs, n_sc, Tagysol, Tagconti_run, Tagcuts, Tagscnumber, Taglistsc, Tagms_send, Tagms_rec, Tagms_ncsend, Tagms_csend, Tagms_run, Tagms_phase;
			  double *c,  **d,  *u, tolerance,  * f, * p;
			  bool M_D_Phase_I,  continue_running;
			  MPI::Status *status_cuts, status_MPsol;
			  MPI::Request **req_rec_cut, requestrec_MPsol;
			  
			  //search variables and params
			  vector <vector <double> > *Y_Sol_Pool, /*//*v_sptheta,*/ *v_cp, *v_cp_copy, *v_sc_obj, * v_cp1; 		
				 										
			  vector <vector <int> > *v_workers_sc, *v_csw;
			  
			  vector <vector <bool> > *v_css, *v_ss, *v_sw; 						
															
			  vector <double> *v_lb, *v_ub, v_aux;					
															
			  vector <int   >  *v_ciw, *v_nscs, *v_sid, *con_copy_binding, *con_binding, *con_copy_iter, *con_iter;
			  bool *v_iw, *w_cut_request, *y_fixed;
			  
			  int    sol_id, worker_id, aloc_memo, *Data_SMP, v_cpl, solsend_strategy, comm_strategy, cut_strategy;	
			  
			  double ***Cuts_recv, *Sol_get, *GlobalUpperBound, *GlobalLowerBound, *Num_of_Equal_LB, LB_best_time, UB_best_time, wtime, tEndLS,first_phase_time, first_phase_obj; 	//I defined this to have the optimality gap as caluclated in the convergance_control function			  
			  
			  
			  
			  double *Sol_send, *Cuts_p_rec;
			   
			  int n_rec_cuts,  *Data_recv, n_cuts;
			  bool **Req_is_Made;
// 			  bool *added_theta;
			  int *global_scenario_id;
			  
			  MPI::Request *requestsend_ysol;
			  double m_objval, tStart, tEnd, Recourse_Obj = 0; //the second stage obj associated to the global scenarios
			  double n_cut_wait;
			  int Potential_Cuts , identical_sol, identical_sol_counter, Equal_LB, Num_newly_added_cuts;		//to check if the solution from current call of the master problem has been stored in the memory (true) or not (false);; how many equal lower bounds we will see in first phase before ending it
			  double Best_Lower_Bound ; 							//to have the best lower bound (and its time)
			  bool Strong_Ineq, stopping_ceriteria, transition_phase, second_phase, end_transition_phase;
 			 
			 		
			  const int num_copy_cut_add=1;
			  
			  //cplex related paramters and variables
// 			  IloEnv env;		//environment 
// 			  IloModel model;	//
// 			  IloCplex cplex;
// 			  IloRangeArray *con, *con_optcut, *con_copy, *con_heuristic;
			  
			  
			  IloNumVarArray3 x;
			  IloNumVarArray  y;
			  IloNumVarArray  theta;
			  
			  
// 			  IloObjective obj;
			  
			  /*
			  double ac_inac_p;
			  int 	 It_devision, optcut_threshold;
			  int 	 Mas_iteration, total_size_of_con_copy;
			  
			  
			  */
			    
			 IloNumArray y_SOL;
			 IloNumArray theta_SOL;
			 IloNumArray3 x_Sol;  
			 IloNumArray best_y;
			 IloNumArray best_theta;
			 IloNumArray3 best_x;
			 IloNumArray y_Incum;
			  
 	  public:
			 class_GSC_Param()
			  {		    

			  }
		
			  void InitializeParam(int Tag_Useless_Sol0, int TagNumYsol0, int num_cut_send0, int**arcs0, double *c0, double **d0, double *u0, int **od0, int n_od0, double tolerance0, bool M_D_Phase_I0, bool continue_running0, int n_workers0, int iteration0, int n_interations0, int n_p0, int master_sc0, int n_nodes0, int n_arcs0, int n_sc0, double* f0, double* p0, int Tagysol0, int Tagconti_run0, int Tagcuts0, int Tagscnumber0, int Taglistsc0, int Tagms_send0, int Tagms_rec0, int Tagms_ncsend0, int Tagms_csend0, int Tagms_run0, int Tagms_phase0/*, MPI::Status *status_cuts0, MPI::Status status_MPsol0, MPI::Request **req_rec_cut0,  MPI::Request requestrec_MPsol0*/)
			  {
			    	Tag_Useless_Sol=Tag_Useless_Sol0;  TagNumYsol=TagNumYsol0;  num_cut_send=num_cut_send0; arcs=arcs0; od=od0;  n_od=n_od0; n_workers=n_workers0; iteration=iteration0; n_interations=n_interations0; n_p=n_p0; master_sc=master_sc0; n_nodes=n_nodes0; n_arcs=n_arcs0; n_sc=n_sc0; Tagysol=Tagysol0; Tagconti_run=Tagconti_run0; Tagcuts=Tagcuts0; Tagscnumber=Tagscnumber0; Taglistsc=Taglistsc0; Tagms_send=Tagms_send0; Tagms_rec=Tagms_rec0; Tagms_ncsend=Tagms_ncsend0; Tagms_csend=Tagms_csend0; Tagms_run=Tagms_run0; Tagms_phase=Tagms_phase0;  c=c0; d=d0; u=u0; tolerance=tolerance0; f=f0;  p=p0;  M_D_Phase_I=M_D_Phase_I0;  continue_running=continue_running0;  /*status_cuts=status_cuts0; status_MPsol=status_MPsol0;			  req_rec_cut=req_rec_cut0; requestrec_MPsol=requestrec_MPsol0;*/

			  }
			
		
			~class_GSC_Param()
			  {		  

			  }
};//end class
#endif