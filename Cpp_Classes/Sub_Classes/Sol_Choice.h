#ifndef CLASS_MAS_SOL_H
#define CLASS_MAS_SOL_H


#include <ilcplex/ilocplex.h>
// ILOSTLBEGIN
#define Comm COMM_WORLD


typedef IloArray<IloNumVarArray> IloNumVarArray2;
typedef IloArray<IloNumVarArray2> IloNumVarArray3;
typedef IloArray<IloNumVarArray3> IloNumVarArray4;



class MasterSolChoice
{ 
	private:

		IloEnv env;		//environment 
		IloModel model;	//
		IloCplex cplex;
		IloRangeArray  *con_optcut, *con_copy, *con_heuristic;
		
		IloNumArray fixed_y;
		
		IloNumVarArray  y;
		IloNumVarArray  theta;
		
		
		IloObjective obj;
		int 	 total_size_of_con_copy ;
		
		
		
        public:
		MasterSolChoice(Data_S *data_S, Search_Param *search_param, int n_sc) //constructor
		{
		        

			model 		=  IloModel(env);
			cplex 		=  IloCplex(model);
			con_optcut   	= new IloRangeArray(env);
			con_copy   	= new IloRangeArray(env);
			con_heuristic   = new IloRangeArray(env);
			obj   		= IloMinimize(env);
			
			
			fixed_y= IloNumArray(env, data_S->getN_arcs());
			
			//***the master only includes y and theta variables****
			y = IloNumVarArray(env, data_S->getN_arcs(), 0, 1);	//ILOINT defining an variable for each arc (on dimensional variable)
			this->y=y;
			
			
			
			
			theta = IloNumVarArray(env, n_sc + search_param->getMaster_sc(), 0, IloInfinity);	//defining variables to approximate the recourse costs
											//I generate a theta for the global scenarios too, but keep their coeficient to zero (just to ease the programming pain :P)
			this->theta = theta;

			(model).add(y); // we dont need to add variables explicitly to the model as they will implicitly be added to the model through the constraints
			(model).add(theta);

			
			
			// ***** OBJECTIVE FUNCTION****
			
			int ss;
			IloExpr expr(env);
			for(int a = 0; a < data_S->getN_arcs(); a++)			
			  expr += data_S->getF(a)*y[a];
			
			
			 
			  for(int s = 0; s < n_sc ; s++)			  
			      expr += data_S->getP(s)* theta[s];//initially I need all theta not be in obj otherwise the master problem runns outbounded
			
			obj.setExpr(expr);

			(model).add(obj);

			expr.end();
			
			
			
			//---------------------solve ------------------------------
			
			(cplex).setParam(IloCplex::TiLim, 3600);
 			(cplex).setParam(IloCplex::ClockType, 1);//to set it equal to cpu time since it is seq
			(cplex).setParam(IloCplex::SimDisplay, 0);			
 			(cplex).setParam(IloCplex::Threads, 1);		

	
		}
		
		
		
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		
//identify which of the solutions in the pool are no longer needed

		void solve_choice_mas(int n_workers,  int n_sc, PoolsandMangers *Managers, Data_S *data_S, Search_Param *search_param)//, vector<int> *v_ac, vector<int> *v_ac_copy, vector<int> *v_inac, vector<int> *v_inac_copy, vector<int> *v_acinac_aux ) 
		{                             
			
			double m_objval, tStart, tEnd;//, *slackvals;
			double Recourse_Obj = 0;			//the second stage obj associated to the global scenarios
			vector <int> useless_y_sol;
			
			
			
			for(int i =0; i <  Managers->getY_Sol_Poolsize(); i++)
			  if(Managers->getY_Sol_Pool(i, data_S->getN_arcs())/*(*Y_Sol_Pool)[i][n_arcs]*/ < 1)
			  for(int s =0; s < n_sc; s++)
			    if( !Managers->getv_css(i, s)/*(*v_css)[i][s]*/ )//if there is at least one sub that has not send back cut for this scenario, we check the validity of the solution
			    {
			      for(int a=0; a < data_S->getN_arcs(); a++)
				fixed_y[a] = Managers->getY_Sol_Pool(i, a)/*(*Y_Sol_Pool)[i][a]*/;
			      
				  y.setBounds( fixed_y, fixed_y);
				  
				  if( !(cplex).solve() ) 
				    {			  
					cout << endl << "its weird it could solve the sol choice master" << endl;
				    }
				  else
				  {
				      m_objval = (cplex).getObjValue() + Managers->getv_sc_obj(i, n_sc);//(*v_sc_obj)[i][n_sc];//I have the recourse cost of global scenarios in this vector location in order to not write them in this auxulary model and take time
//  				      cout << endl << "!!!!!!!: " << m_objval << endl;
				      if(m_objval >= Managers->getGlobalUpperBound(0) )
				      {
					
					cout << endl << "solution id " << i << " is not useful any more" << endl;
					Managers->setY_Sol_Pool( i, data_S->getN_arcs(), 10);//(*Y_Sol_Pool)[i][n_arcs] = 10; //saying this solution is among useless ones and dont check it next time
					useless_y_sol.push_back(i);
					for(int ss =0; ss < n_sc; ss++)
					{
// 					  cout << "theta-" << ss << "  --> " << (cplex).getValue(theta[ss]) << endl;
// 					  (*v_css)[i][ss] = true;
					  Managers->setReq_is_Made(ss, i, false);//   V_css[i][ss] = true; its cut has been received
					}
				      }
				  }
				
				
				break;
			    }
			  
			  
			  int num_useless = useless_y_sol.size();
			  if(num_useless >0)
			  {
			      MPI::Request Usless_sol;
			      int *usless_send;
			      usless_send = new int[num_useless];
			      for(int i=0; i<num_useless; i++)
				usless_send[i] = useless_y_sol[i];
			      for(int w_id=0; w_id< n_workers; w_id++)
			      {
				Usless_sol = MPI::Comm.Isend(usless_send, num_useless, MPI::INT, w_id + 1, search_param->getTag_Useless_Sol());
				Usless_sol.Free();
			      }
			      
			      delete [] usless_send;
			  }
		}


//to cleanup the first phase solutions from the workers memory
		void clean_phase_I_sol(PoolsandMangers *Managers,  int n_workers,  int n_arcs, int n_sc, int Tag_Useless_Sol)//, vector<int> *v_ac, vector<int> *v_ac_copy, vector<int> *v_inac, vector<int> *v_inac_copy, vector<int> *v_acinac_aux ) 
		{                             
			
			double m_objval, tStart, tEnd;//, *slackvals;
			double Recourse_Obj = 0;			//the second stage obj associated to the global scenarios
			vector <int> useless_y_sol;
			
			
			
			for(int i =0; i <  Managers->getY_Sol_Poolsize(); i++)
			  if(Managers->getY_Sol_Pool(i,n_arcs)/*(*Y_Sol_Pool)[i][n_arcs]*/ < 1)
			  {
					cout << endl << "solution id " << i << " will be removed from the solution pool" << endl;
					Managers->setY_Sol_Pool(i,n_arcs, 10);//(*Y_Sol_Pool)[i][n_arcs] = 10; //saying this solution is among useless ones and dont check it next time
					useless_y_sol.push_back(i);
					for(int ss =0; ss < n_sc; ss++)
					  Managers->setReq_is_Made(ss, i, false);//  V_css[i][ss] = true; its cut has been received
					
			  }
			    
			    
			
			  
			  int num_useless = useless_y_sol.size();
			  if(num_useless >0)
			  {
			      MPI::Request Usless_sol;
			      int *usless_send;
			      usless_send = new int[num_useless];
			      for(int i=0; i<num_useless; i++)
				usless_send[i] = useless_y_sol[i];
			      for(int w_id=0; w_id< n_workers; w_id++)
			      {
				Usless_sol = MPI::Comm.Isend(usless_send, num_useless, MPI::INT, w_id + 1, Tag_Useless_Sol);
				Usless_sol.Free();
			      }
			      
			      delete [] usless_send;
			  }
		}


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	
		  //"NOTE: we only have optimality cuts
		void regular_cuts(int n_arcs, int n_sc, PoolsandMangers *Managers) 
		{		  
			
			double  g;
			int sc_id, sol_id;
			


			
			
			for(int i = 0; i < Managers->getv_cpSize()/*(*v_cp).size()*/; i++)
			{
			      sol_id = Managers->getv_cp(i, 0);//(*v_cp)[i][0];	//solution ID
			      sc_id  = Managers->getv_cp(i, 1);;//(*v_cp)[i][1];	//this indicates which scenario this cut belongs to
			      g	     = Managers->getv_cp(i, 2);//(*v_cp)[i][2];
			      
			    
 			   
			      
			      IloExpr expr(env);
			      
			      for(int a = 0; a < n_arcs; a++)
				  expr +=Managers->getv_cp(i, a+3) /*(*v_cp)[i][a+3]*/*y[a];
				    
			      expr += g;
			      expr += -theta[sc_id];			    			   			     
			       
			     			      
			      (*con_optcut).add(IloRange( env, -IloInfinity, expr, 0));//, con_name));	
				
			      expr.end();
			      
			}
			

			(model).add((*con_optcut));
		 			 			
			
			
		}
			
			
			
			
		void copied_cuts( double *Sol_send,  int n_sc, PoolsandMangers *Managers, Data_S *data_S, Search_Param *search_param) 
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
					    if (Managers->getglobal_scenario_id(ss)  == CopyCut )
						add_cut = false;

					  
					  
					    if ( add_cut) //i dont want to add cut for the global scenarios and the current scenairo				  
					    {
						for(int j =0; j <search_param->getNum_copy_cut_add(); j++)
						  cut_to_add_violation[j] = -1e5;
						
						for (int i =0; i <Managers->getv_cp_copySize()/*(*v_cp_copy).size()*/; i++)//here we determine which cuts to add
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
							
						  }
						  
							  //*************************************
						  for(int j=0; j<search_param->getNum_copy_cut_add(); j++)
						      if( cut_to_add_violation[j] >= 10)						 
							  {		
								id = cut_to_add_id[j];
	  // 						      cout << endl << cut_to_add_violation[j];
								IloExpr expr(env);
								expr.clear();
								for(int a = 0; a < data_S->getN_arcs(); a++)
								    expr += Managers->getv_cp_copy(id, 2 + n_sc + CopyCut*data_S->getN_arcs() + a) /*(*v_cp_copy)[id][2 + n_sc + CopyCut*n_arcs + a]*/ * y[a];
							      
								expr += Managers->getv_cp_copy(id, 2 + CopyCut);//(*v_cp_copy)[id][2 + CopyCut];
								
								expr -= theta[CopyCut];
								
 								(*con_copy).add(IloRange( env, -IloInfinity, expr, 0));//, con_name));
								expr.end();
								
								
							  }
						  
						  
					  }//if
					//******************************** end of copy the cut for other scenarios *******************
				  }//for
				  





			//(model).remove((*con_copy));
// 			if( (*con_copy).getSize() - total_size_of_con_copy > 0)
			{
  			  (model).add((*con_copy));
			  total_size_of_con_copy = (*con_copy).getSize();
//  			  cout << endl << "The size of copied Bender cuts are: " <<  (*con_copy).getSize()  << endl;
			}
						 			 			
			
			
			
			
		
		}
				


	


		~MasterSolChoice()
		{
		  
		  
			// delete for new // means pointers
			//delete model;
			//delete cplex;
			delete con_optcut;
			delete con_copy;
			
			
			
			(cplex).end();
			(model).end();
			
			
			(*con_optcut).end();
			(*con_copy).end();
			//(*temp_con).end();
			obj.end();
			env.end();

		}


};
#endif
