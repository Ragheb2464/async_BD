#ifndef CLASS_SUB_OPT_H
#define CLASS_SUB_OPT_H


// ILOSTLBEGIN

#define Comm COMM_WORLD


class sub_manager
{ 
  private:
		
		int count;
		bool sol_has_rec;
		MPI::Request  *Rec_Sol;
		double *Recv_data;
		MPI::Status  Useless_Sol_Status;
		vector <vector <double> > Y_Pool;		          

  public:
		//create the sub-problems
		sub_manager(int n_arcs) 
		{	
		   sol_has_rec=true;
		   Rec_Sol  = new MPI::Request[1];
		   Recv_data = new double [n_arcs + 1 + 1];
		}
//////////////////////////////////////////////////////////////////////////////////////	
		void setInitialY( int n_arcs, int n_sc_worker, double **core_y_SOL)
		{		          
			  vector < double > aux;
// 			  aux.clear();
			  aux.push_back(-1);
			  for (int a =1; a < n_arcs + 1; a++)
			      aux.push_back(core_y_SOL[0][a]);
			  for (int i =0; i < n_sc_worker ; i++)
			      aux.push_back(0);//to indicate if sol is used for local scenario s or not
			  Y_Pool.push_back(aux);//sol_id, y, indicating if it is used for which of local scenarios
// 			  aux.clear();	
		}
///////////////////////////////////////////////////////////////////////////////////////
		bool Y_pool_manager(bool& Phase_I, int n_sc_worker, int s, int id,bool& solve_sub_problems, double *y_SOL, int& sol_id, int n_arcs, Search_Param* search_param) 
		{
		      //check for new incomming solutions from master
		      if(sol_has_rec){
				Rec_Sol[0] = MPI::Comm.Irecv( Recv_data, n_arcs + 1 + 1, MPI::DOUBLE, 0, search_param->getTagysol() );
				sol_has_rec=false;
		      }
		      
		      
		      if( Rec_Sol[0].Test() ){//MPI::Comm.Iprobe(0, Tagysol, status_Y_sol) )		      		      
				  sol_has_rec =true;
				  vector < double > aux;
	// 			  aux.clear();
				  for (int i =0; i < n_arcs + 1; i++)
					  aux.push_back(Recv_data[i]);
				  for (int i =0; i < n_sc_worker ; i++)
					  aux.push_back(0);//to indicate if sol is used for local scenario s or not
				  Y_Pool.push_back(aux);//sol_id, y, indicating if it is used for which of local scenarios
	// 			  aux.clear();	
				  if(Recv_data[n_arcs + 1] < 0.5)
					Phase_I=false;//first phase is finished
		      }
		      
		      
		      //check for useless solutions to be dumped
		      if(MPI::Comm.Iprobe(0, search_param->getTag_Useless_Sol(), Useless_Sol_Status)){		      
	// 			MPI_Get_count(&status, MPI_INT, &number_amount);
				int Useless_Sol_Size = Useless_Sol_Status.Get_count(MPI::INT);
				
				int *Usless_Sol_id;
				Usless_Sol_id = new int[Useless_Sol_Size];
				MPI::Comm.Recv( Usless_Sol_id, Useless_Sol_Size, MPI::INT, 0, search_param->getTag_Useless_Sol() );			
				
				for(int i=0; i < Useless_Sol_Size; i++)  {
	//  			cout << endl << "sol id will not be evaluated!: " << Usless_Sol_id[i]  << endl;
					for(int sss = 0; sss < n_sc_worker; sss++)
					  if(Usless_Sol_id[i] <Y_Pool.size() )//because sometimes it may try to remove a solution that may not arrived
						Y_Pool[Usless_Sol_id[i]][n_arcs + 1 + sss] = 1;			    
				}			  
				delete [] Usless_Sol_id;
		      }
		      
		      
		      //pick a solution from the pool (if any)		      
		      for(int i= Y_Pool.size() -1; i >=0; i--)//LIFO
				if(Y_Pool[i][n_arcs + 1 + s] == 0){
					Y_Pool[i][n_arcs + 1 + s] =1;//it is used for scenario s (sc_id)			  
					sol_id =Y_Pool[i][0];
					for(int a = 0; a < n_arcs; a++)
					  y_SOL[a] = Y_Pool[i][1+ a];				      
					solve_sub_problems =true;	
					break;
				  }
// 		      cout << endl << "Slave: " << id << "  ->" << Y_Pool.size() << " selc " << sol_id << endl;
		      /* for(int i= 0; i < Y_Pool.size() ; i++)//FIFO
			if(Y_Pool[i][n_arcs + 1 + s] == 0)
			{
			    Y_Pool[i][n_arcs + 1 + s] =1;//it is used for scenario s (sc_id)
			  
			    sol_id =Y_Pool[i][0];
			    for(int a = 0; a < n_arcs; a++)
			      y_SOL[a] = Y_Pool[i][1+ a];
				      
			    solve_sub_problems =true;	
			    break;
		      }*/
		      		 
			return sol_has_rec;
		}
///////////////////////////////////////////////////////////////////////////////////////		
		void fill_send_vector(int n_od, bool infPapaSubProblem, bool infSubProblem, double ***Send_cut, int sol_id, int sc_id,  double *g_Values,  double *G, double Sub_Problem_Obj, double tEndworker, int n_arcs, int n_sc, int s)
		{
		    Send_cut[s][sol_id][0] = sol_id;
			Send_cut[s][sol_id][1] = sc_id;				
			Send_cut[s][sol_id][2] = Sub_Problem_Obj;
			Send_cut[s][sol_id][3] = tEndworker;//(*v_aux)[1]; //ss obj val		
			//cut from Papadakos SP
			for(int ss=0; ss<n_sc+1; ss++)
				Send_cut[s][sol_id][4 + ss] = g_Values[ss];//Send_cut[s][sol_id][4] = g_Values[0];//(*v_aux)[0]; // assign g for the actual scenario we have solved 						   											    
			
			
			
			
			
			/*for(int a =0; a < n_arcs ; a++)
				Send_cut[s][sol_id][4 + n_sc +1 + a] = G[a];
			int t=0;
			for(int a = 0; a < n_arcs; a++){
				Send_cut[s][sol_id][4 + n_sc +1 + a] = G[t++];
				for (int k = 0; k < n_od; k++)
					Send_cut[s][sol_id][4 + n_sc + 1 + n_arcs + t] = G[t++];
					t++;
				}
				*/
				
				int t=0, tt=0;
				for(int a = 0; a < n_arcs; a++)  {
					Send_cut[s][sol_id][4 + n_sc + 1 + tt++] = G[t++];// = aux1 - data_S->getU(a) * val_reg[ttt++] ; //p[s]*
					for (int k = 0; k < n_od; k++)
						Send_cut[s][sol_id][4 + n_sc + 1 + tt++] = G[t++];// = -/* data_S->getD(k,s_id) * */ val_reg[ttt++] ;//here we add the dual values associated to the strong inequalities to the coiffiecient of the design y variables										      												      
				}
				
				
				
			//to the end of information regarding the generated cut which I am sending to the master, I add the required information to generate the valid cut for other scenarios too.
			//cut from regular SP
			/*Send_cut[s][sol_id][5+n_arcs] = g_Values[1];	
			for(int a =0; a < n_arcs ; a++)
				Send_cut[s][sol_id][6 + n_arcs + a] = G[n_arcs+a];
						*/		
			//indicator of solution feasibility or infeasibility for this scenario for the regular one
			if(infSubProblem)
				Send_cut[s][sol_id][4 + n_sc + 1 + n_arcs +  n_arcs*n_od] = 0;//infeasible was the solution for this scenario sp
			else
				Send_cut[s][sol_id][4 + n_sc + 1 + n_arcs +  n_arcs*n_od] = 1; //feasible it was
			//indicator of solution feasibility or infeasibility for this scenario for the Papadakos one
			if(infPapaSubProblem)
				Send_cut[s][sol_id][4 + n_sc + 1 + n_arcs +  n_arcs*n_od + 1] = 0;//infeasible was the solution for this scenario sp
			else
				Send_cut[s][sol_id][4 + n_sc + 1 + n_arcs +  n_arcs*n_od + 1] = 1; //feasible it was			
		}
		
////////////////////////////////////////////////////////////////////////////////////////////		
		void UB_Other_Sce(int s, int sol_id, double ***Send_cut,  double *y_SOL, double **X_Value, int sc_id, double Sub_Problem_Obj, int n_sc, Data_S* data_S)
		{
			
		    double  obj, Res_Cap, min_Res_Cap, arc_flow, passed_flow;
		    int source0, source, destination;
		    vector < int >  v_aux, v_aux1;
		    bool any_path;
		    double **aux_X_Value;
		    aux_X_Value = new double*[data_S->getN_arcs()];
		    for(int a =0; a <data_S->getN_arcs(); a++)
		      aux_X_Value[a] = new double[data_S->getN_od()];
		    
		    
		    
//  		    cout << endl << sc_id << " !!!--> " << Sub_Problem_Obj << endl;
		    for(int ss = 0; ss < n_sc; ss++)		// flow conservation constraints for global scenarios
		    {
		      obj = (1/data_S->getP(sc_id))*Sub_Problem_Obj;
		      if(ss != sc_id)
			{		
			  for(int a =0; a <data_S->getN_arcs(); a++)
			    for(int k =0; k <data_S->getN_od(); k++)
			      aux_X_Value[a][k] = X_Value[a][k];
				//for removing the extra flows on the network
				for(int k = 0; k < data_S->getN_od(); k++ ) 
				  if( data_S->getD(k,sc_id) > data_S->getD(k,ss) && data_S->getOd(k,0) != data_S->getOd(k,1))
					{					       					  																				      							     														  
							      v_aux.clear();
							      v_aux1.clear();
							      source = data_S->getOd(k,0);
							      destination = data_S->getOd(k,1);
							      v_aux.push_back(source);
							      any_path =true;//there is a path
							     //finding paths							  							      							
							     while(source != destination )
								  {		
									    source0 = v_aux[v_aux.size() -1];
									    for(int a =0; a < data_S->getN_arcs(); a++ )//it is better to remove flow from dummy arcs first(if they have any)										     
											  if ( (data_S->getArcs(a,0) == source)  && (aux_X_Value[a][k] >1e-2) )//a path that we can decrease flow on it
											      {	
													  source = data_S->getArcs(a,1);
													  v_aux.push_back(source);
													  v_aux1.push_back(a); 												    
													  break;												      
											      }										      
										      
									    if(source0 == source )//a path is not found or the path is not useable since we cannot augment any flow on it
									      any_path = false;	
									    if(!any_path)
									      break;
								  }//while(source != destination )
							      							  							  
							    //augmenting flows ..> increasing/reducing flows on the network to generate feasible one for the current sub-problem s							  
								if( any_path)
									for(int j =0; j < v_aux1.size(); j++)//reduce the additional flow on the pÃ¢th	
									{
									  obj += -data_S->getC(v_aux1[j]) * min(aux_X_Value[v_aux1[j]][k], data_S->getD(k,sc_id) - data_S->getD(k,ss));	//( d[k][s] - d[k][sc_id]);	
									  aux_X_Value[v_aux1[j]][k] -= min(aux_X_Value[v_aux1[j]][k], data_S->getD(k,sc_id) - data_S->getD(k,ss));
									}												          
					}
				//for agumenting the surplus flows of this scenario	
				for(int k = 0; k < data_S->getN_od(); k++ ) 
				  if(data_S->getD(k,ss) > data_S->getD(k,sc_id) && data_S->getOd(k,0) != data_S->getOd(k,1))
					  {								      
							      v_aux.clear();
							      v_aux1.clear();
							      source =data_S->getOd(k,0);
							      destination = data_S->getOd(k,1);
							      v_aux.push_back(source);
							      any_path =true;//there is a path
							      min_Res_Cap = data_S->getD(k,ss) - data_S->getD(k,sc_id);//to find the max flow that we can augment through the path
							     //finding paths							  							      							
							     while(source != destination )
								  {		
									          source0 = v_aux[v_aux.size() -1];									  
										  for(int a =0; a <data_S->getN_arcs(); a++)										    
											  if ( (data_S->getArcs(a,0) == source) && (aux_X_Value[a][k] >1e-2) && (aux_X_Value[a][k] < data_S->getD(k,ss)) && (data_S->getArcs(a,0) != data_S->getArcs(a,1)) )//we want to find a path that can carry extra flow
											  {
												  arc_flow =0;//count total flow on this arc
												  for(int kk = 0; kk < data_S->getN_od(); kk++)
												    arc_flow +=  aux_X_Value[a][kk];
												  if(data_S->getU(a)*y_SOL[a] - arc_flow > 1e-2)//the point is to find a useable path
												    {
													  source = data_S->getArcs(a,1);
													  v_aux.push_back(source);
													  v_aux1.push_back(a);
													  if(data_S->getU(a)*y_SOL[a] - arc_flow < min_Res_Cap)//min cap on the path
													    min_Res_Cap = data_S->getU(a)*y_SOL[a] - arc_flow;
													  break;
												    }	
											  }										    
									    if(source0 == source )//a path is not found or the path is not useable since we cannot augment any flow on it
									      any_path = false;
									    if(!any_path)
									      break;
									        
								  }//while(source != destination )							      							  							  
// 							    //augmenting flows ..> increasing/reducing flows on the network to generate feasible one for the current sub-problem s							  
								
								      if(any_path)
									for(int j =0; j < v_aux1.size(); j++)//pass as much as possible using the residual capacity
									{
									    obj += data_S->getC(v_aux1[j]) * min(min_Res_Cap, data_S->getD(k,ss) - data_S->getD(k,sc_id));
									    aux_X_Value[v_aux1[j]][k] += min(min_Res_Cap, data_S->getD(k,ss) - data_S->getD(k,sc_id));
									}
									  passed_flow=0;
									  for(int aa = 0; aa < data_S->getN_arcs(); aa++)
									    if( (data_S->getArcs(aa,0) == data_S->getOd(k,0)) && (aux_X_Value[aa][k] > 1e-1) )
									      passed_flow +=aux_X_Value[aa][k];
									    
									  obj += data_S->getC(data_S->getN_arcs()-1) * max(0.0, data_S->getD(k,ss) - passed_flow);//the amount that we cannot send through path will be directed twoard dummy arcs						      									    																					
					 }
					
//         				cout << endl << ss << " --> " << p[ss]*obj;
				
			}
			
			Send_cut[s][sol_id][5+data_S->getN_arcs() + n_sc + n_sc * data_S->getN_arcs() + ss] =data_S->getP(ss)*obj;
			
		    }
			
		}
		

		
		
		
		
		~sub_manager()
		{	
		  
		}

};
#endif
