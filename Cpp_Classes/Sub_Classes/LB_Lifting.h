#ifndef CLASS_LB_INEQ_H
#define CLASS_LB_INEQ_H





ILOSTLBEGIN
#define Comm COMM_WORLD



class class_LB_Lifting
{ 
  private:    
		double *min_d; //for the vector of minimum demand of each commodity over all scenarios
		
		IloEnv env;		//environment 
		IloModel model;
		IloCplex cplex;		
 		IloNumVarArray2 x;
 		IloNumVarArray  y;		
		IloNumArray ySol;
		IloRangeArray *con, *con1, *fixCon;
			
		IloNumVarArray var;
		IloNumArray val;
		IloNumArray X;		 		
		IloObjective objFun;
		
		IloNumArray y_SOL;
		IloNumArray2 x_SOL; 
		IloNumArray2 aux_dual_SOL;
		IloNumArray dual_SOL;
		double routing_cost;
		double max_LBF;
		double *scenario_ave;
		int *cluster_id;//scenario s is in which cluster c	
		vector < vector <double > > Cluster_Core;
			  
  public:        		
		class_LB_Lifting(int n_sc, MasterSolChoice *Sol_Choice, Data_S *data_S, Search_Param	*search_param, int sID)		
		{	
			model 		= IloModel(env);
			cplex 		= IloCplex(model);
			con   		= new IloRangeArray(env);
			con1   		= new IloRangeArray(env);
			objFun  		= IloMinimize(env);	
			
			y_SOL = IloNumArray(env, data_S->getN_arcs());
			dual_SOL = IloNumArray(env, data_S->getN_od());
			aux_dual_SOL = IloNumArray2(env, data_S->getN_od());
			for(int k = 0; k < data_S->getN_od(); k++)
				aux_dual_SOL[k] = IloNumArray(env, data_S->getN_nodes());
			x_SOL = IloNumArray2(env, data_S->getN_arcs());
			for(IloInt a = 0; a < data_S->getN_arcs(); a++)
			     x_SOL[a] = IloNumArray(env, data_S->getN_od());		      			
			
			n_sc = n_sc - search_param->getMaster_sc(); // number of scenarios are total scenarios minus global scenarios (those who are not projected)	 			
			min_demand(data_S, n_sc, sID);			
			//***the master only includes y and theta variables****
			y = IloNumVarArray(env, data_S->getN_arcs(), 0, 1);	//ILOINT defining an variable for each arc (on dimensional variable)	
			(model).add(y); // we dont need to add variables explicitly to the model as they will implicitly be added to the model through the constraints
			x = IloNumVarArray2(env, data_S->getN_arcs());//second-stage (x) variable (three dimensional cuz each arc (i-j) is represented with a single name "a"
			for(int a = 0; a < data_S->getN_arcs(); a++){
				x[a] = IloNumVarArray(env, data_S->getN_od(), 0, IloInfinity); 		
				// for(int k = 0; k < data_S->getN_od(); k++)				
					// x[a][k] = IloNumVarArray(env, n_sc + search_param->getMaster_sc(), 0, IloInfinity); 									
			}
// 			this->x=x;
			for(int a = 0; a < data_S->getN_arcs(); a++)
			    (model).add(x[a]);
			// ***** OBJECTIVE FUNCTION****
			IloExpr expr(env);
			for(int a = 0; a < data_S->getN_arcs() /* - data_S->getN_od() */; a++){
			  expr += data_S->getF(a)*y[a];			  
			  for(int k = 0; k < data_S->getN_od(); k++)			//the second stage objective for the global scenarios 
			      // for(int s = 0; s < n_sc + search_param->getMaster_sc(); s++){			     
					expr +=  data_S->getC(a) * x[a][k];//(1/master_sc) * c[a] * x[a][k][s];//
			      
			}
			
			objFun.setExpr(expr);
			(model).add(objFun);
			expr.end();
			
			
			// CONSTRAINT 1				
			for(int k = 0; k < data_S->getN_od(); k++)				  
					for(int i = 0; i < data_S->getN_nodes(); i++){ 
						  IloExpr expr(env);//expr.clear();
						  IloExpr expr1(env);
						  for(int a =0; a < data_S->getN_arcs()/* - data_S->getN_od() */; a++){
							  if(data_S->getArcs(a,0) == i) 
								  expr += x[a][k];								  
							  if(data_S->getArcs(a,1) == i) 
								  expr1 -=  x[a][k];
						  }
						  if (i == data_S->getOd(k,0)){//if origin
							// cout << "origin " << k << " demand: " << maxD << endl;
						    (*con1).add(IloRange(env, min_d[k], expr+expr1, min_d[k] ));//for the origin only outgoing part is neccessary						
						  }
						 else
						   if (i == data_S->getOd(k,1)){//if destination
							(*con1).add(IloRange(env, -min_d[k], expr+expr1, -min_d[k]));						
						   }
						   else
							(*con1).add(IloRange(env , 0, expr + expr1, 0));						     
						  expr.end();
						  expr1.end();
					}			
			// CONSTRAINT 2				
			for(int a = 0; a < data_S->getN_arcs()- data_S->getN_od(); a++){
					IloExpr expr(env);					
					for(int k = 0; k < data_S->getN_od(); k++)
						expr += x[a][k];					
					expr -= data_S->getU(a)*y[a];
					(*con).add(IloRange(env, -IloInfinity, expr, 0));
					expr.end();
			}
			for(int a = 0; a < data_S->getN_arcs(); a++)
				for(int k = 0; k < data_S->getN_od(); k++)
					(*con).add(IloRange(env, -IloInfinity, x[a][k] - min_d[k]*y[a], 0));
							
			
			(model).add((*con1));
			(model).add((*con));
			
			
			(cplex).setParam(IloCplex::TiLim, 3600);
			(cplex).setParam(IloCplex::ClockType, 1);//to set it equal to cpu time since it is seq
			(cplex).setParam(IloCplex::SimDisplay, 0);			 					
			(cplex).setParam(IloCplex::EpGap, 0.03);
		    (cplex).setParam(IloCplex::Param::Threads, 1); 
			
			
			
			solveModel(data_S, n_sc);
		}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		void min_demand( Data_S * data_S, int n_sc, int sID)
		{
		      min_d = new double[ data_S->getN_od()];
		      double minDemand;
		      for(int k =0; k< data_S->getN_od(); k++){
				minDemand = -1e10;
				for(int s=0; s < n_sc; s++)
				  if(data_S->getD(k,s) >= minDemand)
					minDemand = data_S->getD(k,s);
				min_d[k] = data_S->getD(k,sID);//minDemand;
		      }		  
		}
////////////////////////////////////////////////////////////////
	      void solveModel( Data_S * data_S, int n_sc)
	      {
			(cplex).solve();
			cout << endl << "the cost of artificial scenario: " << cplex.getObjValue() <<endl;
			
			//get the dual variables Pi_D(k) associated to the flow conservation constraints
			int t=0;
			for(int k = 0; k < data_S->getN_od(); k++)
				for(int i = 0; i < data_S->getN_nodes(); i++)
					aux_dual_SOL[k][i] = (cplex).getDual( (*con1)[t++]);
	// 		env.out() << "dual values: "<< dual_SOL;
			for(int k = 0; k < data_S->getN_od(); k++){
				dual_SOL[k]=0;
				for(int i = 0; i < data_S->getN_nodes(); i++)
					if (i == data_S->getOd(k,0))
						dual_SOL[k] += aux_dual_SOL[k][i];
					else if (i == data_S->getOd(k,1))
							dual_SOL[k] -= aux_dual_SOL[k][i];
			}
			
			//getting the routing cost
			routing_cost=0;
			for(IloInt a = 0; a < data_S->getN_arcs(); a++)
			   for(int k =0; k< data_S->getN_od(); k++)
				 routing_cost += data_S->getC(a) * cplex.getValue(x[a][k]);
			   
			//getting the y solutions
			for(IloInt a = 0; a < data_S->getN_arcs(); a++)
					y_SOL[a] = cplex.getValue(y[a]);
			
			//printing the approximated cost for each scenario and finding the max among scenarios
			   max_LBF=0;//maximum approximated objective
			   for(int s=0; s< n_sc; s++){
				 double costs=0;
				 costs += routing_cost;
				  for(int k=0; k< data_S->getN_od(); k++)
					 costs += (data_S->getD(k,s) - getMinD(k)) * getDualVal(k);
				  for(int a=0; a < data_S->getN_arcs() - data_S->getN_od(); a++)
					 costs += data_S->getF(a) * getY(a);
				  if(max_LBF <= costs)
					max_LBF= costs;
			     // cout << endl << "approximated cost of scemario " << s << " is  " << costs ;
			   }
	      }
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	     double getY(int a) {return y_SOL[a];}
	     double getDualVal(int k) {return dual_SOL[k];}
	     double getRotuingCost() {return routing_cost;}
	     double getMinD(int k) {return min_d[k];}
	     double getMax_LBF() {return max_LBF;}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

		~class_LB_Lifting()
		{
		}


};//WND class



#endif
