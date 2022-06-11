#ifndef CLASS_OPT_H
#define CLASS_OPT_H

#include <cmath>
#include <ilcplex/ilocplex.h>
ILOSTLBEGIN


typedef IloArray<IloNumVarArray> IloNumVarArray2;
typedef IloArray<IloNumVarArray2> IloNumVarArray3;
typedef IloArray<IloNumVarArray3> IloNumVarArray4;

class OptimalcutModel
{ 
	private:

		IloEnv env, MW_env;
		IloModel model, MW_model;
		IloCplex cplex, MW_cplex;
		IloRangeArray con, con1, MW_con, MW_con1, MW_Eq_Con; //this is for the main constraints (flow conservation and capacity) constraints				
		IloNumVarArray Dual_Capacity_Con, MW_Dual_Capacity_Con;		
		IloNumVarArray2  Dual_Flow_Con, Dual_Strong_Ineq, MW_Dual_Flow_Con, MW_Dual_Strong_Ineq;
			
		IloNumVarArray var;
		IloNumArray val;
		IloNumArray X;		 		
		IloObjective objFun, MW_objFun;
		
		//**************************		
	public:
		//create the sub-problems
		OptimalcutModel(Data_S *data_S, Search_Param *search_param) 
		{
			Build_SP(data_S, search_param)	;//build the regular sub problems
			Build_MW_SP(data_S, search_param);//build the MW sub-problem			
		}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		void Build_SP(Data_S *data_S, Search_Param *search_param)
		{
			model 		=  IloModel(env);
			cplex 		=  IloCplex(model);
			
			con 		=  IloRangeArray(env);
			con1 		=  IloRangeArray(env);
			
			objFun		= IloMaximize(env);						
 			
			//________________Variables definition_____________
			Dual_Flow_Con = IloNumVarArray2(env, data_S->getN_nodes());
			//			 IloNumVarArray2 Dual_Flow_Con(env, n_nodes);//if I define it in this way I need the "this" because I am defining a new one but similar
			for(int i = 0; i < data_S->getN_nodes(); i++)
			  Dual_Flow_Con[i] = IloNumVarArray(env, data_S->getN_od(), -IloInfinity, IloInfinity); // for the flow conversation constraints 			
			//			  this->Dual_Flow_Con = Dual_Flow_Con;
			
 			for(int i = 0; i < data_S->getN_nodes(); i++)
 			    (model).add(Dual_Flow_Con[i]);
			    
			Dual_Capacity_Con = IloNumVarArray(env, data_S->getN_arcs(), 0, IloInfinity);	//for the capacity constraint 
			//			IloNumVarArray Dual_Capacity_Con(env, n_arcs, 0, IloInfinity);	//for the capacity constraint 
			//			this->Dual_Capacity_Con = Dual_Capacity_Con;
			(model).add(Dual_Capacity_Con);
			
						
			Dual_Strong_Ineq = IloNumVarArray2(env, data_S->getN_arcs());
			//			IloNumVarArray2 Dual_Strong_Ineq(env, n_arcs);
			for(int a = 0; a < data_S->getN_arcs(); a++)
			  Dual_Strong_Ineq[a] = IloNumVarArray(env, data_S->getN_od(), 0, IloInfinity); // for the strong ineq 
			//			this->Dual_Strong_Ineq = Dual_Strong_Ineq;
			
 			for(int a = 0; a < data_S->getN_arcs(); a++)
 			    (model).add(Dual_Strong_Ineq[a]);
// 			  			
			//________________Obejective function______________
// 			int Dest_K;
			IloExpr expr(env);
			for(int k = 0; k < data_S->getN_od(); k++) //this first part (associated to the flow constraints)
			  expr += fabs(data_S->getD(k,0)) * (Dual_Flow_Con[data_S->getOd(k,1)][k]);

			
			for (int a = 0; a < data_S->getN_arcs() ; a++)
			  if (data_S->getArcs(a,1) != data_S->getArcs(a,0))//if it is not looping arc
			  {
			    expr += -data_S->getU(a) * Dual_Capacity_Con[a];//the part associated with capacity constraints
			    for (int k = 0; k < data_S->getN_od(); k++)
			      expr += -data_S->getU(a) * Dual_Strong_Ineq[a][k];
			  }
			
			objFun.setExpr(expr);
			
			(model).add(objFun);
			expr.end();						
			//_____________Constraints_______________________
			
// 			int i , j; //this is to have the head and tail of arc a
			 for ( int a = 0; a < data_S->getN_arcs() ; a++)
			   if (data_S->getArcs(a,1) != data_S->getArcs(a,0))//if it is not looping arc
			     for ( int k = 0; k < data_S->getN_od(); k++)			   			   
			      {			     
      // 			     i = arcs[a][0]; //tail of arc a (i)
      // 			     j = arcs[a][1]; // head (end) of arc "a" (j)
// 				  	cout << arcs[a][0] << " - " <<	arcs[a][1] << endl;		
				    IloExpr expr(env);
				    expr +=  Dual_Flow_Con[data_S->getArcs(a,1)][k];
				    expr -=  Dual_Flow_Con[data_S->getArcs(a,0)][k];
				    expr -=  Dual_Capacity_Con[a];
				    expr -=  Dual_Strong_Ineq[a][k];
				    (con).add(IloRange(env, -IloInfinity, expr , data_S->getC(a)));
				    expr.end();				
			      }
			
			//the dual variable associated to the origin is equal to zero
			for(int k = 0; k < data_S->getN_od(); k++)
			  (con1).add(IloRange(env, 0, Dual_Flow_Con[data_S->getOd(k,0)][k] , 0));
			
			
			//adding the strong inequalities to the sub-problem formulation
			if(!search_param->getStrong_Ineq())//if it is false we remove the dual variables associate to strong inequalities (fix them at zero)
			{			  
			   for ( int a = 0; a < data_S->getN_arcs() ; a++)
				for ( int k = 0; k < data_S->getN_od(); k++)
				  (con1).add(IloRange(env, 0,  Dual_Strong_Ineq[a][k] ,0));			 
			}

		        (model).add((con));
   			(model).add((con1));						
			
						
 			X = IloNumArray(env,(con).getSize());
			
			//(*cplex).exportModel("./Outputs/Cuts/OPTCUT_Lp_1.lp");
			(cplex).setParam(IloCplex::TiLim, 3600);
 			(cplex).setParam(IloCplex::SimDisplay, 0);
			(cplex).setParam(IloCplex::Threads, 1);
 			(cplex).setParam(IloCplex::AdvInd, 0);	
			(cplex).setParam(IloCplex::EpGap, 0.00);
		}
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		void Build_MW_SP(Data_S *data_S, Search_Param *search_param)
		{
		  
			MW_model 	=  IloModel(MW_env);
			MW_cplex 	=  IloCplex(MW_model);
			
			MW_con 		=  IloRangeArray(MW_env);
			MW_con1		=  IloRangeArray(MW_env);
			MW_Eq_Con 	=  IloRangeArray(MW_env);
			
			MW_objFun	=  IloMaximize(MW_env);
			
			
 			
			//________________Variables definition_____________
			MW_Dual_Flow_Con = IloNumVarArray2(MW_env, data_S->getN_nodes());
			//			 IloNumVarArray2 Dual_Flow_Con(env, n_nodes);//if I define it in this way I need the "this" because I am defining a new one but similar
			for(int i = 0; i < data_S->getN_nodes(); i++)
			  MW_Dual_Flow_Con[i] = IloNumVarArray(MW_env, data_S->getN_od(), -IloInfinity, IloInfinity); // for the flow conversation constraints 			
			//			  this->Dual_Flow_Con = Dual_Flow_Con;
			
 			for(int i = 0; i < data_S->getN_nodes(); i++)
 			    (MW_model).add(MW_Dual_Flow_Con[i]);
			    
			MW_Dual_Capacity_Con = IloNumVarArray(MW_env, data_S->getN_arcs(), 0, IloInfinity);	//for the capacity constraint 
			//			IloNumVarArray Dual_Capacity_Con(env, n_arcs, 0, IloInfinity);	//for the capacity constraint 
			//			this->Dual_Capacity_Con = Dual_Capacity_Con;
			(MW_model).add(MW_Dual_Capacity_Con);
			
						
			MW_Dual_Strong_Ineq = IloNumVarArray2(MW_env, data_S->getN_arcs());
			//			IloNumVarArray2 Dual_Strong_Ineq(env, n_arcs);
			for(int a = 0; a < data_S->getN_arcs(); a++)
			  MW_Dual_Strong_Ineq[a] = IloNumVarArray(MW_env, data_S->getN_od(), 0, IloInfinity); // for the strong ineq 
			//			this->Dual_Strong_Ineq = Dual_Strong_Ineq;
			
 			for(int a = 0; a < data_S->getN_arcs(); a++)
 			    (MW_model).add(MW_Dual_Strong_Ineq[a]);
// 			  			
		
  
			//________________Obejective function______________
// 			int Dest_K;
			IloExpr expr(MW_env);
			for(int k = 0; k < data_S->getN_od(); k++) //this first part (associated to the flow constraints)
			  expr += fabs(data_S->getD(k,0)) * (MW_Dual_Flow_Con[data_S->getOd(k,1)][k]);

			
			for (int a = 0; a < data_S->getN_arcs() ; a++)
			  if (data_S->getArcs(a,1) != data_S->getArcs(a,0))//if it is not looping arc
			  {
			    expr += -data_S->getU(a) * MW_Dual_Capacity_Con[a];//the part associated with capacity constraints
			    for (int k = 0; k < data_S->getN_od(); k++)
			      expr += -data_S->getU(a) * MW_Dual_Strong_Ineq[a][k];
			  }
			
			MW_objFun.setExpr(expr);
			
			(MW_model).add(MW_objFun);
			expr.end();						
			//_____________Constraints_______________________
			
// 			int i , j; //this is to have the head and tail of arc a
			 for ( int a = 0; a < data_S->getN_arcs() ; a++)
			   if (data_S->getArcs(a,1) != data_S->getArcs(a,0))//if it is not looping arc
			     for ( int k = 0; k < data_S->getN_od(); k++)			   			   
			      {			     
      // 			     i = arcs[a][0]; //tail of arc a (i)
      // 			     j = arcs[a][1]; // head (end) of arc "a" (j)
// 				  	cout << arcs[a][0] << " - " <<	arcs[a][1] << endl;		
				    IloExpr expr(MW_env);
				    expr +=  MW_Dual_Flow_Con[data_S->getArcs(a,1)][k];
				    expr -=  MW_Dual_Flow_Con[data_S->getArcs(a,0)][k];
				    expr -=  MW_Dual_Capacity_Con[a];
				    expr -=  MW_Dual_Strong_Ineq[a][k];
				    (MW_con).add(IloRange(MW_env, -IloInfinity, expr , data_S->getC(a)));
				    expr.end();				
			      }
			
			//the dual variable associated to the origin is equal to zero
			for(int k = 0; k < data_S->getN_od(); k++)
			  (MW_con1).add(IloRange(MW_env, 0, MW_Dual_Flow_Con[data_S->getOd(k,0)][k] , 0));						
			
			//adding the strong inequalities to the sub-problem formulation
			if(!search_param->getStrong_Ineq())//if it is false we remove the dual variables associate to strong inequalities (fix them at zero)
			{			  
			   for ( int a = 0; a < data_S->getN_arcs() ; a++)
				for ( int k = 0; k < data_S->getN_od(); k++)
				  (MW_con1).add(IloRange(MW_env, 0,  MW_Dual_Strong_Ineq[a][k] ,0));			 
			}

		        (MW_model).add((MW_con));
   			(MW_model).add((MW_con1));						
			
						
 			
			//for getting dual values
			var = IloNumVarArray( MW_env, data_S->getN_nodes()*data_S->getN_od() +  data_S->getN_arcs() + data_S->getN_arcs()*data_S->getN_od());
			val = IloNumArray(MW_env, data_S->getN_nodes()*data_S->getN_od() +  data_S->getN_arcs() + data_S->getN_arcs()*data_S->getN_od());
			IloInt t=0;
			
			for(int i = 0; i < data_S->getN_nodes(); i++)
			  for(int k = 0; k < data_S->getN_od(); k++)		//
			     var[t++] =MW_Dual_Flow_Con[i][k];
		        for(int a = 0; a < data_S->getN_arcs(); a++){		//
			   var[t++] =MW_Dual_Capacity_Con[a];
			  for(int k = 0; k < data_S->getN_od(); k++)	//
			    var[t++] =MW_Dual_Strong_Ineq[a][k];			    
			}

  
			//(*cplex).exportModel("./Outputs/Cuts/OPTCUT_Lp_1.lp");
			(MW_cplex).setParam(IloCplex::TiLim, 3600);
  			(MW_cplex).setParam(IloCplex::SimDisplay, 0);
			(MW_cplex).setParam(IloCplex::Threads, 1);
 			(MW_cplex).setParam(IloCplex::AdvInd, 0);	
			(cplex).setParam(IloCplex::EpGap, 0.00);
		}
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		double solve(Search_Param *search_param, Data_S* data_S, double **X_Value, double* y_SOL, double** core_y_SOL, int s, int actual_sc, int n_sc, double* G, double* g_Values) 
		{
			  float ObjectiveOfScenario;
			//****************************** SP ****************************************
			//we first solve the regular Benders sub-problems
			// updating constraints			
			for (int a= data_S->getN_arcs() - data_S->getN_od(); a<data_S->getN_arcs(); a++)
			  y_SOL[a] = 1;
// 			int Dest_K;
			IloExpr expr(env);
			for(int k = 0; k < data_S->getN_od(); k++)
			  expr += data_S->getD(k,actual_sc) * Dual_Flow_Con[data_S->getOd(k,1)][k];//doing it in the expr way is always faster than setcoef			
			for (int a = 0; a < data_S->getN_arcs() ; a++)
			  if (data_S->getArcs(a,1) != data_S->getArcs(a,0))//if it is not looping arc
			    {			      
			      expr -=  data_S->getU(a) * y_SOL[a] * Dual_Capacity_Con[a];
			      if(search_param->getStrong_Ineq())//if we have strong ineq in the sub-problems
				for (int k = 0; k < data_S->getN_od(); k++)
				  expr -= y_SOL[a] * min(data_S->getD(k,actual_sc),data_S->getU(a)) * Dual_Strong_Ineq[a][k];		      
			    }
			objFun.setExpr(expr);			
			 expr.end();//is the obj of regular sub_Problem
				
	
			
			if( !(cplex).solve() ) 
			{				  
				  if ((cplex).getStatus() == IloAlgorithm::Infeasible)
				    cout << endl << "regular sub_Problem infeasibile"<<endl;//(cplex).getStatus() << "        ->" << IloCplex.Status;
	//  			  for (int a=0; a<n_arcs;a++)
	//  			    if (y_SOL[a] !=0)
	//  			      cout << endl << a << "  --> " << y_SOL[a];
				  env.error() << "Failed to optimize regular sub-problem for scenario " << actual_sc << endl;

				  ObjectiveOfScenario= 9e10;
			}
			else//solve M-W
			  if( (cplex).getStatus() == IloAlgorithm::Unbounded )
			    cout  << "Regular sub_Problem unbounded"<< endl;
			  else
			      if ((cplex).getStatus() == IloAlgorithm::Optimal)
			      {			  				  
				      ObjectiveOfScenario = (cplex).getObjValue();				    
//  					cout << endl <<  ObjectiveOfScenario - (cplex).getObjValue() << endl;
				    (cplex).getDuals(X, con);
				      int t =0;
					for ( int a = 0; a < data_S->getN_arcs() ; a++)
					  if (data_S->getArcs(a,1) != data_S->getArcs(a,0))
					  for ( int k = 0; k < data_S->getN_od(); k++)				      
					    X_Value[a][k] = X[t++];//(cplex).getDual(con[t++]);
					
				      solve_MW(ObjectiveOfScenario, search_param, data_S,  y_SOL, core_y_SOL, s, actual_sc, n_sc, G, g_Values); 
			    }//end else 

			return data_S->getP(actual_sc) * ObjectiveOfScenario ;

		}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		void solve_MW(long double ObjectiveOfScenario, Search_Param *search_param, Data_S* data_S, double* y_SOL, double** core_y_SOL, int s, int actual_sc, int n_sc, double* G, double* g_Values)
		{
		  //******************************* Independent M-W based on the Papadakos*****************
				  //to solve the SP for generating the Pareto optimal cut
			  
// 				  for (int a= data_S->getN_arcs() - data_S->getN_od(); a <data_S->getN_arcs(); a++)
// 				    core_y_SOL[s][a] = 1;
				  	
				  //objFun
				  IloExpr expr(MW_env);
				 for(int k = 0; k < data_S->getN_od(); k++)
				    expr += data_S->getD(k,actual_sc) * MW_Dual_Flow_Con[data_S->getOd(k,1)][k];//doing it in the expr way is always faster than setcoef				  
				  for (int a = 0; a < data_S->getN_arcs() ; a++)
				    if (data_S->getArcs(a,1) != data_S->getArcs(a,0))//if it is not looping arc
				      {					
					expr -=  core_y_SOL[s][a] * data_S->getU(a) *  MW_Dual_Capacity_Con[a];
					if(search_param->getStrong_Ineq())
					  for (int k = 0; k < data_S->getN_od(); k++)
					    expr -= core_y_SOL[s][a] * min(data_S->getD(k,actual_sc), data_S->getU(a)) * MW_Dual_Strong_Ineq[a][k];				  					
				      }
				  MW_objFun.setExpr(expr);	
				  
				  expr.end();
					  	
				  
				  //MW constraint				  
				  IloExpr expr1(MW_env);
				  for(int k = 0; k < data_S->getN_od(); k++)
				    expr1 += data_S->getD(k,actual_sc) * MW_Dual_Flow_Con[data_S->getOd(k,1)][k];//doing it in the expr way is always faster than setcoef			
				  for (int a = 0; a < data_S->getN_arcs() ; a++)
				    if (data_S->getArcs(a,1) != data_S->getArcs(a,0))//if it is not looping arc
				      {			      
					expr1 -=  data_S->getU(a) * y_SOL[a] * MW_Dual_Capacity_Con[a];
					if(search_param->getStrong_Ineq())//if we have strong ineq in the sub-problems
					  for (int k = 0; k < data_S->getN_od(); k++)
					    expr1 -= y_SOL[a] * min(data_S->getD(k,actual_sc),data_S->getU(a)) * MW_Dual_Strong_Ineq[a][k];		      
				      }
					    				   				  
				  MW_Eq_Con.endElements();				  
//   				  MW_Eq_Con.clear();
				  
				  (MW_Eq_Con).add(IloRange(MW_env, (cplex).getObjValue()/*ObjectiveOfScenario*/ - 1e-2, expr1 , (cplex).getObjValue()/*ObjectiveOfScenario*/ + 1e-2) );
				  (MW_model).add(MW_Eq_Con);
				  expr1.end();
				  
				  
				  if( !(MW_cplex).solve() ) 
				  {				    				   				    
				    if ((MW_cplex).getStatus() == IloAlgorithm::Infeasible)
				      cout << endl << "M-W sub_Problem infeasibile"<< endl;//(cplex).getStatus() << "        ->" << IloCplex.Status;				  				    				    
				    env.error() << "Failed to optimize M-W sub-problem for scenario " << actual_sc << endl;
				    MPI::Comm.Abort(911);
				  }
				  				  				  				  
				 if ((MW_cplex).getStatus() == IloAlgorithm::Unbounded )			    
				      cout << endl << "M-W sub_Problem unbounded"<< endl;//(cplex).getStatus() << "        ->" << IloCplex.Status;
				  if ((MW_cplex).getStatus() == IloAlgorithm::Optimal) 
				  {		
				    (MW_cplex).getValues(val, var);
				    //____________________________________ g values _______________________
					    //for a cut we need both g anf G; therefore based on the current dual values we calculate these values for each scenario sub-problem					    
					    double aux;
					    int t = 0, tt = 0, ttt;	//a simple counter
					  for(int ss = 0; ss < n_sc; ss++)//I wrote a loop on the scenarios because I want to copy every cut I get for a given scenario to other scenarios
					  {
						    aux = 0;
						    ttt=0;
						    for(int i = 0; i < data_S->getN_nodes(); i++)
						      for(int k = 0; k < data_S->getN_od(); k++)
							if (i == data_S->getOd(k,1))
							  aux += fabs(data_S->getD(k,ss)) * val[ttt++]/*(cplex).getValue(Dual_Flow_Con[i][k])*/;						    
							else
							  ttt++;
						    g_Values[t++] = aux;	//this gives the g for every scenario based on the solution of this scenario												
					//____________________________________ G[a] values _______________________	
					      
						for(int a = 0; a < data_S->getN_arcs(); a++)
						{
						      aux = 0.0;
						      aux = -data_S->getU(a) * val[ttt++] /*(cplex).getValue(Dual_Capacity_Con[a])*/; //p[s]*
						      for (int k = 0; k < data_S->getN_od(); k++)
							      aux += -data_S->getD(k,ss) *val[ttt++] /*(cplex).getValue(Dual_Strong_Ineq[a][k])*/;//here we add the dual values associated to the strong inequalities to the coiffiecient of the design y variables												      
						      G[tt++] = aux; 						      
						}					      
					  }//end of for on ss   					      					    
				  }
			   
			   
			MW_model.remove(MW_Eq_Con);
// 			MW_model.remove(MW_Con);
		}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

		~OptimalcutModel()
		{
			(cplex).end();
			(model).end();			
			(con).end();
			env.end();			
		}

};
#endif
