#ifndef CLASS_AUX_OPT_H
#define CLASS_AUX_OPT_H

#include <cmath>
#include <ilcplex/ilocplex.h>
ILOSTLBEGIN


typedef IloArray<IloNumVarArray> IloNumVarArray2;
typedef IloArray<IloNumVarArray2> IloNumVarArray3;
typedef IloArray<IloNumVarArray3> IloNumVarArray4;

class AuxCutModel
{ 
	private:

		IloEnv env;
		IloModel model;
		IloCplex cplex;
		IloRangeArray con, con1, MW_Eq_Con; //this is for the main constraints (flow conservation and capacity) constraints				
		IloNumVarArray Dual_Capacity_Con;		
		IloNumVarArray2  Dual_Flow_Con, Dual_Strong_Ineq;
			
		IloNumVarArray var, auxVarColumn;
		IloNumArray conLB, conUB, varLB, varUB, colLB, colUB;	 		
		IloObjective objFun;
		double **X_Value;
		int numColums;
		int iteration;
		//**************************		
	public:
		//create the sub-problems
		AuxCutModel(Data_S *data_S, Search_Param *search_param) 
		{
			model 		=  IloModel(env);
			cplex 		=  IloCplex(model);
			con 		=  IloRangeArray(env);
			con1 		=  IloRangeArray(env);
			MW_Eq_Con	=  IloRangeArray(env);
			
			objFun		= IloMaximize(env);
			
			//________________Variables definition_____________
			Dual_Flow_Con = IloNumVarArray2(env, data_S->getN_nodes());
			//			 IloNumVarArray2 Dual_Flow_Con(env, n_nodes);//if I define it in this way I need the "this" because I am defining a new one but similar
			for(int i = 0; i < data_S->getN_nodes(); i++)
			  Dual_Flow_Con[i] = IloNumVarArray(env, data_S->getN_od(), -IloInfinity, IloInfinity); // for the flow conversation constraints 			
			//			  this->Dual_Flow_Con = Dual_Flow_Con;
			 			  			
			Dual_Capacity_Con = IloNumVarArray(env, data_S->getN_arcs(), 0, IloInfinity);	//for the capacity constraint 
			//			IloNumVarArray Dual_Capacity_Con(env, n_arcs, 0, IloInfinity);	//for the capacity constraint 
			//			this->Dual_Capacity_Con = Dual_Capacity_Con;
			// (model).add(Dual_Capacity_Con);
						
			Dual_Strong_Ineq = IloNumVarArray2(env, data_S->getN_arcs());
			//			IloNumVarArray2 Dual_Strong_Ineq(env, n_arcs);
			for(int a = 0; a < data_S->getN_arcs(); a++)
			  Dual_Strong_Ineq[a] = IloNumVarArray(env, data_S->getN_od(), 0, IloInfinity); // for the strong ineq 
			//			this->Dual_Strong_Ineq = Dual_Strong_Ineq;

			
			//for getting extreme ray of deflected sub-problem and regular one respectively 
			auxVarColumn = IloNumVarArray(env);
			var = IloNumVarArray( env, data_S->getN_nodes()*data_S->getN_od() +  data_S->getN_arcs() + data_S->getN_arcs()*data_S->getN_od());			
			IloInt t=0;			
			for(int i = 0; i < data_S->getN_nodes(); i++)
			  for(int k = 0; k < data_S->getN_od(); k++)		//
			     var[t++] =Dual_Flow_Con[i][k];
		    for(int a = 0; a < data_S->getN_arcs(); a++){		//
			   var[t++] =Dual_Capacity_Con[a];
			  for(int k = 0; k < data_S->getN_od(); k++)	//
			    var[t++] =Dual_Strong_Ineq[a][k];			    
			}	
			model.add(var);
			
			
			varLB = IloNumArray(env, data_S->getN_nodes()*data_S->getN_od() +  data_S->getN_arcs() + data_S->getN_arcs()*data_S->getN_od());//for the daul values of the regular SP
			varUB = IloNumArray(env, data_S->getN_nodes()*data_S->getN_od() +  data_S->getN_arcs() + data_S->getN_arcs()*data_S->getN_od());//for the daul values of the regular SP
			colLB = IloNumArray(env);
			colUB = IloNumArray(env);
			
			 
			//________________Obejective function______________
// 			int Dest_K;
			IloExpr expr(env);
			for(int k = 0; k < data_S->getN_od(); k++) //this first part (associated to the flow constraints)
			  expr += fabs(data_S->getD(k,0)) * (Dual_Flow_Con[data_S->getOd(k,1)][k]);
			
			for (int a = 0; a < data_S->getN_arcs() ; a++)
			  if (data_S->getArcs(a,1) != data_S->getArcs(a,0)){//if it is not looping arc
			    expr -= data_S->getU(a) * Dual_Capacity_Con[a];//the part associated with capacity constraints
			    for (int k = 0; k < data_S->getN_od(); k++)
			      expr -= data_S->getU(a) * Dual_Strong_Ineq[a][k];
			  }
			
			objFun.setExpr(expr);			
			(model).add(objFun);
			expr.end();						
			//_____________Constraints_______________________			
// 			int i , j; //this is to have the head and tail of arc a
			 for ( int a = 0; a < data_S->getN_arcs() ; a++)
			   if (data_S->getArcs(a,1) != data_S->getArcs(a,0))//if it is not looping arc
			     for ( int k = 0; k < data_S->getN_od(); k++){			     		
				    IloExpr expr(env);
				    expr +=  Dual_Flow_Con[data_S->getArcs(a,1)][k];
				    expr -=  Dual_Flow_Con[data_S->getArcs(a,0)][k];
				    expr -=  Dual_Capacity_Con[a];
				    expr -=  Dual_Strong_Ineq[a][k];
				    (con).add(IloRange(env, -IloInfinity, expr, data_S->getC(a)));					
				    expr.end();				
			      }
			
			
			
			//the dual variable associated to the origin nodes is equal to zero
			for(int k = 0; k < data_S->getN_od(); k++)
			  (con1).add(IloRange(env, 0, Dual_Flow_Con[data_S->getOd(k,0)][k] , 0));						
			//dual variable associated to capacity constraint of the dummies is equal to zero
			/* t = data_S->getN_arcs() - data_S->getN_od();
			for(int k = 0; k < data_S->getN_od(); k++)
				(con1).add(IloRange(env, 0, Dual_Capacity_Con[t + k], 0)); */ 
			//adding the strong inequalities to the sub-problem formulation
			if(!search_param->getStrong_Ineq())//if it is false we remove the dual variables associate to strong inequalities (fix them at zero)
			{			  
			   for ( int a = 0; a < data_S->getN_arcs() ; a++)
				for ( int k = 0; k < data_S->getN_od(); k++)
				  (con1).add(IloRange(env, 0,  Dual_Strong_Ineq[a][k] ,0));			 
			}

		    (model).add((con));
   			(model).add((con1));						
			
 						
			conLB = IloNumArray(env, (con).getSize());
			conUB = IloNumArray(env, (con).getSize());
			t=0;
			for ( int a = 0; a < data_S->getN_arcs() ; a++)
			   if (data_S->getArcs(a,1) != data_S->getArcs(a,0))//if it is not looping arc
			     for ( int k = 0; k < data_S->getN_od(); k++)
					conUB[t++] = data_S->getC(a);
				 
			X_Value = new double*[data_S->getN_arcs()];
			for(int a=0; a< data_S->getN_arcs(); a++)
				X_Value[a] = new double[data_S->getN_od()]();	
			numColums=0;
			iteration=0;
			//(*cplex).exportModel("./Outputs/Cuts/OPTCUT_Lp_1.lp");
			(cplex).setParam(IloCplex::TiLim, 3600);
 			(cplex).setParam(IloCplex::SimDisplay, 0);
			(cplex).setParam(IloCplex::Threads, 1);
			(cplex).setParam(IloCplex::AdvInd, 0);
			(cplex).setParam(IloCplex::NumericalEmphasis, 1);
			// (cplex).setParam(IloCplex::RootAlg, IloCplex::Dual);
// 			(cplex).setParam(IloCplex::BarDisplay, 2);  						
		}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////						
		bool solve(bool infSubProblem, SubproblemCuts *cutSP, IloNumArray& columnReducedCosts, IloNumArray& varColumnVal, IloNumArray& paretoColumnVal, IloNumArray& val, IloNumArray& X, IloNumArray& reducedCosts, IloNumArray& val_reg, double ObjectiveOfScenario, Search_Param *search_param, Data_S* data_S, double* y_SOL, double** core_y_SOL, int s, int actual_sc, int n_sc) 
		{		
			bool infPapaSubProblem = false;
			iteration++;
			if(infSubProblem && search_param->getUpdateCorePoint())//if sol was infeas dont update the core
				for (int a= 0; a<data_S->getN_arcs(); a++)
				   core_y_SOL[s][a] = 2*core_y_SOL[s][a] - y_SOL[a];
			// try{
			//add the new columns		
			if(search_param->getRunModifiedMW() ){
				// cout << "MW" << endl;
				resetBounds(data_S);
				updateBounds(columnReducedCosts, varColumnVal, data_S, reducedCosts, val_reg, X);				
				updateAuxObj(cutSP, data_S, s, y_SOL, core_y_SOL, actual_sc, search_param);				
			}
			if(search_param->getRunMW()){ 
				updateMWCon(cutSP,search_param, data_S, y_SOL, ObjectiveOfScenario, actual_sc);
				updateAuxObj(cutSP, data_S, s, y_SOL, core_y_SOL, actual_sc, search_param);
			}		
			if ( (!search_param->getRunMW() && !search_param->getRunModifiedMW()) ){//papadakos approach
				updateAuxObj(cutSP, data_S, s, y_SOL, core_y_SOL, actual_sc, search_param);				
				// cout << "Papadakos"<<endl;
			}
							  
			if( !(cplex).solve() ) {
				cplex.exportModel("MW.lp");
				cout << endl << "M-W sub_Problem: " << (cplex).getStatus();				  				    				    
				//env.error() << "Failed to optimize M-W sub-problem for scenario " << actual_sc << endl;
				abort();
			}				    
			if ((cplex).getStatus() == IloAlgorithm::Unbounded ){			    
			      cout << endl << "M-W sub_Problem unbounded";
			      abort();
			}				    
			//determine if it was infeasibile
			if(!search_param->getRunModifiedMW() && !search_param->getRunMW()){
				(cplex).getDuals(X, con);
				infPapaSubProblem = checkFeasOfSol(data_S, X,actual_sc);
			}
 			  				
			(cplex).getValues(val, var);
// 			paretoColumnVal.clear();
			if(auxVarColumn.getSize()>0)
			  (cplex).getValues(paretoColumnVal, auxVarColumn);
			// }
			// catch (IloException& e) {
			    // cerr << e;
 			     // e.end();
 			    // throw;
 			// }
			 // cout << endl << "diff: " << ObjectiveOfScenario - cplex.getObjValue() << endl;
			 
			return infPapaSubProblem;
		}

/////////////////////////////////////////////////////////
		void resetBounds(Data_S* data_S)
		{
			int t=0;
			for ( int a = 0; a < data_S->getN_arcs() ; a++)
				   if (data_S->getArcs(a,1) != data_S->getArcs(a,0))//if it is not looping arc
					 for ( int k = 0; k < data_S->getN_od(); k++)
								conLB[t++] = -IloInfinity;
			// (con).setBounds(conLB, conUB);
			t=0;			
			for(int i = 0; i < data_S->getN_nodes(); i++)
			  for(int k = 0; k < data_S->getN_od(); k++){		//
				 varUB[t] = IloInfinity;//Dual_Flow_Con[i][k];
				 varLB[t++] = -IloInfinity;
			  }
			for(int a = 0; a < data_S->getN_arcs(); a++){		//
				varUB[t] = IloInfinity;//Dual_Flow_Con[i][k];
				 varLB[t++] = 0;//var[t++] =Dual_Capacity_Con[a];
				 for(int k = 0; k < data_S->getN_od(); k++){	//
					varUB[t] = IloInfinity;//Dual_Flow_Con[i][k];
					varLB[t++] = 0;//var[t++] =Dual_Strong_Ineq[a][k];
				 }					 
			}
			//resetting the bound for the new columns
			colLB.clear();
			colUB.clear();
			for(int i=0; i< auxVarColumn.getSize(); i++){
			  colLB.add(0);
			  colUB.add(IloInfinity);
			}
			// (var).setBounds(varLB, varUB);
		}
///////////////////////////////////
		void updateBounds(IloNumArray& columnReducedCosts, IloNumArray& varColumnVal, Data_S* data_S, IloNumArray& reducedCosts, IloNumArray& val_reg, IloNumArray& X)
		{ 
					for(int i=0; i<reducedCosts.getSize(); i++)
						// if()//should not fix the basix variables because they may change their value
						if(reducedCosts[i] < -1e-6 || reducedCosts[i] > 1e-6 ){// reduced cost of zero							 
								varLB[i] = val_reg[i];
								varUB[i] = val_reg[i];
							}
						// else
							// cout << reducedCosts[i] << endl;
					 (var).setBounds(varLB, varUB);
					 //the new columns
					 for(int i=0; i<columnReducedCosts.getSize(); i++)
						if(columnReducedCosts[i] > 1e-6 || columnReducedCosts[i] < -1e-6){
							colLB[i] = varColumnVal[i];
							colUB[i] = varColumnVal[i];
						}
					if(columnReducedCosts.getSize()>0)
						(auxVarColumn).setBounds(colLB, colUB);
					 int t=0;		
					 for( int a = 0; a < data_S->getN_arcs() ; a++)
					   if(data_S->getArcs(a,1) != data_S->getArcs(a,0))//if it is not looping arc
						 for( int k = 0; k < data_S->getN_od(); k++){
								if(X[t] > 1e-9 || X[t] < -1e-9)//note slack reduced cost are negative of the dual values
									conLB[t] = data_S->getC(a);							
								t++;
						}
						(con).setBounds(conLB, conUB);
		}
///////////////////////////////////////////
		void updateAuxObj(SubproblemCuts *cutSP, Data_S* data_S, int s, double* y_SOL, double** core_y_SOL, int actual_sc, Search_Param *search_param)
		{
			 //to know if we are at the LP phase
			 for (int a= data_S->getN_arcs() - data_S->getN_od(); a<data_S->getN_arcs(); a++){
			   core_y_SOL[s][a] = 1;
			   y_SOL[a] 		= 1;
			 }
			 // for (int a= 0; a<data_S->getN_arcs()- data_S->getN_od(); a++)
			   // yCore[a] = 0.3*yCore[a]+0.7*y_SOL[a];
			 IloExpr expr1(env);
			 for(int k = 0; k < data_S->getN_od(); k++)
			   expr1 += fabs(data_S->getD(k,actual_sc)) * Dual_Flow_Con[data_S->getOd(k,1)][k];//doing it in the expr way is always faster than setcoef			 
			 for (int a = 0; a < data_S->getN_arcs() ; a++)
			   if (data_S->getArcs(a,1) != data_S->getArcs(a,0)){//if it is not looping arc			     					
					expr1 -=  fabs(data_S->getU(a)*core_y_SOL[s][a]) * Dual_Capacity_Con[a];
					if(search_param->getStrong_Ineq())
					  for (int k = 0; k < data_S->getN_od(); k++)
						expr1 -= fabs(core_y_SOL[s][a]*data_S->getD(k,actual_sc)) * Dual_Strong_Ineq[a][k];				  					
			     }
			 //the obj of the new variables
			for(int i=0; i< auxVarColumn.getSize(); i++){
			    double coeffi = cutSP->getFixCoeff(i);//fixed part of the coefficient			    
			    for(int a=0; a< data_S->getN_arcs(); a++)
			      coeffi +=  cutSP->getArcCoeff(i,a) * core_y_SOL[s][a];
 			    expr1 +=  coeffi * auxVarColumn[i];			    
			}
						
			 objFun.setExpr(expr1);
			 expr1.end();
		}
///////////////////////////////////////
		bool checkFeasOfSol(Data_S* data_S, IloNumArray& X, int actual_sc)
		{		
			bool infSubProblem=false;
			int t =0;
			for (int a = 0; a < data_S->getN_arcs() ; a++)
			   if (data_S->getArcs(a,1) != data_S->getArcs(a,0))
				for (int k = 0; k < data_S->getN_od(); k++)				      
					X_Value[a][k] = X[t++];//(cplex).getDual(con[t++]);
			//check if the solution was infeasible (any flow on a dummy arc)
			for (int a = data_S->getN_arcs() - data_S->getN_od(); a< data_S->getN_arcs() ; a++){						
			    if (data_S->getArcs(a,1) != data_S->getArcs(a,0))
				  for (int k = 0; k < data_S->getN_od(); k++)
					  if(X_Value[a][k]>1e-5){
						  infSubProblem=true;//sp was infeasible
       				  // cout << data_S->getD(k,actual_sc) << " --inf SP for scenario " << actual_sc << endl;
						  break;
					  }
				if(infSubProblem)
					  break;
			}
			
			return infSubProblem;
		}
//////////////////////////////////////
		void updateMWCon(SubproblemCuts *cutSP, Search_Param *search_param, Data_S* data_S, double* y_SOL, double ObjectiveOfScenario, int actual_sc)
		{		
			 IloExpr expr1(env);
			 for(int k = 0; k < data_S->getN_od(); k++)
				expr1 += fabs(data_S->getD(k,actual_sc)) * Dual_Flow_Con[data_S->getOd(k,1)][k];//doing it in the expr way is always faster than setcoef			
			 for (int a = 0; a < data_S->getN_arcs() ; a++)
				if (data_S->getArcs(a,1) != data_S->getArcs(a,0)){//if it is not looping arc
					expr1 -=  fabs(data_S->getU(a) * y_SOL[a]) * Dual_Capacity_Con[a];
					if(search_param->getStrong_Ineq())//if we have strong ineq in the sub-problems
					  for (int k = 0; k < data_S->getN_od(); k++)
						expr1 -= fabs(y_SOL[a] * min(data_S->getD(k,actual_sc),data_S->getU(a))) * Dual_Strong_Ineq[a][k];		      
				}
				
			//the obj of the new variables
			for(int i=0; i< auxVarColumn.getSize(); i++){
			    double coeffi = cutSP->getFixCoeff(i);//fixed part of the coefficient			    
			    for(int a=0; a< data_S->getN_arcs(); a++)
			      coeffi +=  cutSP->getArcCoeff(i,a) * y_SOL[a];
 			    expr1 +=  coeffi * auxVarColumn[i];			    
			}
			
			MW_Eq_Con.endElements();
			(MW_Eq_Con).add(IloRange(env, ObjectiveOfScenario - 1e-3, expr1 , ObjectiveOfScenario + 1e-3) );
			//(model).remove(MW_Eq_Con);
			(model).add(MW_Eq_Con);
			expr1.end();
		}
/////////////////////////////////
		void addColumns(SubproblemCuts *cutSP, Data_S* data_S)
		{		  
		  for(int i=numColums; i< numColums + cutSP->getNumNewCuts(); i++){
			//obj
			double objCoeff= cutSP->getFixCoeff(i);
			for(int a=0; a< data_S->getN_arcs(); a++)
			      objCoeff +=  cutSP->getArcCoeff(i,a) ;			
			IloNumColumn col = objFun(objCoeff);
			//constraint
			int t=0;
			for ( int a = 0; a < data_S->getN_arcs() ; a++ )
			    if (data_S->getArcs(a,1) != data_S->getArcs(a,0))//if it is not looping arc
			      for ( int k = 0; k < data_S->getN_od(); k++){
				col += (con)[t](cutSP->getConCoeff(i,t));
				t++;
			      }
			    
			  auxVarColumn.add(IloNumVar(col, 0, IloInfinity));	 			 
			  col.end();
		    }	
		    
		    (model).add(auxVarColumn);
		    numColums += cutSP->getNumNewCuts();
		}

////////////////////////////////////////////////////////////////////////

		~AuxCutModel()
		{
			(cplex).end();
			(model).end();			
			(con).end();
			env.end();			
		}

};
#endif
