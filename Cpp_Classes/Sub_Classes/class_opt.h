#ifndef CLASS_OPT_H
#define CLASS_OPT_H

#include <cmath>
#include <ilcplex/ilocplex.h>
#include "./class_SP_cuts.h"
#include "./class_opt_aux.h"
#include "./class_convex_weight.h"
ILOSTLBEGIN


typedef IloArray<IloNumVarArray> IloNumVarArray2;
typedef IloArray<IloNumVarArray2> IloNumVarArray3;
typedef IloArray<IloNumVarArray3> IloNumVarArray4;

class OptimalcutModel
{ 
	private:

		IloEnv env;
		IloModel model;
		IloCplex cplex;
		IloRangeArray con, con1; //this is for the main constraints (flow conservation and capacity) constraints				
		IloNumVarArray Dual_Capacity_Con;		
		IloNumVarArray2  Dual_Flow_Con, Dual_Strong_Ineq;
			
		IloNumVarArray var, varColumn;
		IloNumArray varColumnVal, paretoColumnVal;
		IloNumArray val, val_reg, reducedCosts, columnReducedCosts;
		IloNumArray X;		 		
		IloObjective objFun;
		AuxCutModel *auxSP;
		SubproblemCuts *cutSP;
		// class_ConvW *convW;
		int iteration;
		double mio;//1e-7;//for the maximal nondominated method
		int numColums;
		int numInfeas;
		bool disactiveSIs;
		//**************************		
	public:
		//create the sub-problems
		OptimalcutModel(Data_S *data_S, Search_Param *search_param) 
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
			
// 			  			
			Dual_Capacity_Con = IloNumVarArray(env, data_S->getN_arcs(), 0, IloInfinity);	//for the capacity constraint 
			//			IloNumVarArray Dual_Capacity_Con(env, n_arcs, 0, IloInfinity);	//for the capacity constraint 
			//			this->Dual_Capacity_Con = Dual_Capacity_Con;
			// (model).add(Dual_Capacity_Con);
						
			Dual_Strong_Ineq = IloNumVarArray2(env, data_S->getN_arcs());
			//			IloNumVarArray2 Dual_Strong_Ineq(env, n_arcs);
			for(int a = 0; a < data_S->getN_arcs(); a++)
			  Dual_Strong_Ineq[a] = IloNumVarArray(env, data_S->getN_od(), 0, IloInfinity); // for the strong ineq 
			//			this->Dual_Strong_Ineq = Dual_Strong_Ineq;
			
 			/* for(int a = 0; a < data_S->getN_arcs(); a++)
 			    (model).add(Dual_Strong_Ineq[a]); */
// 			  			


			
			//for getting extreme ray of deflected sub-problem and regular one respectively 
			varColumn = IloNumVarArray(env);
			varColumnVal = IloNumArray(env);
			columnReducedCosts = IloNumArray(env);
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
			val 	= IloNumArray(env, data_S->getN_nodes()*data_S->getN_od() +  data_S->getN_arcs() + data_S->getN_arcs()*data_S->getN_od());//for the Papadakos SP
			val_reg = IloNumArray(env, data_S->getN_nodes()*data_S->getN_od() +  data_S->getN_arcs() + data_S->getN_arcs()*data_S->getN_od());//for the daul values of the regular SP
			reducedCosts = IloNumArray(env, data_S->getN_nodes()*data_S->getN_od() +  data_S->getN_arcs() + data_S->getN_arcs()*data_S->getN_od());//for the daul values of the regular SP
			paretoColumnVal=IloNumArray(env);
			
			 
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
				(con1).add(IloRange(env, 0, Dual_Capacity_Con[t + k], 0));  */
			//adding the strong inequalities to the sub-problem formulation
			if(!search_param->getStrong_Ineq())//if it is false we remove the dual variables associate to strong inequalities (fix them at zero)
			{			  
			   for ( int a = 0; a < data_S->getN_arcs() ; a++)
				for ( int k = 0; k < data_S->getN_od(); k++)
				  (con1).add(IloRange(env, 0,  Dual_Strong_Ineq[a][k] ,0));			 
			}

		      (model).add((con));
		      (model).add((con1));						
			
			//addColumn(data_S);
			
 			X = IloNumArray(env,(con).getSize());			
			auxSP = new AuxCutModel(data_S, search_param);
			cutSP = new SubproblemCuts() ;
			// convW = new class_ConvW(data_S);
			iteration=0;
			mio = search_param->getMio();
			numColums=0;
			numInfeas=0;
			disactiveSIs=false;
			//(*cplex).exportModel("./Outputs/Cuts/OPTCUT_Lp_1.lp");
			(cplex).setParam(IloCplex::TiLim, 3600);
 			(cplex).setParam(IloCplex::SimDisplay, 0);
			(cplex).setParam(IloCplex::Threads, 1);
			(cplex).setParam(IloCplex::NumericalEmphasis, 1);
// 			(cplex).setParam(IloCplex::AdvInd, 0);
			 // (cplex).setParam(IloCplex::RootAlg, IloCplex::Primal);
// 			(cplex).setParam(IloCplex::BarDisplay, 2);  						
		}
////////////////////////////////////////
		void DropSIs(Data_S *data_S)
		{			  
			for ( int a = 0; a < data_S->getN_arcs() ; a++)
				for ( int k = 0; k < data_S->getN_od(); k++)
				  (model).add(IloRange(env, 0,  Dual_Strong_Ineq[a][k] ,0));			 
			
		}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////						
		double solve(bool Phase_I, int  *sc_worker, bool &infSubProblem, bool &infPapaSubProblem, Search_Param *search_param, Data_S* data_S, double **X_Value, double* y_SOL, double** core_y_SOL, int s, int actual_sc, int n_sc, double* G, double* g_Values) 
		{
			// if(!disactiveSIs && !Phase_I && search_param->getStrong_Ineq()){
				// DropSIs(data_S);
				// disactiveSIs=true;
				// cout<< "\n disactivating for scenario " << s << endl;
			// }
			double ObjectiveOfScenario;
						
			iteration++;
		// try{
			//****************************** SP ****************************************
			updateRegObj(data_S, mio, s, y_SOL, core_y_SOL, actual_sc, search_param);
			
			//generating cuts for the SP			
			if( !(cplex).solve() ){				  
				  cout << endl << "regular sub_Problem is " << (cplex).getStatus() << endl;
				  env.error() << "Failed to optimize sub-problem for scenario " << actual_sc << endl;
				  abort();
			}
			
			if(Phase_I){//you can not  generate flow cuts when the y is int
			      (cplex).getDuals(X, con);			      
			      int t =0;
			      for (int a = 0; a < data_S->getN_arcs() ; a++)
					if (data_S->getArcs(a,1) != data_S->getArcs(a,0))
				      for (int k = 0; k < data_S->getN_od(); k++)				      
					      X_Value[a][k] = X[t++];
			      if( !search_param->getStrong_Ineq() &&  search_param->getAddResidualCuts() )//if we have strong ineq, no residual cut can be extracted
					cutSP->genResidualCuts(y_SOL, X_Value, data_S, actual_sc);
			      if(search_param->getFlowPackIneq())
					cutSP->genFPI(n_sc, data_S, y_SOL, X_Value, actual_sc);
			      //solve the SP again if new cuts were generated
			      if( addColumn( data_S, y_SOL) )
					if( !(cplex).solve() ){				  
						cout << endl << "regular sub_Problem is " << (cplex).getStatus() << endl;
						env.error() << "Added columns and failed to optimize sub-problem for scenario " << actual_sc << endl;
						abort();
					}
			}
			
			
			
			
			//if no new column is added, we donot need to resolve the SP
			if( (cplex).getStatus() == IloAlgorithm::Unbounded ){
			    cout  << "Regular sub_Problem unbounded for scenario " << actual_sc << endl;
			    cout << "number of added columns " << varColumn.getSize() << endl;
			    cplex.exportModel("aaa.lp");
			    abort();
			}
			else{
				  ObjectiveOfScenario = (cplex).getObjValue();
				  (cplex).getReducedCosts(reducedCosts, var);
				  (cplex).getValues(val_reg, var);//to get the dual values of regular SP
				  //new col
				  if( varColumn.getSize() >0 ){
				    (cplex).getReducedCosts(columnReducedCosts, varColumn);
				    (cplex).getValues(varColumnVal, varColumn);
				  }
				  //correcting the obj to get real cost
				  int t=0;
 				  if(mio > 0)
 					correctObjValue(s, data_S, core_y_SOL, ObjectiveOfScenario, mio, actual_sc);
				 
				(cplex).getDuals(X, con);
				infSubProblem = checkFeasOfSol(data_S, X_Value, actual_sc);
				if(infSubProblem){
					numInfeas++;
					// cout << "Num of infeasible iterations of scenario  " <<  actual_sc << " is: " << numInfeas << endl;
				}
				
				//*******************************M-W SPs*****************
				//to solve the SP for generating the Pareto optimal cut	
  				if( (mio > 1e-75) || (search_param->getWarmStart() && iteration < search_param->getWarmStartIterations() && Phase_I) ){//if we are using maximal nondominated, Papadakos, the MW or its modification and we are applying the warm start strategy, then we dont need the aux SP to be solved
				    val = val_reg;	
				    paretoColumnVal = varColumnVal;
  				}
  				else {
				    infPapaSubProblem = auxSP->solve(infSubProblem, cutSP, columnReducedCosts, varColumnVal, paretoColumnVal, val, X, reducedCosts, val_reg, ObjectiveOfScenario, search_param, data_S, y_SOL, core_y_SOL, s, actual_sc,  n_sc) ;					
				}
				
				fillCutVector(search_param, n_sc, data_S, G, g_Values, actual_sc);				 
			}			   
		/* }
		  catch (IloException& e) {
		    cerr << e;
		     e.end();
		    throw;
		} */		  
			return data_S->getP(actual_sc) * ObjectiveOfScenario;
		}
/////////////////////////////////////////
		void fillCutVector(Search_Param *search_param, int n_sc, Data_S* data_S, double* G, double* g_Values, int actual_sc)
		{			
				    double aux, aux1;
				    int t = 0;
				    int  tt = 0, ttt;	//a simple counter
				    int s_id = actual_sc;
				    //I wrote a loop on the scenarios because I want to copy every cut I get for a given scenario to other scenarios					   						
				    for(int s=0; s< n_sc; s++){
						aux = 0;
						ttt = 0;
						for(int i = 0; i < data_S->getN_nodes(); i++)
							for(int k = 0; k < data_S->getN_od(); k++)
								if (i == data_S->getOd(k,1))						   
								  aux += fabs(data_S->getD(k,s)) * val[ttt++]/*(cplex).getValue(Dual_Flow_Con[i][k])*/;						    						    
								else
								  ttt++;
						//for the new columns
						for(int i=0; i< varColumn.getSize(); i++)
							aux += cutSP->getScenFixPart(i, s) * paretoColumnVal[i];
						g_Values[s] = aux;	//this gives the g for every scenario based on the solution of this scenario												
					}
					//for art scenario
					aux = 0;
					ttt = 0;
					for(int i = 0; i < data_S->getN_nodes(); i++)
						for(int k = 0; k < data_S->getN_od(); k++)
							if (i == data_S->getOd(k,1))						   
							  aux += fabs(data_S->getD(k,actual_sc)) * val[ttt++]/*(cplex).getValue(Dual_Flow_Con[i][k])*/;						    						    
							else
							  ttt++;
					//for the new columns
					for(int i=0; i< varColumn.getSize(); i++)
						aux += cutSP->getFixCoeff(i) * paretoColumnVal[i];
					g_Values[n_sc] = aux;
				    //____________________________________ G[a] values _______________________											   
					t=0;
					for(int a = 0; a < data_S->getN_arcs(); a++)  {
						aux1 = 0.0;
						//for the new columns
						for(int i=0; i< varColumn.getSize(); i++)
						    aux1 +=  cutSP->getArcCoeff(i,a) * varColumnVal[i];
						G[t++] = aux1 - data_S->getU(a) * val[ttt++] ; //p[s]*
						for (int k = 0; k < data_S->getN_od(); k++)
							if(actual_sc == n_sc+search_param->getMaster_sc())//if its artificial scenario we multiple it to its demand because we are not propagate it
								G[t++] = - data_S->getD(k,actual_sc) * val[ttt++] ;
							else
								G[t++] = -val[ttt++] ;//here we add the dual values associated to the strong inequalities to the coiffiecient of the design y variables										      												      
				    }
					
					
					
					
					
				/*	
					//regular cuts from regular SP
				    aux = 0;
				    ttt =0;
				    for(int i = 0; i < data_S->getN_nodes(); i++)
					for(int k = 0; k < data_S->getN_od(); k++)
					    if (i == data_S->getOd(k,1))						   
					      aux += fabs(data_S->getD(k,s_id)) * val_reg[ttt++];						    						    
					    else
					      ttt++;
				    //for the new columns
				    for(int i=0; i< varColumn.getSize(); i++)
					aux += cutSP->getFixCoeff(i) * varColumnVal[i];
				    g_Values[1] = aux;	//this gives the g for every scenario based on the solution of this scenario					    
				    //____________________________________ G[a] values _______________________							
				    for(int a = 0; a < data_S->getN_arcs(); a++)  {
						aux1 = 0.0;
						aux1 = -data_S->getU(a) * val_reg[ttt++] ; //p[s]*
						for (int k = 0; k < data_S->getN_od(); k++)
							aux1 += -data_S->getD(k,s_id) * val_reg[ttt++] ;//here we add the dual values associated to the strong inequalities to the coiffiecient of the design y variables										      
						//for the new columns
						for(int i=0; i< varColumn.getSize(); i++)
						    aux1 +=  cutSP->getArcCoeff(i,a) * varColumnVal[i];
						G[data_S->getN_arcs() + a] = aux1; 						      
				    }*/
		}
///////////////////////////////////////////////////////
		void updateRegObj(Data_S* data_S, double mio, int s, double* y_SOL, double** core_y_SOL, int actual_sc, Search_Param *search_param)
		{			
			
			for (int a= data_S->getN_arcs() - data_S->getN_od(); a<data_S->getN_arcs(); a++)
			  y_SOL[a] = 1;
// 			int Dest_K;
			IloExpr expr(env), expr1(env);
			for(int k = 0; k < data_S->getN_od(); k++){
			  expr  += abs(data_S->getD(k,actual_sc)) * Dual_Flow_Con[data_S->getOd(k,1)][k];//doing it in the expr way is always faster than setcoef			
 			  expr1 += abs(data_S->getD(k,actual_sc)) * Dual_Flow_Con[data_S->getOd(k,1)][k];
			}
			for (int a = 0; a < data_S->getN_arcs() ; a++)
			  if (data_S->getArcs(a,1) != data_S->getArcs(a,0))//if it is not looping arc
			    {			      
			      expr  -=  fabs(data_S->getU(a)*y_SOL[a]) * Dual_Capacity_Con[a];
 			      expr1 -=  fabs(data_S->getU(a)*core_y_SOL[s][a]) * Dual_Capacity_Con[a];
			      if(search_param->getStrong_Ineq())//if we have strong ineq in the sub-problems
				for (int k = 0; k < data_S->getN_od(); k++){
				  expr  -= fabs(y_SOL[a]*data_S->getD(k,actual_sc)) * Dual_Strong_Ineq[a][k];
 				  expr1 -= fabs(core_y_SOL[s][a]*data_S->getD(k,actual_sc)) * Dual_Strong_Ineq[a][k];
				}
			    }
			//the obj of the new variables
			for(int i=0; i< varColumn.getSize(); i++){
			  double coeffi = cutSP->getFixCoeff(i);//fixed part of the coefficient
 			  double coeffi1 = cutSP->getFixCoeff(i);
			  for(int a=0; a< data_S->getN_arcs(); a++){//the y related part
			    coeffi +=  cutSP->getArcCoeff(i,a) * y_SOL[a];
 			    coeffi1 +=  cutSP->getArcCoeff(i,a) * core_y_SOL[s][a];
			  }
			  expr +=  coeffi * varColumn[i];
 			  expr1 +=  coeffi1 * varColumn[i];
			}
			objFun.setExpr(expr + mio*expr1);			
			expr.end();
 			expr1.end();

		}
/////////////////////////////////////////
		void correctObjValue(int s, Data_S* data_S, double** core_y_SOL, double& ObjectiveOfScenario, double mio, int actual_sc)
		{
			double goalProgrammingCost=0;
			int t=0;			
			for(int k = 0; k < data_S->getN_od(); k++)
			   for(int i = 0; i < data_S->getN_nodes(); i++)
				if(i == data_S->getOd(k,1))
				  goalProgrammingCost += abs(data_S->getD(k,actual_sc)) * val_reg[t++];//var[t++] =Dual_Flow_Con[i][k];
				else
				    t++;
			  for(int a = 0; a < data_S->getN_arcs(); a++)
			     if (data_S->getArcs(a,1) != data_S->getArcs(a,0)){
				  goalProgrammingCost -=  fabs(data_S->getU(a)*core_y_SOL[s][a]) * val_reg[t++];//var[t++] =Dual_Capacity_Con[a];
				  for(int k = 0; k < data_S->getN_od(); k++)	//
				    goalProgrammingCost -= fabs(core_y_SOL[s][a]*data_S->getD(k,actual_sc)) * val_reg[t++];//var[t++] =Dual_Strong_Ineq[a][k];			    
			     }
			     else
					t++;
			   ObjectiveOfScenario -= mio*goalProgrammingCost;
		}
///////////////////////////////////////
		bool checkFeasOfSol(Data_S* data_S, double **X_Value, int actual_sc)
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
					  if(X_Value[a][k]>1e-2){
						  infSubProblem=true;//sp was infeasible
						// cout << data_S->getD(k,actual_sc) << " inf SP for scenario " <<  actual_sc << " -> " << X_Value[a][k] << endl;
						  break;
					  }
				if(infSubProblem)
					  break;
			}
			return infSubProblem;
		}
//////////////////////////////////////
		bool addColumn(Data_S* data_S, double* y_SOL)
		{
		    bool newColAdded=false;
		    for(int i=numColums; i< numColums + cutSP->getNumNewCuts(); i++){
				//obj
				double objCoeff= cutSP->getFixCoeff(i);
				for(int a=0; a< data_S->getN_arcs(); a++)
					  objCoeff +=  cutSP->getArcCoeff(i,a) * y_SOL[a];			
				IloNumColumn col = objFun(objCoeff);
				//constraint
				int t=0;
				for ( int a = 0; a < data_S->getN_arcs() ; a++ )
					if (data_S->getArcs(a,1) != data_S->getArcs(a,0))//if it is not looping arc
					  for ( int k = 0; k < data_S->getN_od(); k++){
						col += (con)[t](cutSP->getConCoeff(i,t));
						t++;
					  }			    
				  varColumn.add(IloNumVar(col, 0, IloInfinity));		 
				  col.end();
		    }
		    
		    numColums += cutSP->getNumNewCuts();
		   		      
		    if(cutSP->getNumNewCuts()>0){
		      (model).add(varColumn);
		       auxSP->addColumns(cutSP, data_S);//add col to the aux SP
			cutSP->cleanCutMemory();
			newColAdded = true;
		    }
		    //IloNumVar var(RollsUsed(1) + Fill(newPatt), 0, MAXCUT);
		    return newColAdded;
		}
////////////////////////////////////
		int getNumInfeas() {return numInfeas;}
		~OptimalcutModel()
		{
			(cplex).end();
			(model).end();			
			(con).end();
			env.end();			
		}

};
#endif
