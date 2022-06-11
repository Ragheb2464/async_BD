#ifndef CLASS_OptRound_H
#define CLASS_OptRound_H


#include <ilcplex/ilocplex.h>
ILOSTLBEGIN



typedef IloArray<IloNumVarArray> IloNumVarArray2;
typedef IloArray<IloNumVarArray2> IloNumVarArray3;
typedef IloArray<IloNumVarArray3> IloNumVarArray4;
typedef IloArray<IloNumVarArray4> IloNumVarArray5;



class OptRound
{ 
  private:		
		IloEnv env;		//environment 
		IloModel model;
		IloCplex cplex;
		IloObjective obj;
		
 		IloNumVarArray3 x;
 		IloNumVarArray  y;
		IloNumArray y_SOL;
		IloNumArray yLB, yUB;
		IloRangeArray con;
		double fracSolSum;
		vector <int> scenarioID;
			
  public:
		OptRound(PoolsandMangers* Managers, Data_S* data_S, Search_Param* search_param, int n_sc) 
		{
			model 		=  IloModel(env);
			cplex 		=  IloCplex(model);
			con   		=  IloRangeArray(env);
			obj   		=  IloMinimize(env);
			
			
			
			
			
			double maxVal=0.0, aux;
			int maxID;
			for(int s=0; s<n_sc; s++){
				aux =0.0;
				for(int k = 0; k < data_S->getN_od(); k++)
					aux += data_S->getD(k,s);
				if(aux > maxVal){
					maxVal= aux ;
					maxID = s;
				}
			}
			n_sc =min(n_sc, 1);//round(n_sc);			
			scenarioID.push_back(maxID);
			// cout << endl << "WWWWWWWWWWWWWWWWWWWWW " << maxID << endl;
			
			
			
			
			//***the master only includes y and theta variables****
			y = IloNumVarArray(env, data_S->getN_arcs(), 0, 1, ILOINT);	//ILOINT defining an variable for each arc (on dimensional variable)									
			(model).add(y); // we dont need to add variables explicitly to the model as they will implicitly be added to the model through the constraints
			//adding x variable to the master for the global scenarios 
			x = IloNumVarArray3(env, data_S->getN_arcs());//second-stage (x) variable (three dimensional cuz each arc (i-j) is represented with a single name "a"
			for(int a = 0; a < data_S->getN_arcs(); a++){
				x[a] = IloNumVarArray2(env, data_S->getN_od());
				for(int k = 0; k < data_S->getN_od(); k++)				
					x[a][k] = IloNumVarArray(env, n_sc, 0, IloInfinity); 									
			}
// 			this->x=x;
			for(int a = 0; a < data_S->getN_arcs(); a++)
			  for(int k = 0; k < data_S->getN_od(); k++)
			    (model).add(x[a][k]);
			// ***** OBJECTIVE FUNCTION****			
			int ss;
			IloExpr expr(env);
			for(int a = 0; a < data_S->getN_arcs(); a++){
			  expr +=  /* fabs(data_S->getF(a)) * */ y[a];			  
			    // for(int k = 0; k < data_S->getN_od(); k++)			//the second stage objective for the global scenarios 
			      // for(int s = 0; s < n_sc; s++)			     
					  // expr +=  fabs(data_S->getC(a)/n_sc) * x[a][k][s];//(1/master_sc) * c[a] * x[a][k][s];//			       
			}			
			obj.setExpr(expr);
			(model).add(obj);
			expr.end();		
			// CONSTRAINT 1		
			for(int s = 0; s < n_sc; s++)	{	// flow conservation constraints for global scenarios			
				ss = scenarioID[s];
				for(int k = 0; k < data_S->getN_od(); k++) 
					for(int i = 0; i < data_S->getN_nodes(); i++){ 
						  IloExpr expr(env);//expr.clear();
						  IloExpr expr1(env);
						  for(int a =0; a < data_S->getN_arcs(); a++) {
							  if(data_S->getArcs(a,0) == i) 
								  expr += x[a][k][s];								  
							  if(data_S->getArcs(a,1) == i) 
								  expr1 -=  x[a][k][s];
						 }
						 if (i == data_S->getOd(k,0))//if origin
						    (con).add(IloRange(env, fabs(data_S->getD(k,ss)), expr+expr1, fabs(data_S->getD(k,ss))));//for the origin only outgoing part is neccessary						
						 else if (i == data_S->getOd(k,1))//if destination
							(con).add(IloRange(env, -fabs(data_S->getD(k,ss)),  expr+expr1, -fabs(data_S->getD(k,ss))));						
						else //(i != data_S->getOd(k,1) && i != data_S->getOd(k,0))
							  (con).add(IloRange(env, 0, expr + expr1, 0));					 							  
						expr.end();
						expr1.end();
					}
			}		
			// CONSTRAINT 2
			for(int s = 0; s < n_sc; s++){		//capacity constraints for global scenarios 				
				ss = scenarioID[s];// Managers->getglobal_scenario_id(s);//global_scenario_id[s];
				for(int a = 0; a < data_S->getN_arcs(); a++){
					IloExpr expr(env);					
					for(int k = 0; k < data_S->getN_od(); k++)
						expr += x[a][k][s];					
					expr -= data_S->getU(a)*y[a];
					(con).add(IloRange(env,-IloInfinity,expr,0));
					expr.end();
				}
			}

			for(int s = 0; s < n_sc; s++)	 	
				for(int k = 0; k < data_S->getN_od(); k++)		
				    for(int a = 0; a < data_S->getN_arcs(); a++){
							double aux= min(fabs(data_S->getU(a)), fabs(data_S->getD(k,scenarioID[s]))) ;								
							(con).add(IloRange(env, -IloInfinity, x[a][k][s] - aux * y[a], 0));													
					}	

			
			// CONSTRAINTS		//I think he is forcing the dummy arcs to be open ... but when their obj is zero ... it is redundent 			
			int t = data_S->getN_arcs() - data_S->getN_od();
			for(int k = 0; k < data_S->getN_od(); k++){
			  IloExpr expr(env);
			  expr = y[t + k];
			  (con).add(IloRange(env, 0, expr, 0));
			  expr.end();
			} 	

			(model).add(con);	
			
			y_SOL = IloNumArray(env, data_S->getN_arcs());
			yLB = IloNumArray(env, data_S->getN_arcs());
			yUB = IloNumArray(env, data_S->getN_arcs());
			for(IloInt a = 0; a < data_S->getN_arcs(); a++){
				yLB[a] = 0 ;
				yUB[a] = 1;			      
			}
			
			
			cplex.setOut(env.getNullStream());
			(cplex).setParam(IloCplex::RandomSeed, 2015);//fix cplex random seed
			(cplex).setParam(IloCplex::TiLim, 600);
			// (cplex).setParam(IloCplex::ClockType, 1);//to set it equal to cpu time since it is seq
			(cplex).setParam(IloCplex::SimDisplay, 0);			 					
			(cplex).setParam(IloCplex::EpGap, 0.005);
			(cplex).setParam(IloCplex::AdvInd, 0);
			(cplex).setParam(IloCplex::Param::Threads, 1); 
			   
		}
//////////////////////////////////////////////
		bool FixVars(IloNumArray& setSolVal, Data_S* data_S)
		{
			// env.out() << setSolVal<<endl;
			bool solveMip=false;
			fracSolSum=0.0;
			// y.setBounds(yLB, yUB);
			// IloExpr expr(env);
			 for(int a = 0; a < data_S->getN_arcs() -data_S->getN_od() ; a++){
				yLB[a] = 0;
				yUB[a] = 1;
				if(setSolVal[a]>0.5)
					yLB[a] = 1;
				else if(setSolVal[a]<1e-7)
					yUB[a] = 0;
				// fracSolSum += setSolVal[a];
				// expr += (1-setSolVal[a])*y[a];				
			}
			 y.setBounds(yLB, yUB); 
			// obj.setExpr(expr);
			// expr.end();	
			return true;
		}
//////////////////////////////////////////////
		bool solve(IloNumArray& setSolVal, Data_S* data_S)
		{
			// (cplex).exportModel("Lp_1.lp");
			// abort();
			if(!cplex.solve()){
				// (cplex).exportModel("Lp_1.lp");
				// abort();
				return false;
			}
			if(cplex.getObjValue() - fracSolSum > 0.001){
				cout << endl << "--------- new feasible sol found = " << cplex.getObjValue()-fracSolSum << endl;
				cplex.getValues(y_SOL, y);
				// IloExpr expr(env);
				for(int a = 0; a < data_S->getN_arcs()/* -data_S->getN_od() */; a++){
					setSolVal[a] = y_SOL[a];
					// if( setSolVal[a] > 0.0001 ){
						// cout << y_SOL[a] << "; ";
						// expr += 1 - y[a];
					// }
					// else
						// expr += y[a];
				}		
				// model.add(IloRange(env, 2, expr, IloInfinity));			
				// cout << endl;
				// abort();
			}
			return true;
		}
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		~OptRound()
		{
		}


};
#endif
