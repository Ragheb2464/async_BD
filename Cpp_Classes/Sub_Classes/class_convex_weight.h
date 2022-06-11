#ifndef CLASS_Conv_W_H
#define CLASS_Conv_W_H


#include <ilcplex/ilocplex.h>

ILOSTLBEGIN
#define Comm COMM_WORLD

class class_ConvW
{ 
  private:    
		IloEnv env;		//environment 
		IloModel model;
		IloCplex cplex;
		IloObjective obj;		
 		IloNumVarArray2 x;
 		IloNumVarArray  lambda;		
 		IloNumArray  lambdaVal;		
		IloRangeArray *con, *con1;
  public:        
		class_ConvW(Data_S *data_S)		
		{	
			model 		= IloModel(env);
			cplex 		= IloCplex(model);
			con   		= new IloRangeArray(env);
			con1   		= new IloRangeArray(env);
			obj   		= IloMinimize(env);			
			//***the master only includes y and theta variables****
			lambda = IloNumVarArray(env, data_S->getN_arcs()- data_S->getN_od(), 0, 1);	//ILOINT defining an variable for each arc (on dimensional variable)	
			(model).add(lambda); // we dont need to add variables explicitly to the model as they will implicitly be added to the model through the constraints
			x = IloNumVarArray2(env, data_S->getN_arcs()- data_S->getN_od());//second-stage (x) variable (three dimensional cuz each arc (i-j) is represented with a single name "a"
			for(int a = 0; a < data_S->getN_arcs()- data_S->getN_od(); a++)
				x[a] = IloNumVarArray(env, data_S->getN_od(), 0, IloInfinity);
// 			this->x=x;
			for(int a = 0; a < data_S->getN_arcs()- data_S->getN_od(); a++)
			    (model).add(x[a]);
			// ***** OBJECTIVE FUNCTION****
			IloExpr expr(env);
			for(int a = 0; a < data_S->getN_arcs()- data_S->getN_od(); a++)
			  expr += lambda[a];
			obj.setExpr(expr);
			(model).add(obj);
			expr.end();
			
			
			// CONSTRAINT 1				// flow conservation constraints for global scenarios			
				for(int k = 0; k < data_S->getN_od(); k++)				  
					for(int i = 0; i < data_S->getN_nodes(); i++){ 
						  IloExpr expr(env);//expr.clear();
						  IloExpr expr1(env);
						  for(int a =0; a < data_S->getN_arcs()- data_S->getN_od(); a++){
							  if(data_S->getArcs(a,0) == i) 
								  expr += x[a][k];								  
							  if(data_S->getArcs(a,1) == i) 
								  expr1 -= x[a][k];
						  }
						  if (i == data_S->getOd(k,0))//if origin
						    (*con).add(IloRange(env, data_S->getD(k,0), expr+expr1, data_S->getD(k,0)));//for the origin only outgoing part is neccessary						
						 else
						   if (i == data_S->getOd(k,1))//if destination
							(*con).add(IloRange(env, -data_S->getD(k,0), expr+expr1, -data_S->getD(k,0)));						
						   else
							(*con).add(IloRange(env , 0, expr + expr1, 0));						     
						  expr.end();
						  expr1.end();
					}
			
			
			// CONSTRAINT 2
			for(int a = 0; a < data_S->getN_arcs()- data_S->getN_od(); a++){
					IloExpr expr(env);					
					for(int k = 0; k < data_S->getN_od(); k++)
						expr += x[a][k];					
					expr -= data_S->getU(a)*lambda[a];
					(*con).add(IloRange(env,-IloInfinity,expr,0));
					expr.end();
				}
			
			
			
			
			(model).add((*con));
			// (model).add((*con1));
			lambdaVal=IloNumArray(env, data_S->getN_arcs()- data_S->getN_od());
			
			(cplex).setParam(IloCplex::TiLim, 3600);
			(cplex).setParam(IloCplex::ClockType, 1);//to set it equal to cpu time since it is seq
			(cplex).setParam(IloCplex::SimDisplay, 0);			 					
			(cplex).setParam(IloCplex::EpGap, 0.03);
		    (cplex).setParam(IloCplex::Param::Threads, 1); 				 
		}
/////////////////////////////////////////////////////////////////////////////
		void solve(int actual_sc, double* yCore, double* y_SOL, Data_S *data_S)
		{
			updateCons(actual_sc, yCore, y_SOL, data_S);
			// updateObj(yCore, y_SOL, data_S);
			if(!cplex.solve())
				cout << "failed to obtain optimal lamdas" << endl;
			else{
				cplex.getValues(lambdaVal, lambda);
				// env.out() << lambdaVal << endl;
				/* for(int a = 0; a < data_S->getN_arcs(); a++)
					for(int k = 0; k < data_S->getN_od(); k++)
						if(cplex.getValue(x[a][k]) > 1e-3)
							cout << cplex.getValue(x[a][k]) << endl; */
				for(int a=0; a<data_S->getN_arcs()- data_S->getN_od(); a++)
					y_SOL[a] = (1-lambdaVal[a]-1e-6)*y_SOL[a] + (lambdaVal[a]+1e-6)*yCore[a];				
			}
		}
////////////////////////////////////////////////////////////////////////////
		void updateCons(int actual_sc, double* yCore, double* y_SOL, Data_S *data_S)
		{
			(*con).endElements();
			// (*con1).endElements();
			// CONSTRAINT 1				// flow conservation constraints for global scenarios			
			for(int k = 0; k < data_S->getN_od(); k++)				  
					for(int i = 0; i < data_S->getN_nodes(); i++){ 
						  IloExpr expr(env);//expr.clear();
						  IloExpr expr1(env);
						  for(int a =0; a < data_S->getN_arcs()- data_S->getN_od(); a++){
							  if(data_S->getArcs(a,0) == i) 
								  expr += x[a][k];								  
							  if(data_S->getArcs(a,1) == i) 
								  expr1 -= x[a][k];
						  }
						  if (i == data_S->getOd(k,0))//if origin
						    (*con).add(IloRange(env, data_S->getD(k,actual_sc), expr+expr1, data_S->getD(k,actual_sc)));//for the origin only outgoing part is neccessary						
						 else
						   if (i == data_S->getOd(k,1))//if destination
							(*con).add(IloRange(env, -data_S->getD(k,actual_sc), expr+expr1, -data_S->getD(k,actual_sc)));						
						   else
							(*con).add(IloRange(env , 0, expr + expr1, 0));						     
						  expr.end();
						  expr1.end();
					}
			
			
			// CONSTRAINT 2
			for(int a = 0; a < data_S->getN_arcs()- data_S->getN_od(); a++){
					IloExpr expr(env);					
					for(int k = 0; k < data_S->getN_od(); k++)
						expr += x[a][k];					
					expr -= data_S->getU(a)*(yCore[a]-y_SOL[a])*lambda[a] ;
					(*con).add(IloRange(env,-IloInfinity,expr,data_S->getU(a)*y_SOL[a]));
					expr.end();
			}		
		(model).add((*con));			
		}
/////////////////////////////////////////////////////////////////////////////
		void updateObj(double* yCore, double* y_SOL, Data_S *data_S)
		{
			IloExpr expr(env);
			for(int a = 0; a < data_S->getN_arcs(); a++)
			  expr += fabs(yCore[a]-y_SOL[a])*lambda[a];
			obj.setExpr(expr);
			expr.end();
		}
////////////////////////////////////////////////////////////////////////////
		~class_ConvW()
		{
		}


};//WND class



#endif