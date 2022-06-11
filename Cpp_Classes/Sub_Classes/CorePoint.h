#ifndef CLASS_CorePoint_H
#define CLASS_CorePoint_H


#include <ilcplex/ilocplex.h>

ILOSTLBEGIN
#define Comm COMM_WORLD

class class_CorePoint
{ 
  private:    
		IloEnv env;		//environment 
		IloModel model;
		IloCplex cplex;
		IloObjective obj;		
 		IloNumVarArray2 x;
 		IloNumVarArray  y;		
		IloNumArray ySol;
		IloRangeArray *con, *fixCon;
  public:        
		class_CorePoint()		
		{			  
		}	
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		void Model(int n_sc, Search_Param	*search_param, Data_S *data_S) //constructor
		{
			model 		= IloModel(env);
			cplex 		= IloCplex(model);
			con   		= new IloRangeArray(env);
			fixCon 		= new IloRangeArray(env);
			obj   		= IloMinimize(env);			
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
			  // for(int k =0; k < data_S->getN_od(); k++)
			    (model).add(x[a]);
			// ***** OBJECTIVE FUNCTION****
			IloExpr expr(env);
			for(int a = 0; a < data_S->getN_arcs() /* - data_S->getN_od() */; a++){
			  expr += data_S->getF(a)*y[a];			  
			  for(int k = 0; k < data_S->getN_od(); k++)			//the second stage objective for the global scenarios 
			      // for(int s = 0; s < n_sc + search_param->getMaster_sc(); s++){			     
					expr +=  data_S->getC(a) * x[a][k];//(1/master_sc) * c[a] * x[a][k][s];//
			      
			}
			
			obj.setExpr(expr);
			(model).add(obj);
			expr.end();
			
			
			// CONSTRAINT 1		
			// for(int s = 0; s < n_sc + search_param->getMaster_sc(); s++)
			{		// flow conservation constraints for global scenarios			
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
							double maxD=0;
							for(int s=0; s< n_sc + search_param->getMaster_sc(); s++)
								if(data_S->getD(k,s) > maxD)
									maxD= data_S->getD(k,s);
							// cout << "origin " << k << " demand: " << maxD << endl;
						    (*con).add(IloRange(env, maxD, expr+expr1, maxD ));//for the origin only outgoing part is neccessary						
						  }
						 else
						   if (i == data_S->getOd(k,1)){//if destination
							double maxD=0;
							for(int s=0; s<n_sc + search_param->getMaster_sc(); s++)
								if(data_S->getD(k,s) > maxD)
									maxD= data_S->getD(k,s);
							(*con).add(IloRange(env, -maxD, expr+expr1, -maxD));						
						   }
						   else
							(*con).add(IloRange(env , 0, expr + expr1, 0));						     
						  expr.end();
						  expr1.end();
					}
			}
			
			// CONSTRAINT 2
			// for(int s = 0; s < n_sc+search_param->getMaster_sc(); s++)
			{		//capacity constraints for global scenarios 				
				for(int a = 0; a < data_S->getN_arcs()- data_S->getN_od(); a++){
					IloExpr expr(env);					
					for(int k = 0; k < data_S->getN_od(); k++)
						expr += x[a][k];					
					expr -= data_S->getU(a)*y[a];
					(*con).add(IloRange(env, -IloInfinity, expr, 0));
					expr.end();
				}
			}
			
			// CONSTRAINTS		//I think he is forcing the dummy arcs to be open ... but when their obj is zero ... it is redundent 			
			// int t = data_S->getN_arcs() - data_S->getN_od();
			// for(int k = 0; k < data_S->getN_od(); k++)			
			  // (*con).add(IloRange(env, 1, y[t + k], 1));
					
			
			
			(model).add((*con));
			(model).add((*fixCon));
			
			
			(cplex).setParam(IloCplex::TiLim, 3600);
			(cplex).setParam(IloCplex::ClockType, 1);//to set it equal to cpu time since it is seq
			(cplex).setParam(IloCplex::SimDisplay, 0);			 					
			(cplex).setParam(IloCplex::EpGap, 0.03);
		    (cplex).setParam(IloCplex::Param::Threads, 1); 
		}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		void ConvertLPtoMIP()
		{		     
		      (model).add(IloConversion(env, y, ILOINT)); 
		}
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////				
		double Solve()
		{          		
			if( !(cplex).solve() ){			  
			  cout << (cplex).getStatus() << endl;
			  env.error() << "Failed to optimize restricted heuristic\n";
			  return 1;
			}		
			ySol = IloNumArray(env);
			double m_objval = (cplex).getObjValue();
			(cplex).getValues(ySol, y);
			// env.out() << y_SOL1 <<endl;
			return m_objval; 
		}
/////////////////////////////////////////////////
		double getySol(int a) {return max(0.5, ySol[a]);}
////////////////////////////////////////////////////////////////////////////////////////
		~class_CorePoint()
		{
		}


};//WND class



#endif