#ifndef CLASS_DEF_H
#define CLASS_DEF_H


#include <ilcplex/ilocplex.h>

ILOSTLBEGIN
#define Comm COMM_WORLD

class class_DEF
{ 
  private:    
		IloEnv env;		//environment 
		IloModel model;
		IloCplex cplex;
		IloObjective obj;		
 		IloNumVarArray3 x;
 		IloNumVarArray  y;		
		IloRangeArray *con, *fixCon;
  public:        
		class_DEF()		
		{			  
		}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		//the master processor main loop
		double main_loop(int n_sc, Search_Param *search_param, Data_S *data_S, class_MP *MP)
		{			  
			  double wtime = MPI::Wtime();
			  MasterModel(n_sc,search_param,data_S); //this creates the master problem	
 			  ConvertLPtoMIP();	//here we turn into integer those who are not
			  FixVariables(MP);
			  double ObjValue = solve_phase_I(); 			
			 cout << endl << ">>>>>>>>>>>>>>>>>>>>>>>> Obj and Time of heuristic are: " << ObjValue << " and " << MPI::Wtime() - wtime << endl;
			 return ObjValue;
		}		
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		void MasterModel(int n_sc, Search_Param	*search_param, Data_S *data_S) //constructor
		{
			model 		= IloModel(env);
			cplex 		= IloCplex(model);
			con   		= new IloRangeArray(env);
			fixCon 		= new IloRangeArray(env);
			obj   		= IloMinimize(env);			
			//***the master only includes y and theta variables****
			y = IloNumVarArray(env, data_S->getN_arcs(), 0, 1);	//ILOINT defining an variable for each arc (on dimensional variable)	
			(model).add(y); // we dont need to add variables explicitly to the model as they will implicitly be added to the model through the constraints
			x = IloNumVarArray3(env, data_S->getN_arcs());//second-stage (x) variable (three dimensional cuz each arc (i-j) is represented with a single name "a"
			for(int a = 0; a < data_S->getN_arcs(); a++){
				x[a] = IloNumVarArray2(env, data_S->getN_od());
				for(int k = 0; k < data_S->getN_od(); k++)				
					x[a][k] = IloNumVarArray(env, n_sc + search_param->getMaster_sc(), 0, IloInfinity); 									
			}
// 			this->x=x;
			for(int a = 0; a < data_S->getN_arcs(); a++)
			  for(int k =0; k < data_S->getN_od(); k++)
			    (model).add(x[a][k]);
			// ***** OBJECTIVE FUNCTION****
			IloExpr expr(env);
			for(int a = 0; a < data_S->getN_arcs() - data_S->getN_od(); a++){
			  expr += data_S->getF(a)*y[a];			  
			  for(int k = 0; k < data_S->getN_od(); k++)			//the second stage objective for the global scenarios 
			      for(int s = 0; s < n_sc + search_param->getMaster_sc(); s++){			     
					expr += data_S->getP(s) * data_S->getC(a) * x[a][k][s];//(1/master_sc) * c[a] * x[a][k][s];//
			      }
			}
			
			obj.setExpr(expr);
			(model).add(obj);
			expr.end();
			
			
			// CONSTRAINT 1		
			for(int s = 0; s < n_sc + search_param->getMaster_sc(); s++){		// flow conservation constraints for global scenarios			
				for(int k = 0; k < data_S->getN_od(); k++)				  
					for(int i = 0; i < data_S->getN_nodes(); i++){ 
						  IloExpr expr(env);//expr.clear();
						  IloExpr expr1(env);
						  for(int a =0; a < data_S->getN_arcs()- data_S->getN_od(); a++){
							  if(data_S->getArcs(a,0) == i) 
								  expr += x[a][k][s];								  
							  if(data_S->getArcs(a,1) == i) 
								  expr1 -=  x[a][k][s];
						  }
						  if (i == data_S->getOd(k,0))//if origin
						    (*con).add(IloRange(env, data_S->getD(k,s), expr+expr1, data_S->getD(k,s)));//for the origin only outgoing part is neccessary						
						 else
						   if (i == data_S->getOd(k,1))//if destination
							(*con).add(IloRange(env, -data_S->getD(k,s), expr+expr1, -data_S->getD(k,s)));						
						   else
							(*con).add(IloRange(env , 0, expr + expr1, 0));						     
						  expr.end();
						  expr1.end();
					}
			}
			
			// CONSTRAINT 2
			for(int s = 0; s < n_sc+search_param->getMaster_sc(); s++){		//capacity constraints for global scenarios 				
				for(int a = 0; a < data_S->getN_arcs()- data_S->getN_od(); a++){
					IloExpr expr(env);					
					for(int k = 0; k < data_S->getN_od(); k++)
						expr += x[a][k][s];					
					expr -= data_S->getU(a)*y[a];
					(*con).add(IloRange(env,-IloInfinity,expr,0));
					expr.end();
				}
			}
			
			// CONSTRAINTS		//I think he is forcing the dummy arcs to be open ... but when their obj is zero ... it is redundent 			
			int t = data_S->getN_arcs() - data_S->getN_od();
			for(int k = 0; k < data_S->getN_od(); k++)			
			  (*con).add(IloRange(env, 1, y[t + k], 1));
					
			
			
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
		double solve_phase_I()
		{             
			if( !(cplex).solve() ){			  
			  cout << (cplex).getStatus() << endl;
			  env.error() << "Failed to optimize restricted heuristic\n";
			  return 1;
			}			
			double m_objval = (cplex).getObjValue();
			(cplex).getValues(y_SOL, y);				
			return m_objval; 
		}
///////////////////////////////////////////////////////////////////////////////////////////
		void FixVariables(class_MP *MP)
		{
			int numFixed=0;
			MP->ReducedCosts(); 
			(*fixCon).removeFromAll();
			for(int a=0; a< MP->getreducedCostsSize(); a++){
				if(MP->getreducedCosts(a) > 0 && MP->gety_SOL(a)<0.05){
					(*fixCon).add(IloRange(env, 0, y[a], 0));
					numFixed++;
				}
				else if(MP->getreducedCosts(a) < -0 && MP->gety_SOL(a)>0.95){
					(*fixCon).add(IloRange(env, 1, y[a], 1));
					numFixed++;
				}
			}			
			// (model).add((*fixCon));
			cout << "\n num of fixed vars: " << numFixed << " out of " << reducedCosts.getSize()<< endl;
		}
////////////////////////////////////////////////////////////////////////////////////////
		~class_DEF()
		{
		}


};//WND class



#endif