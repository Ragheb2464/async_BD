#ifndef CLASS_MH_H
#define CLASS_MH_H


#include <ilcplex/ilocplex.h>
ILOSTLBEGIN



typedef IloArray<IloNumVarArray> IloNumVarArray2;
typedef IloArray<IloNumVarArray2> IloNumVarArray3;
typedef IloArray<IloNumVarArray3> IloNumVarArray4;
typedef IloArray<IloNumVarArray4> IloNumVarArray5;



class MasterHeuristic
{ 
  private:		
		IloEnv env;		//environment 
		IloModel model;
		IloCplex cplex;
		IloObjective obj;
		
 		IloNumVarArray3 x;
 		IloNumVarArray  y;
 		IloNumVarArray  theta;
		IloNumArray y_SOL;
		IloNumArray3 x_Sol;
		IloNumArray theta_SOL;
		IloRangeArray con, *y_con;
                double *Sol_send;
			
  public:
		MasterHeuristic(PoolsandMangers* Managers, Data_S* data_S, Search_Param* search_param, int n_sc) 
		{
			model 		=  IloModel(env);
			cplex 		=  IloCplex(model);
			con   		=  IloRangeArray(env);
			obj   		=  IloMinimize(env);
			y_con		=  new IloRangeArray(env);				
			//***the master only includes y and theta variables****
			y = IloNumVarArray(env, data_S->getN_arcs(), 0, 1, ILOINT);	//ILOINT defining an variable for each arc (on dimensional variable)						
			theta = IloNumVarArray(env, n_sc + search_param->getMaster_sc(), 0, IloInfinity);	//defining variables to approximate the recourse costs
											//I generate a theta for the global scenarios too, but keep their coeficient to zero (just to ease the programming pain :P)
// 			this->theta = theta;

			(model).add(y); // we dont need to add variables explicitly to the model as they will implicitly be added to the model through the constraints
			(model).add(theta);
					
			//cout << endl << "the size of master scenarios " << (*global_scenario_id).size() << endl;
			//adding x variable to the master for the global scenarios 
			x = IloNumVarArray3(env, data_S->getN_arcs());//second-stage (x) variable (three dimensional cuz each arc (i-j) is represented with a single name "a"
			for(int a = 0; a < data_S->getN_arcs(); a++){
				x[a] = IloNumVarArray2(env, data_S->getN_od());
				for(int k = 0; k < data_S->getN_od(); k++)				
					x[a][k] = IloNumVarArray(env, search_param->getMaster_sc(), 0, IloInfinity); 									
			}
// 			this->x=x;
			for(int a = 0; a < data_S->getN_arcs(); a++)
			  for(int k = 0; k < data_S->getN_od(); k++)
			    (model).add(x[a][k]);

			  
			// ***** OBJECTIVE FUNCTION****			
			int ss;
			IloExpr expr(env);
			for(int a = 0; a < data_S->getN_arcs(); a++){
			  expr += data_S->getF(a)*y[a];			  
			  for(int k = 0; k < data_S->getN_od(); k++)			//the second stage objective for the global scenarios 
			      for(int s = 0; s < search_param->getMaster_sc(); s++){			     
				  ss = Managers->getglobal_scenario_id(s);//global_scenario_id[s];//this gives the indicator of the scenario 
				  expr += data_S->getP(ss) * data_S->getC(a) * x[a][k][s];//(1/master_sc) * c[a] * x[a][k][s];//
			      }
			}
			
			 
			for(int s = 0; s < n_sc + search_param->getMaster_sc() ; s++)			  
			      expr += data_S->getP(s)* theta[s];//initially I need all theta not be in obj otherwise the master problem runns outbounded
			
			obj.setExpr(expr);

			(model).add(obj);

			expr.end();
			
			
			// CONSTRAINT 1		
			for(int s = 0; s < search_param->getMaster_sc(); s++)		// flow conservation constraints for global scenarios
			{
				ss = Managers->getglobal_scenario_id(s);//global_scenario_id[s];
				for(int k = 0; k < data_S->getN_od(); k++) 
					for(int i = 0; i < data_S->getN_nodes(); i++) 
					{ 
						  IloExpr expr(env);//expr.clear();
						  IloExpr expr1(env);
						  for(int a =0; a < data_S->getN_arcs(); a++)  
						  {
							  if(data_S->getArcs(a,0) == i) 
								  expr += x[a][k][s];								  
							  if(data_S->getArcs(a,1) == i) 
								  expr1 -=  x[a][k][s];
						  }
						  if (i == data_S->getOd(k,0))//if origin
						    (con).add(IloRange(env, data_S->getD(k,ss), expr+expr1, data_S->getD(k,ss)));//for the origin only outgoing part is neccessary						
						 else
						   if (i == data_S->getOd(k,1))//if destination
							(con).add(IloRange(env, -data_S->getD(k,ss), expr+expr1, -data_S->getD(k,ss)));						
						   else
// 						     if ( (i != data_S->getOd(k,0)) && (i != data_S->getOd(k,1)) )
							  (con).add(IloRange(env , 0, expr + expr1, 0));
						     
						  expr.end();
						  expr1.end();
					}
			}
			
			// CONSTRAINT 2
			for(int s = 0; s < search_param->getMaster_sc(); s++)		//capacity constraints for global scenarios 	
			{
				ss =  Managers->getglobal_scenario_id(s);//global_scenario_id[s];
				for(int a = 0; a < data_S->getN_arcs(); a++)
				{
					IloExpr expr(env);
					
					for(int k = 0; k < data_S->getN_od(); k++)
						expr += x[a][k][s];
					
					expr -= data_S->getU(a)*y[a];
					(con).add(IloRange(env,-IloInfinity,expr,0));
					expr.end();
				}
			}
			
			// CONSTRAINTS		//I think he is forcing the dummy arcs to be open ... but when their obj is zero ... it is redundent 
			
			int t = data_S->getN_arcs() - data_S->getN_od();
			for(int k = 0; k < data_S->getN_od(); k++)
			{
			  IloExpr expr(env);
			  expr = y[t + k];
			  (con).add(IloRange(env, 1, expr, 1));
			  expr.end();
			}							      
			
			(model).add(con);	
			
			y_SOL = IloNumArray(env, data_S->getN_arcs());
			x_Sol = IloNumArray3(env, data_S->getN_arcs());
			for(IloInt a = 0; a < data_S->getN_arcs(); a++){
			     x_Sol[a] = IloNumArray2(env, data_S->getN_od());
			     for(IloInt k = 0; k < data_S->getN_od(); k++)
				x_Sol[a][k] = IloNumArray(env, search_param->getMaster_sc());			      
			}
			theta_SOL = IloNumArray (env, n_sc + search_param->getMaster_sc());
			Sol_send  = new double[data_S->getN_arcs() + n_sc  + 2 + 1]; // y sol and obja val, theta sol, Mp solving time ... 1 is for obj of the global scenarios
			
			(cplex).setParam(IloCplex::EpGap, 0.005);
			   
		}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		void cardinality_cuts(Data_S* data_S, Search_Param* search_param, int n_sc)
		{ 
		      //fixing looping arcs
		      for (int a = 0; a < data_S->getN_arcs() - data_S->getN_od(); a++)
			if(data_S->getArcs(a,0) == data_S->getArcs(a,1))
			  model.add(IloRange(env, 0, y[a], 0));
			//______________Cardinality cuts based on the whole scenarios_________________
 			int Origin, Destination, Card_Cut_O, Card_Cut_D, Card_Count=0;// ; ; cardinality cut for the origins (out going arcs); cardinality cuts for the destination (in comming arcs)
 			double Origin_Flow, Destination_Flow, Max_Origin_Flow, Max_Destination_Flow;
			 	
			 //for destination nodes
 			for (int k = 0; k < data_S->getN_od(); k++)	//master_sc + ----> later when you further developed the code, look into all scenatios to find the most violated one, for now it is not necessary to check all scenarios	
 			{
 			    Max_Origin_Flow =0;
			    Max_Destination_Flow =0;
			    
 			    Origin 	   = data_S->getOd(k,0);		//= origin[k]-1;
  			    Destination    = data_S->getOd(k,1);	//= destination[k]-1;
 				 
 			    for(int ss = 0; ss < n_sc + search_param->getMaster_sc(); ss++)//here we find the maximum flow going out (into) from origin (destination) among all scenarios nodes
 			    {	
 			          Origin_Flow 	   = 0; 
  				  Destination_Flow = 0;
 			    
 				  for (int kk = 0; kk < data_S->getN_od(); kk++)//this is to find all flow that originates - destinates from a node since one node can be the o or d fro several k
 				  {
 				    if (data_S->getOd(kk,0) == Origin)				    
 				      Origin_Flow += data_S->getD(kk,ss);
 				    if (data_S->getOd(kk,1) == Destination)
 				      Destination_Flow += data_S->getD(kk,ss); 
 				  }
 				  if (Origin_Flow > Max_Origin_Flow)
				    Max_Origin_Flow = Origin_Flow;
				  if(Destination_Flow > Max_Destination_Flow)
				    Max_Destination_Flow = Destination_Flow;
			    }
 				    
 				  IloExpr expr(env);
 				  IloExpr expr1(env);
 				  for (int a = 0; a < data_S->getN_arcs() - data_S->getN_od(); a++)//I am only writing this inequality on the regular arcs because aux arcs are already open (they are the last arcs added to the model in the numbering sense)
 				  {				    
 				      if ( data_S->getArcs(a,0) == Origin ) //for the originations 
 					expr += y[a];
 				      if ( data_S->getArcs(a,1) == Destination )//for the destinations
 					expr1 += y[a];
 				  }
 				   
 				  Card_Cut_O = ceil(Max_Origin_Flow/data_S->getU(1));
 				  Card_Cut_D = ceil(Max_Destination_Flow/data_S->getU(1));
//   				  cout << endl << "Origin_Flow/u[1] is: " <<  Max_Origin_Flow/u[1] << "  while Card_Cut is " << Card_Cut_O << endl;
 				  if( Card_Cut_O - (Max_Destination_Flow/data_S->getU(1)) > 0.001 )
				  {
 				    (con).add(IloRange(env, Card_Cut_O, expr, IloInfinity));
				    Card_Count++;
				  }
 			      	  if( Card_Cut_D - (Max_Destination_Flow/data_S->getU(1)) > 0.001 )
				  {
 				    (con).add(IloRange(env, Card_Cut_D, expr1, IloInfinity));
				    Card_Count++;
				  }
 			      				  
 				expr.end();
 				expr1.end(); 			    
 			}
 			 			
 			cout << endl << "I have added " << Card_Count << " Cardinality cuts" <<endl; 
			
 			(model).add((con));			
		}				
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////		
		//this function convert the first-stage variable to integer 
		void y_fixing(PoolsandMangers* Managers, Data_S* data_S, int i)
		{
		      int sol_id = Managers->getBest_UB_sol_id(i);
		      for(int a =0; a <data_S->getN_arcs(); a++)
		      {
			if(Managers->getY_Sol_Pool(sol_id, a) < 0.4 )
			   (*y_con).add(IloRange(env, 0, y[a], 0));//y[a].setUB(0);
			if(Managers->getY_Sol_Pool(sol_id, a) > 0.6 )
			  (*y_con).add(IloRange(env, 1, y[a], 1));
		      }	
		      model.add(*y_con);
		}
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////		
		double* heuristic_solve(PoolsandMangers* Managers, Data_S* data_S, Search_Param* search_param, int n_sc)
		{                             			
			double m_objval, cplex_start_time;
//  			cplex_start_time = (cplex).getTime() ;
			(cplex).setOut(env.getNullStream());//to not print the solution ssteps on the my file
			if( !(cplex).solve() ) 
			{			  
			  env.error() << "Failed to optimize Restircted heuristic problem\n";
// 			  return 9999999999999;
			}
			else 
			{
// 			  Heuristic_Sol =(cplex).getObjValue();	
			    m_objval = (cplex).getObjValue();
			    int ss;				    
			    (cplex).getValues(y_SOL, y);
			    double Recourse_Obj = 0;
			    for(int a = 0; a < data_S->getN_arcs(); a++)  			      
			      for(int k = 0; k < data_S->getN_od(); k++){
				 (cplex).getValues( x_Sol[a][k], x[a][k]);//the second stage objective for the global scenarios 
				for(int s = 0; s < search_param->getMaster_sc(); s++){			     
				      ss = Managers->getglobal_scenario_id(s);//ss = global_scenario_id[s];
				      Recourse_Obj +=  data_S->getP(ss) * data_S->getC(a) * x_Sol[a][k][s];
				  }				  
				}			    
			    (cplex).getValues(theta_SOL , theta);		    
			    
			//----------put the information we obtain from solving the master problem into an array
			//to be send to the processor-0 in the main			
			  for(int a = 0; a < data_S->getN_arcs(); a++)
			      Sol_send[a] = y_SOL[a];			
			  Sol_send[data_S->getN_arcs()]= m_objval;			
			  for(int s = 0; s < n_sc ; s++)
			      Sol_send[data_S->getN_arcs() + 1 + s] = theta_SOL[s];			
			  Sol_send[data_S->getN_arcs() + n_sc + 1] = 1;	
			  Sol_send[data_S->getN_arcs() + n_sc + 1 + 1] = Recourse_Obj;		
			}		
			
  			cout << "  obj is  " << m_objval << endl;
			model.remove(*y_con);
			// delete y_con;
			y_con->clear();
			(*y_con).end();
			y_con		=  new IloRangeArray(env);
			return  Sol_send; 			
		}	
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		~MasterHeuristic()
		{
		}


};
#endif
