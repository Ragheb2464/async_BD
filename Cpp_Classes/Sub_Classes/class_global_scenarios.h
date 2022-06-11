#ifndef CLASS_GLOBAL_Scenario_H
#define CLASS_GLOBAL_Scenario_H


ILOSTLBEGIN
#define Comm COMM_WORLD
typedef IloArray<IloNumVarArray> IloNumVarArray2;
class GlobalScenarios
{ 

  private: 
	IloEnv env;
	IloModel model;
	IloCplex cplex;
	IloObjective objFun;
	IloRangeArray con; //this is for the main constraints (flow conservation and capacity) constraints				
	IloNumVarArray alpha, r, auxD;		
	IloNumVarArray2  b, bBar, e;
	//
	IloNumArray rVal, dVal, alphaVal;
	int K,A;
	double coeff;
	
  public:
	GlobalScenarios() 
	{	  	      
	}
/////////////////////////////
	void creatModel(Search_Param* search_param, int n_od, double  **ins_scenario, int n_sc)
	{
	      K = search_param->getMaster_sc();
	      A = search_param->getNumArtificialScenarios();
	      rVal	= IloNumArray(env,n_sc);
	      dVal	= IloNumArray(env,n_od);
	      alphaVal	= IloNumArray(env,n_sc);
		  coeff =1.0;
		  coeff /= n_sc-K;//0.0384615384;
		  //cout << endl <<1/62 << " ; K= " << K << "; " << "A= " << A << "; coeff= " << coeff << endl;
	      //
	      model		=  IloModel(env);
	      cplex		=  IloCplex(model);
	      objFun		= IloMaximize(env);
	      con		= IloRangeArray(env);
	      //var
	      alpha 	= IloNumVarArray(env, n_sc, 0, 1);
	      r	 	= IloNumVarArray(env, n_sc, 0, 1, ILOINT);
	      b 		= IloNumVarArray2(env, n_sc);
	      for(int s=0; s< n_sc; s++)
			b[s] 	= IloNumVarArray(env,  n_od, 0, 1, ILOINT);
	      bBar 		= IloNumVarArray2(env, n_sc);
	      for(int s=0; s< n_sc; s++)
			bBar[s] = IloNumVarArray(env,  n_od, 0, 1, ILOINT);
	      e 		= IloNumVarArray2(env, n_sc);
	      for(int s=0; s< n_sc; s++)
			e[s] 	= IloNumVarArray(env,  n_od, 0, IloInfinity);
	      auxD 		= IloNumVarArray(env, n_od, 0, IloInfinity);
	      //obj
	      IloExpr expr(env);
	      for(int s=0; s< n_sc; s++)
			for(int l=0; l< n_od; l++)
				expr += b[s][l];
	      objFun.setExpr(expr);
	      model.add(objFun);
	      expr.end();
	      //con
	      for(int s=0; s< n_sc; s++)
			for(int l=0; l< n_od; l++){
				IloExpr expr(env);
				for(int ss=0; ss< n_sc; ss++)
					if(ins_scenario[ss][l] >= ins_scenario[s][l])
						expr += r[ss];
				expr += bBar[s][l];
				expr -= b[s][l];
				con.add(IloRange(env, 0, expr, IloInfinity));
				expr.end();
			}
	      //
	      for(int l=0; l< n_od; l++){
			IloExpr expr(env);
			for(int s=0; s< n_sc; s++)
			  expr += ins_scenario[s][l] * alpha[s];
			expr -= auxD[l];
			con.add(IloRange(env, -1e-75, expr, 1e-75));
			expr.end();	    
	      }
	      //
	      for(int s=0; s< n_sc; s++)
			for(int l=0; l< n_od; l++){
				con.add(IloRange(env, ins_scenario[s][l], auxD[l] + e[s][l], IloInfinity));
			}
	      //
	      for(int s=0; s< n_sc; s++)
			for(int l=0; l< n_od; l++)
				if(ins_scenario[s][l]>1e-17)
					con.add(IloRange(env, -IloInfinity, bBar[s][l] + (e[s][l]/ins_scenario[s][l]), 1));
			
	      //
	      IloExpr expr1(env);
	      for(int s=0; s< n_sc; s++)
			expr1 += r[s];
	      con.add(IloRange(env, -IloInfinity, expr1, K));
	      expr1.end();
	      //
	      for(int s=0; s< n_sc; s++){
			con.add(IloRange(env, A*coeff, alpha[s] + A*coeff*r[s], A*coeff));
		}
	      //
	      /* for(int s=0; s< n_sc; s++)
			con.add(IloRange(env, A*coeff, alpha[s] + A*r[s], IloInfinity)); */
	      //
	      IloExpr expr11(env);
	      for(int s=0; s< n_sc; s++)
			expr11 += alpha[s];
	      con.add(IloRange(env, A, expr11, A));
	      expr11.end();
	      //
		  for(int s=0; s< n_sc; s++)
			  for(int l=0; l< n_od; l++)
				  con.add(IloRange(env, 0, e[s][l], ins_scenario[s][l]));
		/////////////////////////////////////////////////////////////////////
	      model.add(con);
	     (cplex).setParam(IloCplex::Threads, 1);
	     cplex.setOut(env.getNullStream());
	}
////////////////////////////
	void solveModel(double *artificialScenario, vector <int>& gScID, double *alphaWeights)
	{
// 	      if(A == 0 && K ==1)
	      if(!cplex.solve()){
			cout << "failed to solve the retension and creation problem" << endl;
			abort();
	      }
	      //		  
	      cplex.getValues(rVal, r);
		  // env.out() << rVal << endl;
	      cplex.getValues(dVal, auxD);
	      cplex.getValues(alphaVal, alpha);
	      //
	      for(int s=0; s< alphaVal.getSize(); s++)
			alphaWeights[s] = alphaVal[s];
	     for(int l=0; l<dVal.getSize(); l++)
			artificialScenario[l] = dVal[l];
	     gScID.clear();
	     for(int s=0; s<rVal.getSize(); s++)
			if(rVal[s] > 0)
			  gScID.push_back(s);
	}
////////////////////////////
	~GlobalScenarios()
	{
	}


};
#endif
