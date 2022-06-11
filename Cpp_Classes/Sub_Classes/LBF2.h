#ifndef CLASS_LB_INEQ_H
#define CLASS_LB_INEQ_H





ILOSTLBEGIN
#define Comm COMM_WORLD



class class_LB_Lifting
{ 
  private:    
		double *min_d; //for the vector of minimum demand of each commodity over all scenarios
		
		IloEnv env;
		IloModel model;
		IloCplex cplex;
		IloRangeArray con, con1; //this is for the main constraints (flow conservation and capacity) constraints				
		IloNumVarArray Dual_Capacity_Con;		
		IloNumVarArray2  Dual_Flow_Con, Dual_Strong_Ineq;
			
		IloNumVarArray var;
		IloNumArray val;
		IloNumArray X;		 		
		IloObjective objFun;
		
		IloNumArray y_SOL;
		IloNumArray2 x_SOL; 
		IloNumArray dual_SOL;
		double routing_cost;
		double max_LBF;
		double *scenario_ave;
		int *cluster_id;//scenario s is in which cluster c	
		vector < vector <double > > Cluster_Core;
			  
  public:        		
		class_LB_Lifting(int n_sc, MasterSolChoice *Sol_Choice, Data_S *data_S, Search_Param	*search_param)		
		{	
			model 		=  IloModel(env);
			cplex 		=  IloCplex(model);
			con 		=  IloRangeArray(env);
			con1 		=  IloRangeArray(env);			
			objFun		= IloMaximize(env);
			
			y_SOL = IloNumArray(env, data_S->getN_arcs());
			dual_SOL = IloNumArray(env, data_S->getN_od());
			x_SOL = IloNumArray2(env, data_S->getN_arcs());
			for(IloInt a = 0; a < data_S->getN_arcs(); a++)
			     x_SOL[a] = IloNumArray(env, data_S->getN_od());		      			
			
			n_sc = n_sc - search_param->getMaster_sc(); // number of scenarios are total scenarios minus global scenarios (those who are not projected)	 			
			min_demand( data_S, n_sc);
						
			
			//________________Variables definition_____________
			Dual_Flow_Con = IloNumVarArray2(env, data_S->getN_nodes());
			//			 IloNumVarArray2 Dual_Flow_Con(env, n_nodes);//if I define it in this way I need the "this" because I am defining a new one but similar
			for(int i = 0; i < data_S->getN_nodes(); i++)
			  Dual_Flow_Con[i] = IloNumVarArray(env, data_S->getN_od(), -IloInfinity, IloInfinity); // for the flow conversation constraints 			
			//			  this->Dual_Flow_Con = Dual_Flow_Con;
			
 			for(int i = 0; i < data_S->getN_nodes(); i++)
 			    (model).add(Dual_Flow_Con[i]);
// 			  			
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
			  expr += fabs(min_d[k]) * (Dual_Flow_Con[data_S->getOd(k,1)][k]);

			
			for (int a = 0; a < data_S->getN_arcs() ; a++)
			  if (data_S->getArcs(a,1) != data_S->getArcs(a,0))//if it is not looping arc
			  {
			    expr += -data_S->getU(a) * Dual_Capacity_Con[a];//the part associated with capacity constraints
			    for (int k = 0; k < data_S->getN_od(); k++)
			      expr += -min(data_S->getU(a), fabs(min_d[k])) * Dual_Strong_Ineq[a][k];
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
				    (con).add(IloRange(env, -IloInfinity, expr , data_S->getC(a)+ (data_S->getF(a)/ data_S->getU(a)) ) );
				    expr.end();				
			      }
			
			//the dual variable associated to the origin is equal to zero
			for(int k = 0; k < data_S->getN_od(); k++)
			  (con1).add(IloRange(env, 0, Dual_Flow_Con[data_S->getOd(k,0)][k] , 0));						
			
			//adding the strong inequalities to the sub-problem formulation
			if(false /*/!search_param->getStrong_Ineq()*/)//if it is false we remove the dual variables associate to strong inequalities (fix them at zero)
			{			  
			   for ( int a = 0; a < data_S->getN_arcs() ; a++)
				for ( int k = 0; k < data_S->getN_od(); k++)
				  (con1).add(IloRange(env, 0,  Dual_Strong_Ineq[a][k] ,0));			 
			}

		        (model).add((con));
   			(model).add((con1));						
			// clustering(data_S, search_param, n_sc);
			solveModel(data_S, n_sc);
		}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		void min_demand( Data_S * data_S, int n_sc)
		{
		      min_d = new double[ data_S->getN_od()];
		      double minDemand;
		      for(int k =0; k< data_S->getN_od(); k++){
				minDemand = 1e10;
				for(int s=0; s < n_sc; s++)
				  if(data_S->getD(k,s) <= minDemand)
					minDemand = data_S->getD(k,s);
				min_d[k] = minDemand;
		      }		  
		}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		double get_distance(Data_S * data_S, int cluter_id, int s_id)
		{		  
		    double dist=0;
			for(int k =0; k< data_S->getN_od(); k++)
				  dist += pow(Cluster_Core[cluter_id][k] - data_S->getD(k,s_id),2);
			return sqrt(dist);
		}		
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		 void clustering(Data_S * data_S, Search_Param  *search_param, int n_sc)
		{								
			cluster_id = new int[n_sc]();			
			initial_cluster(search_param, data_S,n_sc);
			double dist;
			int closeest_cluster;
			bool stop =false;
			
			while(!stop){
				stop = true;
				for( int c =0 ; c < search_param->getLBF_clusters(); c++)
					get_cluster_core(c,data_S,n_sc);
				for(int s= 0; s< n_sc; s++){
					dist = 1e10;
					for( int c =0 ; c < search_param->getLBF_clusters(); c++)
						if(get_distance(data_S, c, s) < dist){
							closeest_cluster = c;
							dist = get_distance(data_S, c, s);
						}																
					if(cluster_id[s] != closeest_cluster){//move
						cluster_id[s] = closeest_cluster;
						stop =false;
						break;
					}										
				}
			}
			
			print_cluter(search_param,n_sc);		
		} 
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		void initial_cluster(Search_Param  *search_param, Data_S * data_S, int n_sc)
		{
			vector < double> aux;
			int t=0;
			//assign scenarios to clusters
			while(t < n_sc){
				for( int c=0 ; c < search_param->getLBF_clusters(); c++){
					cluster_id[t] = c;
					t++;
				}
			}
			//calculate the center of each cluster
			for( int c=0 ; c < search_param->getLBF_clusters(); c++){
				aux.clear();
				for(int k =0; k< data_S->getN_od(); k++)
					aux.push_back(data_S->getD(k,c));
				Cluster_Core.push_back(aux);
			}
		}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		void get_cluster_core(int cluter_id, Data_S * data_S, int n_sc)
		{
			double *cluster_ave;
			int cluseterSize=0;
			cluster_ave = new double[data_S->getN_od()]();
			for( int s=0 ; s < n_sc; s++)
				if(cluster_id[s] == cluter_id)	{
					cluseterSize++;					
					for(int k =0; k< data_S->getN_od(); k++)
						cluster_ave[k] += data_S->getD(k,s);
				}					
			
			for(int k =0; k< data_S->getN_od(); k++)
				Cluster_Core[cluter_id][k]= cluster_ave[k]/cluseterSize;
			
			delete [] cluster_ave;					
		}

////////////////////////////////////////////////////////////////////////////////////////////////////		
			void print_cluter(Search_Param  *search_param, int n_sc)
			{//print clusters				
				for( int c=0 ; c < search_param->getLBF_clusters(); c++){
// 					cout << endl << "scenarios in the " << c << "'th cluster: ";
					for( int s=0 ; s < n_sc; s++)
						if(cluster_id[s] == c)
							cout << s << "; ";
				}
				cout << endl;
			}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	      void solveModel( Data_S * data_S, int n_sc)
	      {
		(cplex).solve();
		cout << endl << "the cost of artificial scenario: " << cplex.getObjValue() <<endl;
		for ( int a = 0; a < data_S->getN_arcs() ; a++)
		  if(cplex.getValue(Dual_Capacity_Con[a]))
		    cout << cplex.getValue(Dual_Capacity_Con[a]) << endl;
		//get the dual variables Pi_D(k) associated to the flow conservation constraints
		for(int k = 0; k < data_S->getN_od(); k++)
		  dual_SOL[k] = (cplex).getValue( Dual_Flow_Con[data_S->getOd(k,1)][k]);
// 		env.out() << "dual values: "<< dual_SOL;
		
		//get the X variables
		X = IloNumArray(env,(con).getSize());
	    	(cplex).getDuals(X, con);
		int t =0;
		for ( int a = 0; a < data_S->getN_arcs() ; a++){
		    if (data_S->getArcs(a,1) != data_S->getArcs(a,0)){
		      for ( int k = 0; k < data_S->getN_od(); k++)				      
				x_SOL[a][k] = X[t++];//(cplex).getDual(con[t++]);
		    }
		    else{
		      for ( int k = 0; k < data_S->getN_od(); k++)				      
				x_SOL[a][k] = 0;//(cplex).getDual(con[t++]);
		    }
		}      
		X.end();
		
		//getting the routing cost
		routing_cost=0;
		for(IloInt a = 0; a < data_S->getN_arcs(); a++)
		   for(int k =0; k< data_S->getN_od(); k++)
		     routing_cost += data_S->getC(a) * x_SOL[a][k];
		   
		//getting the y solutions
		for(IloInt a = 0; a < data_S->getN_arcs(); a++){
		  double arcFlow=0;
		  for(int k =0; k< data_S->getN_od(); k++)
		    arcFlow += x_SOL[a][k];
		  y_SOL[a] = arcFlow /data_S->getU(a);
//  		  if( y_SOL[a] )
//  		    cout << endl << a << "   " <<  y_SOL[a];
		}
		
		//printing the approximated cost for each scenario and finding the max among scenarios
		   max_LBF=0;//maximum approximated objective
		   for(int s=0; s< n_sc; s++)
		   {
		     double costs=0;
		     costs += routing_cost;
		      for(int k=0; k< data_S->getN_od(); k++)
				 costs += (data_S->getD(k,s) - getMinD(k)) * getDual(k);
		      for(int a=0; a < data_S->getN_arcs() - data_S->getN_od(); a++)
				 costs += data_S->getF(a) * getY(a);
		      if(max_LBF <= costs)
			max_LBF= costs;
// 		     cout << endl << "approximated cost of scemario " << s << " is  " << costs <<endl;
		   }
	      }
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	     double getY(int a) {return y_SOL[a];}
	     double getDual(int k) {return dual_SOL[k];}
	     double getRotuingCost() {return routing_cost;}
	     double getMinD(int k) {return min_d[k];}
	     double getMax_LBF() {return max_LBF;}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

		~class_LB_Lifting()
		{
		}


};//WND class



#endif
