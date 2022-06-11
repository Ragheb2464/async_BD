// #ifndef CLASS_CB_H
// #define CLASS_CB_H
// 
// 
// #include <ilcplex/ilocplex.h>
// ILOSTLBEGIN
// #define Comm COMM_WORLD

	IloNumArray y_SOL;
	IloNumArray theta_SOL;
	IloNumArray3 x_Sol;
	IloNumVarArray setSolVar;
	IloNumArray setSolVal;
	IloConstraintArray branchConst;
					  
	
	class_CI *coverCuts;
	class_MCI  *cardinalityCuts;
	OptRound  *optRound;
		
	int n_workers, n_sc;			
	double *Sol_send, *GlobalLowerBound;
	bool stopping_ceriteria;
	IloInt LazyNodeID=-1, UserNodeID=-1, UserNodeID1=-1, UserNodeID2=-1;
	vector <vector <double> > node_prune_Pool;
	IloNumArray reducedCosts;
	double glob_LB_in_tree=0;
	int *arry_y_fixed;
	int total_fix_Y;
	IloRangeArray interPhaseCuts;//to keep the benders cuts we generate during the intermidiary phase
	vector < vector <double> >  stored_nodes, node_avail_to_select;
	bool rootBranch =true;
ILOLAZYCONSTRAINTCALLBACK7(BendersLazyCallback, IloNumVarArray3, x, IloNumVarArray, y, IloNumVarArray,  theta, PoolsandMangers*, Managers, MasterSolChoice*, Sol_Choice, Data_S*, data_S, Search_Param*, search_param)		
{		 		
		if( LazyNodeID != getNnodes() ){
		    LazyNodeID = getNnodes();
			cout << " **BendersLazyCallback-" /* << getNodeId() */ << endl;
		    IloEnv env = getEnv();					    
  			// cout << endl << "NodeID: " << 	LazyNodeID <<endl;	    
		    //get the values we need(plus those for IloIncumbentCallBack)
		    double Recourse_Obj=0;
			int solId =1e-2;//the solution id which we add the combinatorial cut for										    	
			int t = data_S->getN_arcs() + n_sc;
			getValues(setSolVal, setSolVar);
		    for(IloInt a = 0; a < data_S->getN_arcs(); a++)	
			  for(IloInt k = 0; k < data_S->getN_od(); k++)
				 for(IloInt s = 0; s < search_param->getMaster_sc(); s++)
					Recourse_Obj +=  data_S->getP(Managers->getglobal_scenario_id(s)) * data_S->getC(a) * setSolVal[t++];//x_Sol[a][k][s];						  
		    //send the information which is required by workers		
			t=0;			
		    for(IloInt a = 0; a < data_S->getN_arcs(); a++)
				Sol_send[a] = setSolVal[t++];//y_SOL[a];					    
		    Sol_send[data_S->getN_arcs()]= getObjValue();					    
		    for(IloInt s = 0; s < n_sc ; s++)
				Sol_send[data_S->getN_arcs() + 1 + s] = setSolVal[t++];//theta_SOL[s];					    
		    Sol_send[data_S->getN_arcs() + n_sc + 1] = 0;	
		    Sol_send[data_S->getN_arcs() + n_sc + 1 + 1] = Recourse_Obj;
														  
									  
			
			Managers->Sol_Pool_Manager(true, search_param,data_S, search_param->getMaster_sc(), n_sc, n_workers, Sol_send,  GlobalLowerBound[0], GlobalLowerBound[1]) ;
			GlobalLowerBound[1] = GlobalLowerBound[0];
			Managers->Sol_Pool_controler(true, n_workers, n_sc , data_S,  search_param);	//sends solutions to the workers "asynchronously" based on the given strategy ... he wont send if there is nothing to send					
			if(search_param->getCreate_art_SPs()){
					Managers->art_cut_pool_manager(n_sc, n_workers, GlobalLowerBound[0], data_S, search_param);			 				            				      							  						  
					int sc_id;//the last one is to check we are not adding the same CBC cut several times
					for(int i = 0; i < Managers->getv_cp_artSize(); i++){
								sc_id  = Managers->getv_cp_art(i, 1);//(*v_cp)[i][1];	//this indicates which scenario this cut belongs to
								IloExpr expr(env);				    
								for(int a = 0; a < data_S->getN_arcs(); a++)
									expr += Managers->getv_cp_art(i, a+3)*y[a];									
								expr += Managers->getv_cp_art(i, 2);
								int w_id = sc_id-n_sc;
								double prob_sum=0.0;
								for(int s=0; s< Managers->get_worker_sc_size(w_id); s++)
									prob_sum += data_S->getP(Managers->get_worker_sc_ids(w_id, s));
								for(int s=0; s< Managers->get_worker_sc_size(w_id); s++)
										expr -= (data_S->getP(Managers->get_worker_sc_ids(w_id, s))/prob_sum) * theta[Managers->get_worker_sc_ids(w_id, s)];
								add(expr <= 0, IloCplex::UseCutPurge).end();					      
								expr.end();
					}							  	 			 			
					Managers->Cleanv_cp_art();
					// Managers->Cleanv_cp_copy();
			}
			Managers->cut_pool_manager( false, stopping_ceriteria, n_sc, n_workers,  GlobalLowerBound[0], data_S,  search_param); //v_sc_time,v_w_time, return number of cuts that we can receive			 				            				      											
			//add regular cuts
			// double viol=0.0;
			
			
			
			
			
			
			
			
			bool addCC=false;//to see if we need to add combinatorial cut when a subproblem is infeasible
			for(int j = 0; j < search_param->getNumAggrCluster(); j++){
				IloExpr expr(env), expr1(env);
				for(int k = 0; k < Managers->getClusterSize(j); k++)
					for(int i = 0; i < Managers->getv_cpSize(); i++)
						if(Managers->getv_cp(i, 0)== Managers->getY_Sol_Poolsize()-1 && Managers->getv_cp(i, 1) == Managers->getClusterIDs(j, k)  ){
							// sol_id = Managers->getv_cp(i, 0);//(*v_cp)[i][0];	//solution ID
							// cout <<  j << " , " << Managers->getv_cp(i, 1) << " - " << Managers->getClusterIDs(j, k) << endl;
							// Num_added_cuts++;
							int sc_id  = Managers->getv_cp(i, 1);//(*v_cp)[i][1];	//this indicates which scenario this cut belongs to
							for(int a = 0; a < data_S->getN_arcs(); a++)
									expr += Managers->getv_cp(i, a+3)*y[a];
							expr += Managers->getv_cp(i, 2);							
							expr -= theta[sc_id];
							if(!addCC && Managers->getv_cp(i, 3+data_S->getN_arcs())<0.5){
								addCC=true;
								solId=Managers->getv_cp(i, 0);
							}
						}
						else if(Managers->getv_cp(i, 0) != Managers->getY_Sol_Poolsize()-1 && Managers->getv_cp(i, 1) == Managers->getClusterIDs(j, k) ){
							int sc_id  = Managers->getv_cp(i, 1);//(*v_cp)[i][1];	//this indicates which scenario this cut belongs to
							for(int a = 0; a < data_S->getN_arcs(); a++)
									expr1 += Managers->getv_cp(i, a+3)*y[a];
							expr1 += Managers->getv_cp(i, 2);							
							expr1 -= theta[sc_id];
							if(!addCC && Managers->getv_cp(i, 3+data_S->getN_arcs())<0.5){
								addCC=true;
								solId=Managers->getv_cp(i, 0);
							}
						}
				add(expr <= 0).end();//(*con_optcut).add(IloRange( env, -IloInfinity, expr, 0));//, con_name));					      
				add(expr1 <= 0, IloCplex::UseCutPurge ).end();//(*con_optcut).add(IloRange( env, -IloInfinity, expr, 0));//, con_name));					      
				expr.end(); 
				expr1.end(); 
			}
			if(addCC && search_param->getAddCombinatorialCuts()){//if we want to add CBC, and the solution has been infeasible, and it has not been added previousely
							cout << "adding CBC" << endl;
							// solId = Managers->getv_cp(i, 0);
							IloExpr expr11(env);
							for(int a = 0; a < data_S->getN_arcs() - data_S->getN_od(); a++)
								if(Managers->getY_Sol_Pool(solId, a) < 1e-3)
									expr11 += y[a];
							add( expr11 >= 1 , IloCplex::UseCutPurge).end();//, con_name));	
							// (interPhaseCuts).add(IloRange(env, 1, expr, IloInfinity));
							expr11.end();
			}
			
				//free memory //clean up the cut pool			 			 															
				Managers->Cleanv_cp();	
				// Managers->Cleanv_cp_copy();	
		}
		return;
} // END BendersLazyCallback				

//for the frac solutions 
ILOUSERCUTCALLBACK7(BendersUserCallback, IloNumVarArray3, x, IloNumVarArray, y, IloNumVarArray,  theta, PoolsandMangers*, Managers, MasterSolChoice*, Sol_Choice, Data_S*, data_S, Search_Param*, search_param)		
{
	
	if(false && UserNodeID2 != getNnodes() && UserNodeID2 <= 100){
			UserNodeID2 = getNnodes();
			//cout << " **BendersLazyCallback-" /* << getNodeId() */ << endl;
		    IloEnv env = getEnv();					    
  			// cout << endl << "NodeID: " << 	LazyNodeID <<endl;	    
		    //get the values we need(plus those for IloIncumbentCallBack)
		    double Recourse_Obj=0;
			int solId =1e-2;//the solution id which we add the combinatorial cut for										    	
			int t = data_S->getN_arcs() + n_sc;
			getValues(setSolVal, setSolVar);
		    for(IloInt a = 0; a < data_S->getN_arcs(); a++)	
			  for(IloInt k = 0; k < data_S->getN_od(); k++)
				 for(IloInt s = 0; s < search_param->getMaster_sc(); s++)
					Recourse_Obj +=  data_S->getP(Managers->getglobal_scenario_id(s)) * data_S->getC(a) * setSolVal[t++];//x_Sol[a][k][s];						  
		    //send the information which is required by workers		
			t=0;			
		    
			
			
			if(optRound->FixVars(setSolVal, data_S))
				if(!optRound->solve(setSolVal, data_S))//if couldnt find a feasible sol go out
						return;
			
			
			
			for(IloInt a = 0; a < data_S->getN_arcs(); a++){
				setSolVal[a]=round(setSolVal[a]+0.1);				
				Sol_send[a] = setSolVal[t++];//y_SOL[a];
			}				
		    Sol_send[data_S->getN_arcs()]= getObjValue();					    
		    for(IloInt s = 0; s < n_sc ; s++)
				Sol_send[data_S->getN_arcs() + 1 + s] = setSolVal[t++];//theta_SOL[s];					    
		    Sol_send[data_S->getN_arcs() + n_sc + 1] = 0;	
		    Sol_send[data_S->getN_arcs() + n_sc + 1 + 1] = Recourse_Obj;
														  
									  
			
			Managers->Sol_Pool_Manager(true, search_param,data_S, search_param->getMaster_sc(), n_sc, n_workers, Sol_send,  GlobalLowerBound[0], GlobalLowerBound[1]) ;
			GlobalLowerBound[1] = GlobalLowerBound[0];
			Managers->Sol_Pool_controler(true, n_workers, n_sc , data_S,  search_param);	//sends solutions to the workers "asynchronously" based on the given strategy ... he wont send if there is nothing to send					
			
			//Managers->cut_pool_manager( false, stopping_ceriteria, n_sc, n_workers,  GlobalLowerBound[0], data_S,  search_param); //v_sc_time,v_w_time, return number of cuts that we can receive			 				            				      											
			
	}
	
	
	
	
	
	
	
	
	
	
	
	
	if( UserNodeID != getNnodes() && Managers->auxcut_pool_manager(false, n_sc, n_workers,  GlobalLowerBound[0], data_S,  search_param) > 0) {
		UserNodeID = getNnodes();
		IloEnv env = getEnv();
		int solId =1e-2;      
		bool addCC=false;//to see if we need to add combinatorial cut when a subproblem is infeasible
			for(int j = 0; j < search_param->getNumAggrCluster(); j++){
				IloExpr expr(env);
				for(int k = 0; k < Managers->getClusterSize(j); k++)
					for(int i = 0; i < Managers->getv_cpSize(); i++)
						if(Managers->getv_cp(i, 1) == Managers->getClusterIDs(j, k)  ){
							int sc_id  = Managers->getv_cp(i, 1);//(*v_cp)[i][1];	//this indicates which scenario this cut belongs to
							for(int a = 0; a < data_S->getN_arcs(); a++)
									expr += Managers->getv_cp(i, a+3)*y[a];
							expr += Managers->getv_cp(i, 2);							
							expr -= theta[sc_id];
							if(!addCC && Managers->getv_cp(i, 3+data_S->getN_arcs())<0.5){
								addCC=true;
								solId=Managers->getv_cp(i, 0);
							}
						}				
				add(expr <= 0, IloCplex::UseCutPurge ).end();//(*con_optcut).add(IloRange( env, -IloInfinity, expr, 0));//, con_name));					      
				expr.end(); 
			}			  		 			 															
		    Managers->Cleanv_cp();	
	}
	





	
	return;
} // END BendersUserCallback 

ILOINCUMBENTCALLBACK1(ControlIncumSol, PoolsandMangers*, Managers)//it will be called whenever a new potential incumbent is found
{
	      cout << "**IncumbentCallBack" << endl;
	      //we want to have the upper bound only from benders which is on the original formulation becasue the upper bound 
	      //it calculates inside the tree might not be a right one due to asynchronous behaviour of the algorithm which 
	      //does not wait fro all the cuts
	      //so we reject any new incumbent if its is not the one we have give it to him in ILOHEURISTICCALLBACK0 or 
	      //the incumbent solution that finds and number of waiting cuts are zero
		      
		      
	      /* if(getSolutionSource()	== UserSolution)
			cout << "UserSolution" << endl;
	      else
			if(getSolutionSource() == NodeSolution)
			  cout << "NodeSolution" << endl;
			else
			  if(getSolutionSource() == HeuristicSolution)
				cout << "HeuristicSolution" << endl; */
	    	    
		 if(getSolutionSource() != UserSolution /*&& Managers->getPotential_Cuts() !=0*/)//do not reject the sol if we dont have any more cuts to receive
			reject();	      
}//END ILOINCUMBENTCALLBACK

ILOHEURISTICCALLBACK6(MyHeuristic, IloNumVarArray3, x, IloNumVarArray, y, IloNumVarArray,  theta, PoolsandMangers*, Managers, Data_S*, data_S, Search_Param*, search_param)
{
//  	      cout <<  "**IloHeuristicCallBack" << endl;  
         
	    if(getIncumbentObjValue()-0.001 > Managers->getGlobalUpperBound(0) /*|| hasIncumbent()*/)//if we got the upper bound updated
		{
			if(Managers->getGlobalUpperBound(1) <= Managers->getGlobalUpperBound(0) ){
				cout << "\n Abort because the feasible solution has been rejected\n";
				//exit(-1);
				//abort();
			}
			  // IloEnv env = getEnv();
			  Managers->setGlobalUpperBound(1, Managers->getGlobalUpperBound(0));//GlobalUpperBound[1] = GlobalUpperBound[0];
			  cout << endl << "*********new upper bound with obj " <<  Managers->getGlobalUpperBound(0) << " is found " << endl;
			  //fill the vector of new incumbent sol
			  /* IloInt t=0;
			  for(IloInt a=0; a<data_S->getN_arcs(); a++){
			      setSolVal[t] = best_y[a];
			      t++;
			  }
			  for(IloInt s=0; s<n_sc;s++){
			      setSolVal[t] = best_theta[s] ;
			      t++;
			  }
			  for(IloInt a=0; a<data_S->getN_arcs(); a++)
			    for(IloInt k=0; k< data_S->getN_od(); k++)
			      for(IloInt s=0; s<search_param->getMaster_sc();s++){
					  setSolVal[t] = best_x[a][k][s];
					  t++;
			      } */
			 //propose the new incumbent sol
			  setSolution(setSolVar, setSolVal, Managers->getGlobalUpperBound(0) - 1); 
		}
}		  

ILOSOLVECALLBACK1(mySolve, IloNumVarArray, y)
{
  cout << endl << "MySolve" << endl;
// reducedCostFixing();
}

ILOCONTINUOUSCALLBACK1(myContinous, IloNumVarArray, y)
{  
}

/*ILOMIPINFOCALLBACK0(loggingCallback)
{
	
	if ( hasIncumbent() )
	{
    // 	  reject();
	      cout << "there is new incumbent" << endl;
	      getIncumbentValues(y_Incum, y);
	      
	    for(int a =0; a <n_arcs; a++)//reject it if its not what we have gave it to him (we gave him a solution where every arcs is open but its value is the benders upper bound)
	    if(y_Incum[a] >0 )
	      cout << a << " -> " << y_Incum[a] << "  " ;
	    cout << endl;
	    for(int s=0; s<n_sc; s++)
	      cout << s << " -> " << getIncumbentValue(theta[s])<< "  " ;
	    cout << endl;
	    
	    
	    
	    cout << "the cutoff value is: " << getCutoff() << endl ;
	}
	
//   	(cplex).setParam(IloCplex::CutUp, GlobalUpperBound[0]);  						
	
	

}*/	    

ILONODECALLBACK0(MySelect) 
{
//    cout << "**NodeSelectionCallBack" << endl;
//    IloInt remainingNodes = getNremainingNodes();
   IloInt /*maxdepth = 1e5,*/ Sel_Nod =-1/*, NodeObjHigh =-1*/, NodeObjLow =1e10;
//    IloInt depth;
  
   /*
   //select a node at the top of the tree
   for (IloInt i = 0; i < remainingNodes; i++) {
     cout << endl << "Active node with obj: " << getObjValue(i) << endl;
       depth = getDepth(i);
      if(depth < maxdepth){
	maxdepth =depth;
	Sel_Nod =i;
      }      
   }*/
   
   
   //based on the highest obj value
//    for (IloInt i = 0; i < remainingNodes; i++)
//      if( getObjValue(i) > NodeObjHigh){
//        NodeObjHigh = getObjValue(i);
//        Sel_Nod =i;
//      }
//    

   /*   for (IloInt i = 0; i < node_avail_to_select.size(); i++)
	if(node_avail_to_select[i][1] < 1)
	  for (IloInt j = i+1; j < node_avail_to_select.size(); j++)
	    if( fabs( node_avail_to_select[j][1] - node_avail_to_select[i][1]) <1)
	      NodeObjLow = node_avail_to_select[i][2] = 10;//this is not our periority to invistigate 


   //based on the lowest obj value
//    cout << endl << "Active node with obj: " << getObjValue(i) << endl;
    for (IloInt i = 0; i < node_avail_to_select.size(); i++)
       if(node_avail_to_select[i][2] < 1)      		//if node has been examined 
	if(  node_avail_to_select[i][1] < NodeObjLow){
	  NodeObjLow = node_avail_to_select[i][1];	//node obj
	  Sel_Nod =node_avail_to_select[i][0];		//node id
	}
    
  */
   
   if(Sel_Nod >= 0){
     selectNode(Sel_Nod);
     node_avail_to_select[Sel_Nod][2] =10;//
     cout << endl << "Selected node is: " << Sel_Nod << " with obj: " << node_avail_to_select[Sel_Nod][1] << endl;
   }
}
				
ILOBRANCHCALLBACK3(MyBranch, IloNumVarArray, y, PoolsandMangers*, Managers, Data_S*, data_S) 
{
	// cout << endl << "Branch" << endl;
  if(rootBranch){
	  IloEnv env = getEnv();
	  rootBranch=false;
	  //fix the variables
		int counter =0;
		int oneFixSize=0;
		// getValues(y_SOL,y);
		IloExpr expr0(env), expr1(env), expr01(env);//if the fixation yielded infeasible problem, then at least one of those fixed to zero has to be opened
		for(int a =0; a < data_S->getN_arcs() /* - data_S->getN_od() */; a++){
			if(y_SOL[a] <= 0.1){
			   expr0  += y[a];
			   expr01 += y[a];
			   counter++;
			}
			if(y_SOL[a] >= 0.9){
			   expr1  += y[a];
			   oneFixSize++;
			   expr01 += 1 - y[a];
			   counter++;
			}
		}
	   cout << "\nN of fixed variable in the left branch: " << counter << endl;
	  (branchConst).add(IloRange(env , 0, expr0, 0));	
	  (branchConst).add(IloRange(env , oneFixSize, expr1, oneFixSize));	
	  makeBranch(branchConst, 0, 0);
	  (branchConst).endElements();
	  (branchConst).add(IloRange(env , 1, expr01, IloInfinity));
	  makeBranch(branchConst, 1, 0);
	  (branchConst).endElements();
	  expr01.end();
	  expr0.end();
	  expr1.end();
	}
	// else
	// {			 
		/* getValues(y_SOL,y);
		for(int a=0; a< data_S->getN_arcs() - data_S->getN_od(); a++)
			if(y_SOL[a] > 0.001 && y_SOL[a] <0.9999)
				cout << endl << a << "  -: " << y_SOL[a]; */
		/* getValues(y_SOL,y);
		double sumInfes=0;
		for(int a=0; a< data_S->getN_arcs() - data_S->getN_od(); a++)
			if(y_SOL[a] > 0.1 && y_SOL[a] <0.9)
				sumInfes += y_SOL[a];
		if( fabs(sumInfes - round(sumInfes)) >= 0.1 ){
			IloEnv env = getEnv();
			IloExpr expr(env);
			for(int a=0; a< data_S->getN_arcs() - data_S->getN_od(); a++)
				if(y_SOL[a] > 0.1 && y_SOL[a] <0.9)
					expr += y[a];
			(branchConst).add(IloRange(env , 0, expr0, 0));
			makeBranch(expr, ceil(sumInfes), IloCplex::BranchUp,   getObjValue() );
			makeBranch(expr, floor(sumInfes), IloCplex::BranchDown, getObjValue() );
	
		} */
	// }
	  
	  
	  
	  
	  
	  
	  
	  
	  
	  

	  // 
	  /* IloNumArray y_temp(env, data_S->getN_arcs());
	  
	  vector <double> aux;
	  
	  aux.push_back(getObjValue());  
	  getValues(y_temp, y);
	  for(IloInt a=0; a<data_S->getN_arcs() ; a++)
	      aux.push_back(y_temp[a]);
	  aux.push_back((getNnodes()));
	  
	  cout << endl << getNnodes() << "  "<< aux[0] << "  ->";
	  for(IloInt a=0; a<data_S->getN_arcs() ; a++)
	    if(y_temp[a] > 1e-5)
	      cout << a << "- " << y_temp[a] << "; ";
	    cout << endl;
	    
	  y_temp.end(); */
	  
// 	  cout << endl << aux[0] << endl;

	  /* bool prunning =false, same_node_exist;
	  if( node_prune_Pool.size() >0)	  
	    for(IloInt i=0; i < node_prune_Pool.size(); i++)
	    {
// 	      cout << i << " " << node_prune_Pool[i][0] << endl;
		same_node_exist = true;
		if(node_prune_Pool[i][1+data_S->getN_arcs()] != aux[1+data_S->getN_arcs()])
		  if( fabs(node_prune_Pool[i][0] - aux[0]) >= 0.1 )
		  {
		    for(IloInt a=0; a <data_S->getN_arcs() ; a++)
		      if(fabs(node_prune_Pool[i][a+1] - aux[a+1]) >= 0.1)
			same_node_exist = false;
		  }
		  else
		    same_node_exist = false;
		  if(same_node_exist)
		  {
		    prunning = true;
		    cout << endl << "PrUUUUUUnnnnningg " << aux[1+data_S->getN_arcs()] << " since it is the same as node-" << i << endl;
		    break;
		  }
	    } */
	  
	/*if(prunning)
	  prune();
	else
	{	  
	  node_prune_Pool.push_back(aux);
	  aux.clear();
	}
	*/
	
	
	
	
	
	
	
	
//     cout << endl << "+++++++Info of node " << getNodeId()/*getNodeNumber(getNodeId(i))*/ << endl;
//     cout << "Obj: " << getObjValue() << endl;
//       for(IloInt a=0; a<data_S->getN_arcs() - data_S->getN_od(); a++)
//         if( getValue(y[a]) > 1e-5)
// 	  cout << a << "--" << getValue(y[a]) << ";  ";
// 	cout << endl;
//     
//   }
  
  /* if ( getBranchType() != BranchOnVariable )
      return;
    cout << endl << "**VarSelecCallBack" ;
   // Branch on var with largest objective coefficient
   // among those with largest infeasibility

   IloNumArray x;
   IloNumArray obj;
   IntegerFeasibilityArray feas;

   try {
      x    = IloNumArray(getEnv());
      obj  = IloNumArray(getEnv());
      feas = IntegerFeasibilityArray(getEnv());
      getValues(x, y);
      getObjCoefs(obj, y);
      getFeasibilities(feas, y);

      IloInt bestj  = -1;
      IloNum maxinf = 0.0;
      IloNum maxobj = 0.0;
      IloInt cols = y.getSize();
      for (IloInt j = 0; j < cols; j++) //based on max infeasibility
      {
         if ( feas[j] == Infeasible ) {
            IloNum xj_inf = x[j] - IloFloor (x[j]);
            if ( xj_inf > 0.5 )
               xj_inf = 1.0 - xj_inf;
            if ( xj_inf >= maxinf                  &&           (xj_inf > maxinf || IloAbs (obj[j]) >= maxobj)  ) {
               bestj  = j;
               maxinf = xj_inf;
               maxobj = IloAbs (obj[j]);
            }
         }
      }
//       IloNum max_diff=0;
//       for (IloInt j = 0; j < cols; j++)//based max difference between getUpPseudoCost and getDownPseudoCost
// 	if ( feas[j] == Infeasible )
// 	{
// 	  IloNum diff = fabs(getUpPseudoCost(y[j]) - getDownPseudoCost(y[j]));
// 	  if(diff > max_diff)
// 	  {
// 	    max_diff = diff;
// 	    bestj =j;
// 	  }
// 	}
      if ( bestj >= 0 ) {
	 cout << endl << "branching on a variable with value: " <<  x[bestj] <<endl;
         makeBranch(y[bestj], 1, IloCplex::BranchUp,   getObjValue() + getUpPseudoCost(y[bestj]));
         makeBranch(y[bestj], 0, IloCplex::BranchDown, getObjValue() + getDownPseudoCost(y[bestj]));
	 cout << endl <<getObjValue() << " Up: " << getObjValue() + getUpPseudoCost(y[bestj]) << " Down: " << getObjValue() + getDownPseudoCost(y[bestj]) << endl;
      }
   }
   catch (...) {
      x.end();
      obj.end();
      feas.end();
      throw;
   }
   x.end();
   obj.end();
   feas.end();*/
}



// #endif
