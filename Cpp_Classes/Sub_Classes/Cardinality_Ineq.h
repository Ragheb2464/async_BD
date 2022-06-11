#ifndef CLASS_MCI_INEQ_H
#define CLASS_MCI_INEQ_H

//#include <ilcplex/ilocplex.h>
ILOSTLBEGIN
#define Comm COMM_WORLD

/* bool PairCompare(const std::pair<double, int>& firstElem, const std::pair<double, int>& secondElem) 
{
  return firstElem.first < secondElem.first;
} */

class class_MCI
{ 
  private:  
			//for the lifting problems
			IloEnv env;
			IloModel model;
			IloCplex cplex;
			IloRangeArray con; 				
			IloNumVarArray  y;		 		
			IloObjective objFun;
			//required arraies and vectors
			int *openArcs, *closeArcs;//the set of open and close arcs
			double upEpsilon, downEpsilon; //epsilons for openCloseArc algorithm
			double violationThreshhold; // an inequality how should be violated to accept it
			vector< int> orderedSolution;//the solution ordered based on the current value and capacity
			int *coverSet; //the set of arcs in the minimal cover set
			double *coverCoeff;//cefficeint of the variable in the CI ineq
			int *liftedVar;//variables lifted
			int *cutSet; //the arcs in the cutset
			double maxCutFlow; //the flow over the cutset
			int *baseNodes;//the selected nodes to be the cutset
			vector <vector< double> > IneqPool;// to keep the generated violated cuts
			double bigU, bigD;//for the residual demand and capacity 
			int *capBasedOrder;//to order arcs based on their capacity
			double minCardinality;
			vector<int > liftList;//the list containg the arcs to be lifted
  public:        		
		class_MCI(Data_S *data_S, int n_sc)		
		{	
			//cplex objects
			model 		=  IloModel(env);
			cplex 		=  IloCplex(model);
			con 		=  IloRangeArray(env);			
			objFun		=  IloMinimize(env);
			//selection variable
			y = IloNumVarArray(env, data_S->getN_arcs(), 0, 1, ILOINT);
			(model).add(y);			
			//build objFun
			IloExpr expr(env);
			for(int a =0; a < data_S->getN_arcs() - data_S->getN_od(); a++)
				expr += 0.0 * y[a];
			objFun.setExpr(expr);			
			(model).add(objFun);
			expr.end();
			//add some shitty constraint for the sake of safety 
			(con).add(IloRange(env, 0, y[0] , 1));
			(model).add(con);
			//cplex parameters
			(cplex).setParam(IloCplex::MIPDisplay, 0);			 					
			(cplex).setParam(IloCplex::EpGap, 0.003);
			(cplex).setParam(IloCplex::Param::Threads, 1); 
 			//initialization of the required arraies and vectors 
			openArcs  =new int[data_S->getN_arcs()]();
			closeArcs =new int[data_S->getN_arcs()]();
			upEpsilon = 0.4999;
			downEpsilon = 0.4999;
			violationThreshhold= 1.0e-3;
			coverSet = new int[data_S->getN_arcs()]();
			coverCoeff = new double[data_S->getN_arcs()](); 
			liftedVar = new int[data_S->getN_arcs()]();
			cutSet = new int[data_S->getN_arcs()]();
			baseNodes = new int[2];//
			capBasedOrder = new int[data_S->getN_arcs()]();
			SortingCapacities(data_S);//order design variables based on their capacity in non-increasing order
		}		
//////////////////////////////////////////////////////////////////////////////////////////////////
		void MainLoop(Data_S *data_S, int cutSetCardinality, IloNumArray &y_SOL, int n_sc)
		{
			bool outArc, inArc;
			//clean the cut pool
			IneqPool.clear();
			//for the single note cutset 
			if(cutSetCardinality == 1){
				for(int i=0; i< data_S->getN_nodes(); i++){	
					if( IsOriginNode(data_S, i)){
						baseNodes[0] = i;
						baseNodes[1] = i;//our cutset nodes
						//cutset
						outArc = true;
						inArc  = false;						
						CutSetGenerator(data_S, outArc, inArc);//find the arcs across cut set (going out from cutset nodes)
						CutSetMaxFlow(data_S, outArc, inArc, n_sc);//find the flow across cutset
						OpenCloseArcs(data_S, y_SOL);//find the set of open and close arcs
						if(FindMinCardinality(data_S, y_SOL)){//lift violated ineq
							SortingSol(data_S, y_SOL);		
							LiftingProcedure(data_S);
							ViolatedCoverIneqPool(data_S, y_SOL);
						}						
					}
					if(IsDestNode(data_S, i) ){
						baseNodes[0] = i;
						baseNodes[1] = i;//our cutset nodes
						//cutset
						outArc = false;
						inArc  = true;						
						CutSetGenerator(data_S, outArc, inArc);//find the arcs across cut set (going out from cutset nodes)
						CutSetMaxFlow(data_S, outArc, inArc, n_sc);//find the flow across cutset
						OpenCloseArcs(data_S, y_SOL);//find the set of open and close arcs
						if(FindMinCardinality(data_S, y_SOL)){
							SortingSol(data_S, y_SOL);		
							LiftingProcedure(data_S);
							ViolatedCoverIneqPool(data_S, y_SOL);
						}
					}
				}				
				if(getNumCoverIneq()>0)
				  cout << "Number of violated cardinality inequalities with single cutset: " << getNumCoverIneq() << endl;
			}			
			//for 2 we should make sure they are not a OD pair 
			if(cutSetCardinality == 2){
				for(int i=0; i< data_S->getN_nodes(); i++)
					for(int j =i+1; j< data_S->getN_nodes(); j++)
						if(IJConnected(data_S,i,j,y_SOL) /*&& !ODPair(data_S,i,j)*/ ){//if there is at least one arc between them, they are a OD piar, at least one is origin or destination							
							baseNodes[0] = i;
							baseNodes[1] = j;
							/*if( IsOriginNode(data_S, i) || IsOriginNode(data_S, j))*/{//if at least one in origin node (to have non zero outgoing flow)
								//cutset
								outArc = true;
								inArc  = false;						
								CutSetGenerator(data_S, outArc, inArc);//find the arcs across cut set (going out from cutset nodes)
								CutSetMaxFlow(data_S, outArc, inArc, n_sc);//find the flow across cutset
								if(maxCutFlow>0){
								  OpenCloseArcs(data_S, y_SOL);//find the set of open and close arcs
								  if(FindMinCardinality(data_S, y_SOL)){
									  SortingSol(data_S, y_SOL);		
									  LiftingProcedure(data_S);
									  ViolatedCoverIneqPool(data_S, y_SOL);
								  }
								}
							}
							/*if( IsDestNode(data_S, i) || IsDestNode(data_S, j))*/{
								//cutset complement
								outArc = false;
								inArc  = true;						
								CutSetGenerator(data_S, outArc, inArc);//find the arcs across cut set (going out from cutset nodes)
								CutSetMaxFlow(data_S, outArc, inArc, n_sc);//find the flow across cutset
								if(maxCutFlow>0){
								  OpenCloseArcs(data_S, y_SOL);//find the set of open and close arcs
								  if(FindMinCardinality(data_S, y_SOL)){
									  SortingSol(data_S, y_SOL);		
									  LiftingProcedure(data_S);
									  ViolatedCoverIneqPool(data_S, y_SOL);
								  }
								}
							}
						}
				if(getNumCoverIneq()>0)
				  cout << "Number of violated cardinality inequalities with cardinality of two: " << getNumCoverIneq() << endl;
			} 
		}
//////////////////////////////////////////////////////////////////////////////////////////////////
		void SortingCapacities(Data_S *data_S)
		{
			//ordering the variables in ascending fashion first
			vector <pair <double,int>  > aux ;
			for(int a = 0; a < data_S->getN_arcs() - data_S->getN_od(); a++)
				aux.push_back(make_pair(data_S->getU(a),a));
			sort(aux.begin(), aux.end(), PairCompare);//sort based on the double value
			//record variables order in descending fashion
			int t=0;			
			for(int a = aux.size()-1; a >=0; a--)
				capBasedOrder[t++] =aux[a].second;
			
		}
//////////////////////////////////////////////////////////////////////////////////////////////////
		bool FindMinCardinality(Data_S *data_S, IloNumArray &y_SOL)
		{
			bool cardinalityFound =false;
			//reset restricted cardinality size
			minCardinality = 0;
			//count minimum number of arcs to cover the demand in the restricted set
			if(/*minCardinality>0 &&*/ bigD>0){
				int arcId;
				double cumulativeCap=0;
				for(int a = 0; a < data_S->getN_arcs() - data_S->getN_od(); a++){
					arcId = capBasedOrder[a];
					if(cutSet[arcId] == 1 && closeArcs[arcId] == 0 && openArcs[arcId] == 0){					
							if(cumulativeCap + data_S->getU(arcId) < bigD){
								cumulativeCap += data_S->getU(arcId);
								minCardinality++;
							}
							else
								break;
						}					
				}		
				//if there is some demand left to cover
				if(cumulativeCap < bigD)
					minCardinality++;	
				cardinalityFound =true;
			}				
			//to check if the restricted cardinality is violated(is found or not); 
			/*double violCheck=0.0;
			for(int a = 0; a < data_S->getN_arcs() - data_S->getN_od(); a++)
				if(cutSet[arcId] == 1 && closeArcs[arcId] == 0 && openArcs[arcId] == 0)
					violCheck += y_SOL[a];	*/									
// 			if(violCheck > minCardinality)//if found restricted ineq is violated	
// 				cardinalityFound =false;
			//function return method	
			return cardinalityFound;
		}
//////////////////////////////////////////////////////////////////////////////////////////////////
		bool IJConnected(Data_S *data_S, int i, int j, IloNumArray &y_SOL)//see if between i and j an active arc exist 
		{
			bool twoNodesConnected=false;
			//check if two nodes i and j have an arc in between
			for(int a = 0; a < data_S->getN_arcs() - data_S->getN_od(); a++)
				if( (i == data_S->getArcs(a,0) && j == data_S->getArcs(a,1)) || (j == data_S->getArcs(a,0) && i == data_S->getArcs(a,1)) )
				  //if(y_SOL[a] > 0)
				  {
					  twoNodesConnected =true;
					  break;
				  }
			return twoNodesConnected;
		}
//////////////////////////////////////////////////////////////////////////////////////////////////
		bool ODPair(Data_S *data_S, int i, int j)
		{
			bool originDestinationPair =false;
			//check if they are OD pair
			for(int k = 0; k < data_S->getN_od(); k++)
				if( (i ==data_S->getOd(k,0) && j ==data_S->getOd(k,1)) || (j ==data_S->getOd(k,0) && i ==data_S->getOd(k,1))){
					originDestinationPair=true;
					break;
				}
			return originDestinationPair;
		}
//////////////////////////////////////////////////////////////////////////////////////////////////
		bool IsOriginNode(Data_S *data_S, int nodeID)
		{
			bool itIs=false;
			//check if the node is origin of at least one commodity 
			for(int k = 0; k < data_S->getN_od(); k++)
				if(nodeID == data_S->getOd(k,0))
					itIs =true;
			return itIs;
		}
/////////////////////////////////////////////		
		bool IsDestNode(Data_S *data_S, int nodeID)
		{
			bool itIs=false;
			for(int k = 0; k < data_S->getN_od(); k++)
				if(nodeID == data_S->getOd(k,1))
					itIs =true;
			return itIs;
		}
//////////////////////////////////////////////////////////////////////////////////////////////////
		void CutSetGenerator(Data_S *data_S, bool outArc, bool inArc)
		{
			//reinitialize the cutset array
			for(int a = 0; a < data_S->getN_arcs(); a++)
				cutSet[a] = 0;//arc a is not in the cut set
			//find the cut set
			for(int a = 0; a < data_S->getN_arcs() - data_S->getN_od(); a++)//dummies are not part of cutset
				{//not for the arcs in between			
					if(outArc)
					  if( (baseNodes[0] == data_S->getArcs(a,0) && baseNodes[1] != data_S->getArcs(a,1)) || (baseNodes[1] == data_S->getArcs(a,0) && baseNodes[0] != data_S->getArcs(a,1)) )//if we want the cutset
						cutSet[a] =1;
					if(inArc)
					  if( (baseNodes[0] == data_S->getArcs(a,1) && baseNodes[1] != data_S->getArcs(a,0) ) || (baseNodes[1] == data_S->getArcs(a,1) && baseNodes[0] != data_S->getArcs(a,0)) )//if we want its complement
						cutSet[a] =1;
				}	
		}
//////////////////////////////////////////////////////////////////////////////////////////////////
		void CutSetMaxFlow(Data_S *data_S, bool outArc, bool inArc, int n_sc)
		{
			//reinitialize the max flow across cutset
			maxCutFlow =0;
			double auxMaxCutFlow;
			//find the maximum flow over the cutset among all scenarios
			for(int s =0; s< n_sc; s++){
				auxMaxCutFlow =0;
				for(int k =0; k< data_S->getN_od(); k++){
					if(outArc)//for the outflow
					  if( (data_S->getOd(k,0) == baseNodes[0] && data_S->getOd(k,1) != baseNodes[1]) || (data_S->getOd(k,0) == baseNodes[1] && data_S->getOd(k,1) != baseNodes[0]) )//if one of the nodes is origin of this commodity, //if the destination is not inside the cutset 
						auxMaxCutFlow += data_S->getD(k,s);
					if(inArc)//for inward flow
					  if( (data_S->getOd(k,1) == baseNodes[0] && data_S->getOd(k,0) != baseNodes[1]) || (data_S->getOd(k,1) == baseNodes[1] && data_S->getOd(k,0) != baseNodes[0]) )//if one is destination of this commodity
						auxMaxCutFlow += data_S->getD(k,s);
				}
				if(auxMaxCutFlow > maxCutFlow)
					maxCutFlow = auxMaxCutFlow;
			}
		}
//////////////////////////////////////////////////////////////////////////////////////////////////
		void OpenCloseArcs(Data_S *data_S, IloNumArray &y_SOL)
		{
			bigU=0.0, bigD = maxCutFlow;
			for(int a = 0; a < data_S->getN_arcs() - data_S->getN_od(); a++)
				if(cutSet[a] == 1)
					bigU += data_S->getU(a);
			for(int a = 0; a < data_S->getN_arcs()  - data_S->getN_od() ; a++){
				closeArcs[a] =0;
				openArcs[a]  =0;//it is member of neither 
				if(cutSet[a] == 1){
					if( (y_SOL[a] <= downEpsilon) && (bigU - data_S->getU(a) >= bigD) ){
						closeArcs[a] =1;
						bigU -= data_S->getU(a);
					}
					if( (y_SOL[a] >= 1- upEpsilon) && (bigD - data_S->getU(a) > 0) ){
						openArcs[a]  =1;
						bigD -= data_S->getU(a);
						bigU -= data_S->getU(a);
					}
				}
			}
			if(bigU < bigD)
				cout << "\n Cardinality condition not satisfied!! \n";
		}
//////////////////////////////////////////////////////////////////////////////////////////////////
		void SortingSol(Data_S *data_S, IloNumArray &y_SOL)
		{		
			orderedSolution.clear();
			//ordering the variables in ascending fashion
			vector <pair <double,int>  > aux ;
			for(int a = 0; a < data_S->getN_arcs() - data_S->getN_od(); a++)
				aux.push_back(make_pair(y_SOL[a],a));
			sort(aux.begin(), aux.end(), PairCompare);//sort based on the double value
			for(int i = 0; i < aux.size(); i++)
				orderedSolution.push_back(aux[i].second);
			//refining the orders based on non-increasing orders of capacities				
			for(int a = 0; a < data_S->getN_arcs() - data_S->getN_od(); a++)
				for(int aa = a; aa < data_S->getN_arcs() - data_S->getN_od(); aa++)
					if(y_SOL[a] == y_SOL[aa] && data_S->getU(a) < data_S->getU(aa) )							
						swap(orderedSolution[a],orderedSolution[aa]);							
		}
//////////////////////////////////////////////////////////////////////////////////////////////////
		/*void LiftingProcedure(Data_S *data_S)
		{
			//update the coefficient and lifted sets  
			for(int a = 0; a < data_S->getN_arcs()- data_S->getN_od(); a++){
				liftedVar[a] =0;//initially no var is lifted
				if(cutSet[a] == 1 && closeArcs[a] == 0 && openArcs[a] == 0)
					coverCoeff[a] =1;
				else
					coverCoeff[a] =0;
			}
			//apply lifting down procedure (variable with smaller value, first)
			int arcId;
			for(int a = 0; a < data_S->getN_arcs() - data_S->getN_od(); a++){
				arcId = orderedSolution[a];
				if(cutSet[arcId] == 1 && closeArcs[arcId] == 0 && openArcs[arcId] == 1)//if it is in the cutset but not in the close and cover set
					coverCoeff[arcId] = LiftingDown(data_S, arcId);
			}
			//apply lifting up procedure (variable with larger value first)
			for(int a = data_S->getN_arcs() - data_S->getN_od() -1; a>=0; a--){
				arcId = orderedSolution[a];
				if( closeArcs[arcId] == 1 )
					coverCoeff[arcId] = LiftingUp(data_S, arcId); 			
			}					
		}*/
		void LiftingProcedure(Data_S *data_S)
		{
			//list of arcs for lifting
			liftList.clear();
			//update the coefficient and lifted sets  
			for(int a = 0; a < data_S->getN_arcs(); a++){
				liftedVar[a] =0;//initially no var is lifted
				if(cutSet[a] == 1 && closeArcs[a] == 0 && openArcs[a] == 0){
					coverCoeff[a] =1;
					liftList.push_back(a);
				}
				else
					coverCoeff[a] =0;
			}
			//apply lifting down procedure (variable with smaller value, first)
			int arcId;
			for(int a = 0; a < data_S->getN_arcs() - data_S->getN_od(); a++){
				arcId = orderedSolution[a];
				if(cutSet[arcId] == 1 && closeArcs[arcId] == 0 && openArcs[arcId] == 1){//if it is in the cutset but not in the close and cover set
					coverCoeff[arcId] = LiftingDown(data_S, arcId);
					liftList.push_back(arcId);
				}
			}
			//apply lifting up procedure (variable with larger value first)
			for(int a = data_S->getN_arcs() - data_S->getN_od() -1; a>=0; a--){
				arcId = orderedSolution[a];
				if( closeArcs[arcId] == 1 ){
					coverCoeff[arcId] = LiftingUp(data_S, arcId);
					liftList.push_back(arcId);
				}
			}					
		}
//////////////////////////////////////////////////////////////////////////////////////////////////
		double DynamicProg(Data_S *data_S, double Cap, int arcID)
		{
			if(arcID == liftList.size() )
			{
			  if(Cap > 0)
			    return 1.0e15;
			  else
			    return 0.0;
			}
			else
			{
			  if(Cap >= 0)
			    return min(coverCoeff[liftList[arcID]] + DynamicProg(data_S, Cap - data_S->getU(liftList[arcID]), arcID+1) , 
				     DynamicProg(data_S, Cap, arcID+1)
				     );
			  else
			    return 0.0;
			}
		}
//////////////////////////////////////////////////////////////////////////////////////////////////
		double LiftingDown(Data_S *data_S, int yVar)
		{
			//constraint right hand side
			double Cap = maxCutFlow + data_S->getU(yVar);
			for(int a = 0; a < data_S->getN_arcs() - data_S->getN_od(); a++)
				if(openArcs[a] == 1 && liftedVar[a] == 0)
					Cap -= data_S->getU(a);
			double objValue = DynamicProg(data_S, Cap, 0);
			//solve the problem and calculated the lifting coefficient
			double liftCoefficient= -minCardinality;			
			for(int a = 0; a < data_S->getN_arcs() - data_S->getN_od(); a++)
				if( liftedVar[a] == 1 && closeArcs[a] == 0)
					liftCoefficient -= coverCoeff[a];
			
			if (objValue < 1.0e15){
				liftCoefficient += objValue;
			}						
			else{			
				for(int a = 0; a < data_S->getN_arcs() - data_S->getN_od(); a++)
					if((cutSet[a] == 1 && closeArcs[a] == 0 && openArcs[a] == 0) || (liftedVar[a]==1) )
						liftCoefficient += coverCoeff[a];
				liftCoefficient++;
					
			}					
			//update the list of lifted variables so far
			liftedVar[yVar] =1;
			//return the lifting coefficient of the lifted variable
			return liftCoefficient;
		}
/////////////////////////////////////
		double LiftingUp(Data_S *data_S, int yVar)
		{	
			//constraint right hand side
			double Cap = maxCutFlow - data_S->getU(yVar);
			for(int a = 0; a < data_S->getN_arcs() - data_S->getN_od(); a++)
				if(openArcs[a] == 1 && liftedVar[a] == 0)
					Cap -= data_S->getU(a);
			double objValue = DynamicProg(data_S, Cap, 0);
			//solve the problem 
			double liftCoefficient =minCardinality;			
			for(int a = 0; a < data_S->getN_arcs() - data_S->getN_od(); a++)
				if(liftedVar[a] == 1 && closeArcs[a] == 0)
					liftCoefficient += coverCoeff[a];
			if( objValue < 1.0e15){
				liftCoefficient -= objValue ;
			}
			else{
				for(int a = 0; a < data_S->getN_arcs() - data_S->getN_od(); a++)
					if((cutSet[a] == 1 && closeArcs[a] == 0 && openArcs[a] == 0) || (liftedVar[a]==1) )
						liftCoefficient -= coverCoeff[a];
				liftCoefficient--;
			}					
			//update the list of lifted variables so far
			liftedVar[yVar] =1;
			//return the lifting coefficient of the lifted variable
			return liftCoefficient;
		}
		
		/*double LiftingDown(Data_S *data_S, int yVar)
		{
			//build objFun
			IloExpr expr(env);
			for(int a =0; a < data_S->getN_arcs() - data_S->getN_od(); a++)
				if((cutSet[a] == 1 && closeArcs[a] == 0 && openArcs[a] == 0) || (liftedVar[a]==1) )
					expr += coverCoeff[a] * y[a];
			objFun.setExpr(expr);
			expr.end();
			//create constraint
			(model).remove(con);
			(con).endElements();
			IloExpr expr1(env);			
			for(int a = 0; a < data_S->getN_arcs() - data_S->getN_od(); a++)
				if((cutSet[a] == 1 && closeArcs[a] == 0 && openArcs[a] == 0) || (liftedVar[a]==1) )
					expr1 += data_S->getU(a) * y[a];
			//constraint right hand side
			double RHS = maxCutFlow + data_S->getU(yVar);
			for(int a = 0; a < data_S->getN_arcs() - data_S->getN_od(); a++)
				if(openArcs[a] == 1 && liftedVar[a] == 0)
					RHS -= data_S->getU(a);
			(con).add(IloRange(env, RHS, expr1 , IloInfinity));
			(model).add(con);
			expr1.end();
			//solve the problem and calculated the lifting coefficient
			double liftCoefficient= -minCardinality;			
			for(int a = 0; a < data_S->getN_arcs() - data_S->getN_od(); a++)
				if( liftedVar[a] == 1 && closeArcs[a] == 0)
					liftCoefficient -= coverCoeff[a];
			
			if ((cplex).solve()){
				liftCoefficient += (cplex).getObjValue();
			}						
			else{			
				for(int a = 0; a < data_S->getN_arcs() - data_S->getN_od(); a++)
					if((cutSet[a] == 1 && closeArcs[a] == 0 && openArcs[a] == 0) || (liftedVar[a]==1) )
						liftCoefficient += coverCoeff[a];
				liftCoefficient++;
					
			}					
			//update the list of lifted variables so far
			liftedVar[yVar] =1;
			//return the lifting coefficient of the lifted variable
			return liftCoefficient;
		}
/////////////////////////////////////
		double LiftingUp(Data_S *data_S, int yVar)
		{	
			//build objFun
			IloExpr expr(env);
			for(int a =0; a < data_S->getN_arcs() - data_S->getN_od(); a++)
				if((cutSet[a] == 1 && closeArcs[a] == 0 && openArcs[a] == 0) || (liftedVar[a]==1) )
					expr += coverCoeff[a] * y[a];
			objFun.setExpr(expr);
			expr.end();
			//create constraint
			(model).remove(con);
			(con).endElements();
			IloExpr expr1(env);			
			for(int a = 0; a < data_S->getN_arcs() - data_S->getN_od(); a++)
				if((cutSet[a] == 1 && closeArcs[a] == 0 && openArcs[a] == 0) || (liftedVar[a]==1) )
					expr1 += data_S->getU(a) * y[a];
			//constraint right hand side
			double RHS = maxCutFlow - data_S->getU(yVar);
			for(int a = 0; a < data_S->getN_arcs() - data_S->getN_od(); a++)
				if(openArcs[a] == 1 && liftedVar[a] == 0)
					RHS -= data_S->getU(a);
			(con).add(IloRange(env, RHS, expr1 , IloInfinity));
			(model).add(con);
			expr1.end();
			//solve the problem 
			double liftCoefficient =minCardinality;			
			for(int a = 0; a < data_S->getN_arcs() - data_S->getN_od(); a++)
				if(liftedVar[a] == 1 && closeArcs[a] == 0)
					liftCoefficient += coverCoeff[a];
			if((cplex).solve()){
				liftCoefficient -= (cplex).getObjValue() ;
			}
			else{
				for(int a = 0; a < data_S->getN_arcs() - data_S->getN_od(); a++)
					if((cutSet[a] == 1 && closeArcs[a] == 0 && openArcs[a] == 0) || (liftedVar[a]==1) )
						liftCoefficient -= coverCoeff[a];
				liftCoefficient--;
			}					
			//update the list of lifted variables so far
			liftedVar[yVar] =1;
			//return the lifting coefficient of the lifted variable
			return liftCoefficient;
		}*/
//////////////////////////////////////////////////////////////////////////////////////////////////
		void ViolatedCoverIneqPool(Data_S *data_S, IloNumArray &y_SOL)
		{
			//check for violated of the cut
			double violation =1;
			double RHS = minCardinality;
			for(int a = 0; a < data_S->getN_arcs() - data_S->getN_od(); a++)
				if(cutSet[a] == 1){
					violation -= coverCoeff[a]*y_SOL[a];
					if(openArcs[a] == 1){
						violation += coverCoeff[a];
						RHS += coverCoeff[a];
					}
				}				
			//store the cut if violated
			if(violation > violationThreshhold){
// 				cout << "\n Violation is: " << violation <<endl;
				//record the coefficients
				vector <double> aux;
				for(int a = 0; a < data_S->getN_arcs() ; a++)
					aux.push_back(coverCoeff[a]);
				//record the right hand side of the ineq
				aux.push_back(RHS);
				//shot it into the pool
				IneqPool.push_back(aux);
			}
		}
//////////////////////////////////////////////////////////////////////////////////////////////////
		int getNumCoverIneq() {return IneqPool.size();}
		vector <double > getCoverIneqCoeff(int q) {return IneqPool[q];}
//////////////////////////////////////////////////////////////////////////////////////////////////
		~class_MCI()
		{
		}


};//WND class



#endif
