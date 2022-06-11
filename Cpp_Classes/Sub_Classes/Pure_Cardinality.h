#ifndef CLASS_PMCI_INEQ_H
#define CLASS_PMCI_INEQ_H

//#include <ilcplex/ilocplex.h>
ILOSTLBEGIN
#define Comm COMM_WORLD


class class_PMCI
{ 
  private:  
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
			vector <int> baseNodes;//the selected nodes to be the cutset
			vector <vector< double> > IneqPool;// to keep the generated violated cuts
			double bigU, bigD;//for the residual demand and capacity 
			int *capBasedOrder;//to order arcs based on their capacity
			double minCardinality;
			vector<int > liftList;//the list containg the arcs to be lifted			
			vector <int > PureOriginrNode;//set of nodes which are origin but not destination
			vector <int > PureDestinationNode;//the opposite same
  public:        		
		class_PMCI(Data_S *data_S, int n_sc)		
		{	
 			//initialization of the required arraies and vectors 
			openArcs  =new int[data_S->getN_arcs()]();
			closeArcs =new int[data_S->getN_arcs()]();
			upEpsilon = 0.4999;
			downEpsilon = 0.4999;
			violationThreshhold= 1.0e-3;
			coverSet 	= new int[data_S->getN_arcs()]();
			coverCoeff 	= new double[data_S->getN_arcs()](); 
			liftedVar	= new int[data_S->getN_arcs()]();
			cutSet 		= new int[data_S->getN_arcs()]();
// 			baseNodes	= new int[2];//
			capBasedOrder 	= new int[data_S->getN_arcs()]();
			SortingCapacities(data_S);//order design variables based on their capacity in non-increasing order
			OriginSet(data_S);//derive the origin set
			DestinationSet(data_S);
		}
		
//////////////////////////////////////////////////////////////////////////////////////////////////
		void MainLoop(Data_S *data_S, int n_sc, int cutSetCardinality)
		{
			bool outArc, inArc;
			//clean the cut pool
			IneqPool.clear();
			//for the single note cutset 
				for(int i=0; i< data_S->getN_nodes(); i++){	
					if( IsOriginNode(data_S, i)){
						baseNodes.clear();
						baseNodes.push_back(i);
						baseNodes.push_back(i);
// 						baseNodes[0] = i;
// 						baseNodes[1] = i;//our cutset nodes
						//cutset
						outArc = true;
						inArc  = false;						
						CutSetGenerator(data_S, outArc, inArc);//find the arcs across cut set (going out from cutset nodes)
						CutSetMaxFlow(data_S, outArc, inArc, n_sc);//find the flow across cutset
						// OpenCloseArcs(data_S);//find the set of open and close arcs
						if(FindMinCardinality(data_S)){//lift violated ineq		
							LiftingProcedure(data_S); 
							ViolatedCoverIneqPool(data_S);
						}						
					}
					if(IsDestNode(data_S, i) ){
						baseNodes.clear();
						baseNodes.push_back(i);
						baseNodes.push_back(i);
// 						baseNodes[0] = i;
// 						baseNodes[1] = i;//our cutset nodes
						//cutset
						outArc = false;
						inArc  = true;						
						CutSetGenerator(data_S, outArc, inArc);//find the arcs across cut set (going out from cutset nodes)
						CutSetMaxFlow(data_S, outArc, inArc, n_sc);//find the flow across cutset
						// OpenCloseArcs(data_S);//find the set of open and close arcs
						if(FindMinCardinality(data_S)){//lift violated ineq		
							LiftingProcedure(data_S); 
							ViolatedCoverIneqPool(data_S);
						}
					} 
				}
				double singleCarCut = getNumCarIneq();
				if(singleCarCut>0)
				  cout << "Number of initial cardinality inequalities with single cutset: " << singleCarCut << endl;					
			//for 2 we should make sure they are not a OD pair 
			if(cutSetCardinality == 2 ) {
				for(int i=0; i< data_S->getN_nodes(); i++)
					for(int j =i+1; j< data_S->getN_nodes(); j++)
						if(IJConnected(data_S,i,j) /*&& !ODPair(data_S,i,j)*/ ){//if there is at least one arc between them, they are a OD piar, at least one is origin or destination							
// 							baseNodes[0] = i;
// 							baseNodes[1] = j;
							baseNodes.clear();
							baseNodes.push_back(i);
							baseNodes.push_back(j);
							/*if( IsOriginNode(data_S, i) || IsOriginNode(data_S, j))*/{//if at least one in origin node (to have non zero outgoing flow)
								//cutset
								outArc = true;
								inArc  = false;						
								CutSetGenerator(data_S, outArc, inArc);//find the arcs across cut set (going out from cutset nodes)
								CutSetMaxFlow(data_S, outArc, inArc, n_sc);//find the flow across cutset
								if(maxCutFlow>0){
								  // OpenCloseArcs(data_S);//find the set of open and close arcs
								  if(FindMinCardinality(data_S)){
									  //SortingSol(data_S, y_SOL);		
									  LiftingProcedure(data_S);
									  ViolatedCoverIneqPool(data_S);
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
								  // OpenCloseArcs(data_S);//find the set of open and close arcs
								  if(FindMinCardinality(data_S)){
									  // SortingSol(data_S);		
									  LiftingProcedure(data_S);
									  ViolatedCoverIneqPool(data_S);
								  }
								}
							}
						}
			}
			double doubleCarCut= getNumCarIneq() - singleCarCut;
			if(doubleCarCut>0)
				  cout << "Number of initial cardinality inequalities with cardinality of two: " << doubleCarCut<< endl;
			//for the set of origin nodes
			if(PureOriginrNode.size() >=3){
				baseNodes.clear();
				for(int i=0; i< PureOriginrNode.size(); i++)
				  baseNodes.push_back(PureOriginrNode[i]);
				outArc = true;
				inArc  = false;						
				CutSetGenerator(data_S, outArc, inArc);//find the arcs across cut set (going out from cutset nodes)
				CutSetMaxFlow(data_S, outArc, inArc, n_sc);//find the flow across cutset
				if(maxCutFlow>0){
				// OpenCloseArcs(data_S);//find the set of open and close arcs
					if(FindMinCardinality(data_S)){
					  //SortingSol(data_S, y_SOL);		
						  LiftingProcedure(data_S);
						  ViolatedCoverIneqPool(data_S);
					}
				}
			}
			//for the set of destination nodes
			if(PureDestinationNode.size() >=3){
				baseNodes.clear();
				for(int i=0; i< PureDestinationNode.size(); i++)
				  baseNodes.push_back(PureDestinationNode[i]);
				outArc = false;
				inArc  = true;						
				CutSetGenerator(data_S, outArc, inArc);//find the arcs across cut set (going out from cutset nodes)
				CutSetMaxFlow(data_S, outArc, inArc, n_sc);//find the flow across cutset
				if(maxCutFlow>0){
				// OpenCloseArcs(data_S);//find the set of open and close arcs
					if(FindMinCardinality(data_S)){
					  //SortingSol(data_S, y_SOL);		
						  LiftingProcedure(data_S);
						  ViolatedCoverIneqPool(data_S);
					}
				}
			}
			if( getNumCarIneq() - singleCarCut - doubleCarCut >0)
			  cout << "number of Origin/destination set cuts: " << getNumCarIneq() - singleCarCut - doubleCarCut << endl;
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
		bool FindMinCardinality(Data_S *data_S)
		{		  
			bool cardinalityFound =false;
			//reset restricted cardinality size
			minCardinality = 0;
			//count minimum number of arcs to cover the demand in the restricted set
			if(/*minCardinality>0 &&*/ maxCutFlow>0){
				int arcId;
				double cumulativeCap =0.0;
				for(int a =0; a < data_S->getN_arcs() - data_S->getN_od(); a++){
					arcId = capBasedOrder[a];
					if(cutSet[arcId] == 1 )					
						if(cumulativeCap + data_S->getU(arcId) < maxCutFlow){
							cumulativeCap += data_S->getU(arcId);
							minCardinality++;
						}
						else
							break;									
				}		
				//if there is some demand left to cover
				if(cumulativeCap < maxCutFlow)
					minCardinality++;	
				cardinalityFound =true;
			}				
			//function return method	
			return cardinalityFound;
		}
//////////////////////////////////////////////////////////////////////////////////////////////////
		bool IJConnected(Data_S *data_S, int i, int j)//see if between i and j an active arc exist 
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
				cutSet[a] =0;//arc a is not in the cut set
			//find the cut set
			for(int a = 0; a < data_S->getN_arcs() - data_S->getN_od(); a++)
			    if( AcceptArc(data_S, a, outArc, inArc) ) 
				cutSet[a] =1;	
		}
//////////////////////////////////////////////////////////////////////////////////////////////////
		bool AcceptArc(Data_S *data_S, int a, bool outArc, bool inArc)
		{
		    bool accept=false;
		    //outward arc
		    if(outArc)
		    {
			for(int i =0; i< baseNodes.size(); i++)
			  if(baseNodes[i] == data_S->getArcs(a,0)){
			    accept=true;
			    for(int j =0; j< baseNodes.size(); j++)
			      if( data_S->getArcs(a,1) == baseNodes[j])
				accept =false;
			  }
		    }
		    //inward arc
		    if(inArc)
		    {
			for(int i =0; i< baseNodes.size(); i++)
			  if(baseNodes[i] == data_S->getArcs(a,1)){
			    accept=true;
			    for(int j =0; j< baseNodes.size(); j++)
			      if( data_S->getArcs(a,0) == baseNodes[j])
				accept =false;
			  }
		    }
		    
		    return accept;
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
				for(int k =0; k< data_S->getN_od(); k++)
				  if(AcceptCommod(data_S, outArc, inArc, k))
				    auxMaxCutFlow += data_S->getD(k,s);
				  
				if(auxMaxCutFlow > maxCutFlow)
					maxCutFlow = auxMaxCutFlow;
			}
		}
//////////////////////////////////////////////////////////////////////////////////////////////////
		bool AcceptCommod(Data_S *data_S, bool outArc, bool inArc, int k)
		{
		  bool accept =false;
		  //origin
		  if(outArc)
		  {
			for(int i =0; i< baseNodes.size(); i++)
			  if(baseNodes[i] == data_S->getOd(k,0)){
			    accept=true;
			    for(int j =0; j< baseNodes.size(); j++)
			      if( data_S->getOd(k,1) == baseNodes[j])
				accept =false;
			  }		    
		  }
		  //destination
		  if(inArc)
		  {
			for(int i =0; i< baseNodes.size(); i++)
			  if(baseNodes[i] == data_S->getOd(k,1)){
			    accept=true;
			    for(int j =0; j< baseNodes.size(); j++)
			      if( data_S->getOd(k,0) == baseNodes[j])
				accept =false;
			  }		    
		  
		  }
		    return accept;
		}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		void OriginSet(Data_S *data_S)
		{
		  bool ItIs;
		  for(int i = 0; i < data_S->getN_nodes(); i++){
		    ItIs =true;
		    for(int k = 0; k < data_S->getN_od(); k++)//if i is destination for at least one commodity, it is not in the set
		      if(i == data_S->getOd(k,1))
				  ItIs =false;
		      if(ItIs)
			PureOriginrNode.push_back(i);
		  }
		}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		void DestinationSet(Data_S *data_S)
		{
		  bool ItIs;
		  for(int i = 0; i < data_S->getN_nodes(); i++){
		    ItIs =true;
		    for(int k = 0; k < data_S->getN_od(); k++)//if i is origin for at least one commodity, it is not in the set
		      if(i == data_S->getOd(k,0))
				  ItIs =false;
		      if(ItIs)
			PureDestinationNode.push_back(i);
		  }
		}
//////////////////////////////////////////////////////////////////////////////////////////////////
		void LiftingProcedure(Data_S *data_S)
		{
			//update the coefficient and lifted sets  
			for(int a = 0; a < data_S->getN_arcs(); a++){
				if(cutSet[a] == 1 )
					coverCoeff[a] =1;
				else
					coverCoeff[a] =0;
			}					
		}
//////////////////////////////////////////////////////////////////////////////////////////////////
		void ViolatedCoverIneqPool(Data_S *data_S)
		{			
			double RHS = minCardinality;
			vector <double> aux;
			for(int a = 0; a < data_S->getN_arcs() ; a++)
				aux.push_back(coverCoeff[a]);
			//record the right hand side of the ineq
			aux.push_back(RHS);
			//shot it into the pool
			IneqPool.push_back(aux);
			
		}
//////////////////////////////////////////////////////////////////////////////////////////////////
		int getNumCarIneq() {return IneqPool.size();}
		vector <double > &getCarIneqCoeff(int q) {return IneqPool[q];}
//////////////////////////////////////////////////////////////////////////////////////////////////
		~class_PMCI()
		{
		}


};//WND class



#endif
