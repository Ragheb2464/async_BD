#ifndef CLASS_NC_INEQ_H
#define CLASS_NC_INEQ_H





ILOSTLBEGIN
#define Comm COMM_WORLD



class class_Network_Connectivity 
{ 
  private:    
		vector <int > PureTransferNode;
		vector <int > PureOriginrNode;
		vector <int > PureDestinationNode;
		vector <vector< double> > cutPool;
			  
  public:        		
		class_Network_Connectivity(int n_sc, MasterSolChoice *Sol_Choice, Data_S *data_S, Search_Param	*search_param)		
		{
		  TransferSet(data_S);
		  OriginSet(data_S);
		  DestinationSet(data_S);
		}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		void TransferSet(Data_S *data_S)
		{
		  bool ItIs;
		  for(int i = 0; i < data_S->getN_nodes(); i++){
		    ItIs =true;
		    for(int k = 0; k < data_S->getN_od(); k++)
		      if(i ==data_S->getOd(k,0) || i ==data_S->getOd(k,1))
				  ItIs =false;
		      if(ItIs)
			  PureTransferNode.push_back(i);
		  }
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
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		void SingleCuts(Data_S *data_S, int n_sc)
		{
		    vector <double> aux;
		    int nodeID;
		    //for the origin set
		    for(int i=0; i< getPureOriginNodeSize(); i++){
			aux.clear();
			nodeID=getOriginNode(i);
			for(int a =0; a < data_S->getN_arcs() - data_S->getN_od() ; a++){
			    if(data_S->getArcs(a,0) == nodeID)
				aux.push_back(data_S->getU(a));			
			    else if(data_S->getArcs(a,1) == nodeID)
				aux.push_back(-data_S->getU(a));			
			    else/*(data_S->getArcs(a,0) != nodeID && data_S->getArcs(a,1) != nodeID)*/
				aux.push_back(0.0);
			}
			//get the max demand originating from this node
			aux.push_back( ceil(OriginMaximumFlow(nodeID, data_S, n_sc)) );
			cutPool.push_back(aux);	
		    }
		    //for the destination set
		   for(int i=0; i< getPureDestinationNodeSize() ; i++){
			  aux.clear();
			  for(int a =0; a < data_S->getN_arcs() - data_S->getN_od() ; a++){
			      if(data_S->getArcs(a,1) == getdestinationNode(i))
				  aux.push_back(data_S->getU(a));			
			      else if(data_S->getArcs(a,0) == getdestinationNode(i))
				  aux.push_back(-data_S->getU(a));			
			      else
				  aux.push_back(0);
			  }
			  //get the max demand originating from this node
			  aux.push_back( -DestinationMaximumFlow(getdestinationNode(i), data_S, n_sc) );
			  cutPool.push_back(aux);			
		    }		      		    
		}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		double OriginMaximumFlow(int nodeID, Data_S *data_S, int n_sc)
		{
		    double maxDemand=0.0;
		    double auxMaxD=0.0;
		    for(int s =0; s < n_sc; s++){
		      auxMaxD =0.0;
		      for(int k =0; k< data_S->getN_od(); k++)
			    if(data_S->getOd(k,0) == nodeID)
				auxMaxD += data_S->getD(k,s);
			  
		      if(auxMaxD > maxDemand)
			  maxDemand = auxMaxD;
		    }
 		    cout << nodeID << "  " << maxDemand << endl;
		    return maxDemand;		  
		}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		double DestinationMaximumFlow(int nodeID, Data_S *data_S, int n_sc)
		{
		    double maxDemand=0;
		    double auxMaxD=0;
		    for(int s=0; s<n_sc; s++){
		      auxMaxD=0;
		      for(int k=0; k< data_S->getN_od(); k++)
			  if(data_S->getOd(k,1) == nodeID)
			    auxMaxD += data_S->getD(k,s);
		      if(auxMaxD > maxDemand)
			  maxDemand = auxMaxD;
		    }
		    
		    return maxDemand;		  
		}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		int getPureTransferNodeSize() 		{return PureTransferNode.size();}
		int getTransferNode(int i) 		{return PureTransferNode[i];}
		
		int getPureOriginNodeSize()		{return PureOriginrNode.size();}
		int getOriginNode(int i) 		{return PureOriginrNode[i];}
		
		int getPureDestinationNodeSize() 	{return PureDestinationNode.size();}
		int getdestinationNode(int i) 		{return PureDestinationNode[i];}
		
		int getCutPoolSize()			{return cutPool.size();}
		vector<double>& getCut(int q)		{return cutPool[q];}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

		~class_Network_Connectivity()
		{
		}


};//WND class



#endif
