#ifndef CLASS_FLOW_CUT_H
#define CLASS_FLOW_CUT_H


ILOSTLBEGIN



class FlowCuts
{ 
	private:
		int *baseNodes;//the selected nodes to be the cutset
		vector <int> cutSet, complementoryCutSet, commSet, setC1, setC2, setD1, setD2, FPItoFCIset; //the arcs in the cutset
		double demandSum;
		vector <double> scenDemandSum;
		int numFPI;
		vector <double> fixPart;
		vector <vector  <double> > scenFixPart;
		vector < vector <double > >  arcCoeff, conCoeff;
		vector < vector <double > >  scenArcCoeff;
	public:
		FlowCuts() 
		{	
		  baseNodes = new int[2];
		}
//////////////////////////////////
		void flowPack(int n_sc, Data_S* data_S, double* y_SOL, double **X_Value, int sID)
		{
			bool outArc, inArc;
			numFPI=0;
			fixPart.clear();
			scenFixPart.clear();
			arcCoeff.clear();
			conCoeff.clear();
			//cardinality one
			for(int i=0; i< data_S->getN_nodes(); i++){	
				if( IsOriginNode(data_S, i) /*|| IsDestNode(data_S, i)*/){//for outgoing arcs from cutset
					baseNodes[0] = i;
					baseNodes[1] = i;//our cutset nodes
					//cutset
					outArc = true;
					inArc  = false;						
					CutSetGenerator(data_S, outArc, inArc);//find the arcs across cut set (going out from cutset nodes)
					for(int a=0; a <cutSet.size(); a++){
					    commoditySet(n_sc, cutSet[a], sID, X_Value, data_S);
					    findSetC1andC2(cutSet[a], X_Value, y_SOL, data_S);
					    double uViol=vioCut(cutSet[a],sID, y_SOL, data_S);
					    if(uViol < -1e-17){
					      liftCutCase1(cutSet[a], uViol, y_SOL, X_Value, data_S);
					      recordFPI(n_sc, uViol, sID, cutSet[a], data_S, y_SOL, X_Value);
					    }
 					    else if (uViol + min(demandSum,data_S->getU(cutSet[a])) > 1e-17 ){
 					      genFCIfromFPI(cutSet[a], uViol, data_S, y_SOL, X_Value);
 					      recordFCI(n_sc, uViol, sID, cutSet[a], data_S, y_SOL, X_Value);
 					    }
					}  
				} 				
			}
			
		}
///////////////////////
		void genFCIfromFPI(int arcID, double& uViol, Data_S *data_S, double* y_SOL, double **X_Value)
		{
		  setD2.clear();
		  FPItoFCIset.clear();
		  setC1.push_back(arcID);
		  uViol += min(demandSum,data_S->getU(arcID));
		  
//  		  for(int a=0; a < setC1.size(); a++)
//  		    if( min(demandSum,data_S->getU(setC1[a])) > 1e-17 + uViol ){
//  		      FPItoFCIset.push_back(setC1[a]);
// 		    }
		  
		  for(int a=0; a<complementoryCutSet.size(); a++)
		    if(!(find(setC2.begin(), setC2.end(), complementoryCutSet[a]) != setC2.end()))
		      if(sumX(complementoryCutSet[a], X_Value) > min(uViol,min(demandSum,data_S->getU(complementoryCutSet[a]))) * y_SOL[complementoryCutSet[a]]  )
			  setD2.push_back(complementoryCutSet[a]);
		  
		}
//////////////////////
		void recordFCI(int n_sc, double uViol, int sID, int arcID, Data_S *data_S, double* y_SOL, double **X_Value)
		{
		    double aux=0, auxaux;
		    double viol=0;
		    vector <double> auxVec (data_S->getN_arcs(),0);
		    //fix part of the cut
		    for(int a=0; a< FPItoFCIset.size(); a++){
		      auxaux 					= min(data_S->getU(FPItoFCIset[a]), demandSum);
		      aux 						+= auxaux;
		      auxVec[FPItoFCIset[a]] 	-= auxaux;
		      viol 						-= auxaux * y_SOL[FPItoFCIset[a]];
		    }
		    for(int a=0; a< setC2.size(); a++)
		      aux -= min(data_S->getU(setC2[a]), demandSum);
		    aux -= demandCutSetL(sID,data_S);
		    for(int a=0; a< setC1.size(); a++){
		      auxaux 			= max(0.0, min(data_S->getU(setC1[a]), demandSum) - uViol);
		      aux 				+=  auxaux;
		      auxVec[setC1[a]] 	-= auxaux;
		      viol 				-= auxaux * y_SOL[setC1[a]];
		    }
		    viol += aux;
		    
			vector <double> auxScenFixPart;
			for(int s=0; s<n_sc; s++){//for scenario when we propagate cuts
				double auxS=0.0;
				for(int a=0; a< FPItoFCIset.size(); a++)
					auxS += min(data_S->getU(FPItoFCIset[a]), scenDemandSum[s]);
				for(int a=0; a< setC2.size(); a++)
					auxS -= min(data_S->getU(setC2[a]), scenDemandSum[s]);
				auxS -= demandCutSetL(s,data_S);
				 for(int a=0; a< setC1.size(); a++)
					auxS += max(0.0, min(data_S->getU(setC1[a]), scenDemandSum[s]) - uViol);
				auxScenFixPart.push_back(auxS);
			}
			
			//arc coeff
 		    for(int a=0; a< setD2.size(); a++){
		      auxaux 			= min(uViol, min(data_S->getU(setD2[a]), demandSum) );
 		      auxVec[setD2[a]] 	-= auxaux;
 		      viol 				-= auxaux * y_SOL[setD2[a]];
 		    }
			
			vector <vector <double> > auxauxVec;
 		    for(int s=0; s<n_sc; s++){//for scenario when we propagate cuts
				vector <double> auxauxauxVec (data_S->getN_arcs(),0);
				for(int a=0; a< FPItoFCIset.size(); a++)
					auxauxauxVec[FPItoFCIset[a]] 	-= min(data_S->getU(FPItoFCIset[a]), scenDemandSum[s]);
				for(int a=0; a< setC1.size(); a++)
					auxauxauxVec[setC1[a]] 	-= max(0.0, min(data_S->getU(setC1[a]), scenDemandSum[s]) - uViol);
				for(int a=0; a< setD2.size(); a++)
					auxauxauxVec[setD2[a]] 	-= min(uViol, min(data_S->getU(setD2[a]), scenDemandSum[s]) );
				auxauxVec.push_back(auxauxauxVec);
			}
			
 		    if(viol > 1e-3) 
		      cout << "something is wrong when deriving the FCI: Danger " << viol << endl;
		    
		    //con coeff
		    vector <double> auxVec1;
		    int t=0;
		    for ( int a = 0; a < data_S->getN_arcs() ; a++)
			if (data_S->getArcs(a,1) != data_S->getArcs(a,0))//if it is not looping arc
				for ( int k = 0; k < data_S->getN_od(); k++){
				    auxVec1.push_back(0);
				    if( find(commSet.begin(), commSet.end(), k) != commSet.end() ){
						if(  find(setC1.begin(), setC1.end(), a) != setC1.end() ){
							auxVec1[t] = -1;
							viol += X_Value[a][k];					    
						}
					else if( find(complementoryCutSet.begin(), complementoryCutSet.end(), a) != complementoryCutSet.end() )
					      if(!(find(setC2.begin(), setC2.end(), a) != setC2.end()) )
							if(! (find(setD2.begin(), setD2.end(), a) != setD2.end())  ){
								auxVec1[t] = 1;
								viol -= X_Value[a][k];
							}						
				    }
				    t++;
				}
		    if(viol > 1e-1){
	//  			cout << sID << " violation of the derived FCI: " << viol << endl;
				conCoeff.push_back(auxVec1);
				arcCoeff.push_back(auxVec);
				fixPart.push_back(aux);
				scenFixPart.push_back(auxScenFixPart);
				for(int s=0; s<n_sc; s++)
					scenArcCoeff.push_back(auxauxVec[s]);
				numFPI++;
		    }		  
		}
//////////////////////
		void recordFPI(int n_sc, double uViol, int sID, int arcID, Data_S *data_S, double* y_SOL, double **X_Value)
		{
		    double aux=0, auxaux;
		    double viol=0;
		    vector <double> auxVec (data_S->getN_arcs(),0);
		    //fix part of the cut
		    for(int a=0; a< setC1.size(); a++){
		      aux -= min(data_S->getU(setC1[a]), demandSum);
		      /*if(y_SOL[setC1[a]] < 1e-5){
			  double auxLift = liftCoeff(setC1[a], uViol, data_S) - min(data_S->getU(setC1[a]), demandSum);
			  aux += auxLift;
			  cout << "lift coeff is " << auxLift << endl;
			  auxVec[setC1[a]] -= auxLift;
			  viol -= auxLift*y_SOL[setC1[a]];
		      }*/
		    }
			vector <double> auxScenFixPart;
			for(int s=0; s<n_sc; s++){//for scenario when we propagate cuts
				double auxS=0.0;
				for(int a=0; a< setC1.size(); a++)
					auxS -= min(data_S->getU(setC1[a]), scenDemandSum[s]);
				auxScenFixPart.push_back(auxS);
			}
		    for(int a=0; a< setC2.size(); a++){
		      auxaux 			= max(0.0, uViol + min(data_S->getU(setC2[a]), demandSum) );
		      aux				+= auxaux;
		      auxVec[setC2[a]] 	-= auxaux;
		      viol 				-= auxaux * y_SOL[setC2[a]];
		    }
		    viol += aux;
		    for(int s=0; s<n_sc; s++)//for other scenarios when we propagate 
				for(int a=0; a< setC2.size(); a++)
				  auxScenFixPart[s] += max(0.0, uViol + min(data_S->getU(setC2[a]), scenDemandSum[s]) );				
			
			//arc coeff
 		    for(int a=0; a< setD1.size(); a++){
		      auxaux 		= min(-uViol, min(data_S->getU(setD1[a]), demandSum) );
 		      auxVec[setD1[a]] 	-= auxaux;
 		      viol 		-= auxaux * y_SOL[setD1[a]];
 		    }
			vector <vector <double> > auxauxVec;			
			for(int s=0; s<n_sc; s++){//for cut propagation 
				vector <double> auxauxauxVec (data_S->getN_arcs(),0);
				for(int a=0; a< setC2.size(); a++)
					auxauxauxVec[setC2[a]] -= max(0.0, uViol + min(data_S->getU(setC2[a]), scenDemandSum[s]) );
				for(int a=0; a< setD1.size(); a++)
					auxauxauxVec[setD1[a]] -= min(-uViol, min(data_S->getU(setD1[a]), scenDemandSum[s]) );
				auxauxVec.push_back(auxauxauxVec);
			}
		    //con coeff
		    vector <double> auxVec1;
		    int t=0;
		    for ( int a = 0; a < data_S->getN_arcs() ; a++)
			if (data_S->getArcs(a,1) != data_S->getArcs(a,0))//if it is not looping arc
				for ( int k = 0; k < data_S->getN_od(); k++){
				    auxVec1.push_back(0);
				    if( find(commSet.begin(), commSet.end(), k) != commSet.end() ){
					if(  find(setC1.begin(), setC1.end(), a) != setC1.end() ){
					    auxVec1[t] = -1;
					    viol += X_Value[a][k];					    
					}
					else if( (find(complementoryCutSet.begin(), complementoryCutSet.end(), a) != complementoryCutSet.end())  && !(find(setC2.begin(), setC2.end(), a) != setC2.end()) ){
					    auxVec1[t] = 1;
					    viol -= X_Value[a][k];
					}	
					else if(find(setD1.begin(), setD1.end(), a) != setD1.end()){
					  auxVec1[t] = -1;
					  viol += X_Value[a][k];
					}
				    }
				    t++;
				}
		    if(viol > 1e-1){
	//  			cout << sID << " violation of the FPI: " << viol << endl;
				conCoeff.push_back(auxVec1);
				arcCoeff.push_back(auxVec);
				fixPart.push_back(aux);
				scenFixPart.push_back(auxScenFixPart);
				for(int s=0; s<n_sc; s++)
					scenArcCoeff.push_back(auxauxVec[s]);
				numFPI++;
		    }
		}
//////////////////////////
		double vioCut(int arcID, int sID, double* y_SOL, Data_S *data_S)
		{
		  double uViol=0.0;
		  for(int a=0; a < setC1.size(); a++)
		     uViol += min(data_S->getU(setC1[a]), demandSum);
		  for(int a=0; a < setC2.size(); a++)
		    uViol -= min(data_S->getU(setC2[a]), demandSum);
		  uViol -= demandCutSetL(sID,data_S);
		  return uViol;
		}
////////////////////
		double demandCutSetL(int sID, Data_S *data_S)
		{
		  double sumVal=0;
		  for(int k=0; k<commSet.size(); k++)
		    if( (baseNodes[0]== data_S->getOd(commSet[k],0) && baseNodes[1]!= data_S->getOd(commSet[k],1)) || (baseNodes[1] == data_S->getOd(commSet[k],0) &&  baseNodes[0]!= data_S->getOd(commSet[k],1))  )
		      sumVal += data_S->getD(commSet[k],sID);
		  return sumVal;	      
		}
/////////////////////
		void liftCutCase1(int arcID, double uViol, double* y_SOL, double **X_Value, Data_S *data_S)
		{
 		  setD1.clear();
 		  setD1.push_back(arcID);
		  for(int a=0; a< cutSet.size(); a++)
 		      if(cutSet[a] != arcID && !(find(setC1.begin(), setC1.end(), cutSet[a]) != setC1.end()) )
 			if( sumX(cutSet[a],X_Value) -  min(-uViol, min(data_S->getU(cutSet[a]), demandSum))*y_SOL[cutSet[a]]  > 0  )
 			  setD1.push_back(cutSet[a]);
		}
///////////////////////
		double liftCoeff(int arcID, double uViol, Data_S *data_S)
		{
		    double Fz=0;		  
		    double z = min(data_S->getU(arcID), demandSum);
		    //
		    vector <pair <double,int>  > aux ;
		    for(int a=0; a< setC2.size(); a++){
		      double aux1 =min(data_S->getU(setC2[a]), demandSum);
		      if(aux1 >  uViol)
			aux.push_back(make_pair(aux1,setC2[a]));
		    }
		    for(int a=0; a< setD1.size(); a++){
		      double aux1 =min(data_S->getU(setD1[a]), demandSum);
		      if(aux1 >  uViol)
			aux.push_back(make_pair(aux1,setD1[a]));
		    }
		    //we sort increasingly but we want the reverse
		    sort(aux.begin(), aux.end(), PairCompare);
		    //get the w_k
		    vector <double> auxW;
		    auxW.push_back(0);
		    int t=0;
		    for(int k=aux.size()-1; k >=0; k--)
		      auxW.push_back(auxW[t++] + aux[k].first);
		    //check the lifting function
		    for(int k= 0; k < auxW.size()-1; k++)
		      if(auxW[k] <= z && z < auxW[k+1]-uViol)
			  return k*uViol;
		    for(int k= 1; k < auxW.size(); k++)
		      if(auxW[k]-uViol <= z && z < auxW[k])
			  return k*uViol-auxW[k] + z;
		    if(auxW[auxW.size()-1] <= z)
		      return (auxW.size()-1)*uViol - auxW[auxW.size()-1] + z;
		    
		    cout << "none of the conditions are satisfied ---- Something is wrong" << endl;
		}
///////////////////////
		bool IsOriginNode(Data_S *data_S, int nodeID)
		{
			bool itIs=false;
			//check if the node is origin of at least one commodity 
			for(int k = 0; k < data_S->getN_od(); k++)
				if(nodeID == data_S->getOd(k,0))
					itIs =true;
			return itIs;
		}
//////////////////////		
		bool IsDestNode(Data_S *data_S, int nodeID)
		{
			bool itIs=false;
			for(int k = 0; k < data_S->getN_od(); k++)
				if(nodeID == data_S->getOd(k,1))
					itIs =true;
			return itIs;
		}
/////////////////////////
		void CutSetGenerator(Data_S *data_S, bool outArc, bool inArc)
		{
			//reinitialize the cutset array
			cutSet.clear();
			complementoryCutSet.clear();
			//find the cut set
			for(int a = 0; a < data_S->getN_arcs()-data_S->getN_od(); a++){//dummies are not part of cutset
				if(outArc){
				  if( (baseNodes[0] == data_S->getArcs(a,0) && baseNodes[1] != data_S->getArcs(a,1)) || (baseNodes[1] == data_S->getArcs(a,0) && baseNodes[0] != data_S->getArcs(a,1)) )//if we want the cutset
					cutSet.push_back(a);
				  else if( (baseNodes[0] == data_S->getArcs(a,1) && baseNodes[1] != data_S->getArcs(a,0)) || (baseNodes[1] == data_S->getArcs(a,1) && baseNodes[0] != data_S->getArcs(a,0)) )
					complementoryCutSet.push_back(a);
				}
				if(inArc){
				  if( (baseNodes[0] == data_S->getArcs(a,1) && baseNodes[1] != data_S->getArcs(a,0) ) || (baseNodes[1] == data_S->getArcs(a,1) && baseNodes[0] != data_S->getArcs(a,0)) )//if we want its complement
					cutSet.push_back(a);
				  else if( (baseNodes[0] == data_S->getArcs(a,0) && baseNodes[1] != data_S->getArcs(a,1) ) || (baseNodes[1] == data_S->getArcs(a,0) && baseNodes[0] != data_S->getArcs(a,1)) )
					complementoryCutSet.push_back(a);
				}
			}
		}
/////////////////////////
		void commoditySet(int n_sc, int arcID, int sID, double **X_Value, Data_S *data_S)
		{
		  commSet.clear();
		  demandSum=0;
		  for(int k=0; k<data_S->getN_od(); k++)
		    if(X_Value[arcID][k] > 1e-17){
		      commSet.push_back(k);
		      demandSum += data_S->getD(k,sID);
		    }
		  for(int s=0; s<n_sc; s++){
			double aux=0.0;
			for(int k=0; k<data_S->getN_od(); k++)
				if(X_Value[arcID][k] > 1e-17)
				  aux += data_S->getD(k,s);				
			scenDemandSum.push_back(aux);
		  }
		  
		}
////////////////////////
		void findSetC1andC2(int arcID, double **X_Value, double* y_SOL, Data_S *data_S)
		{
		  setC1.clear();
		  setC2.clear();
		  for(int aa=0 ; aa < cutSet.size(); aa++)
		    if( cutSet[aa] != arcID && sumX(cutSet[aa],X_Value) > 1e-17 + (1-y_SOL[arcID]) * min(data_S->getU(cutSet[aa]), demandSum) )
		      setC1.push_back(cutSet[aa]);
		  for(int aa=0 ; aa < complementoryCutSet.size(); aa++)		  
		    if( min(data_S->getU(complementoryCutSet[aa]), demandSum) * y_SOL[arcID] + 1e-17 < sumX(complementoryCutSet[aa],X_Value) )
		      setC2.push_back(complementoryCutSet[aa]);
		}
////////////////////////
		double sumX(int arcID, double **X_Value)
		{
		    double X_a_L=0;
		    for(int k=0; k < commSet.size(); k++)
				X_a_L += X_Value[arcID][commSet[k]];
		    return  X_a_L;
		}
//////////////////////
		int getNumFPI()    											{return numFPI;}
		double getFixPart(int i) 									{return fixPart[i];}
		vector <double>& getScenFixPart(int i) 						{return scenFixPart[i];}
		vector <double>& getScenArcCoeff(int i, int s, int n_sc) 	{return scenArcCoeff[i*n_sc + s];}
		vector <double>& getArcCoeff(int i) 						{return arcCoeff[i];}
		vector <double>& getConCoeff(int i) 						{return conCoeff[i];}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////						
		~FlowCuts()
		{			
		}

};
#endif
