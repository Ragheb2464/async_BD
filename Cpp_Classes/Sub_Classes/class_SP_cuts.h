#ifndef CLASS_SP_CUT_H
#define CLASS_SP_CUT_H


#include "./class_flow_cuts.h"
ILOSTLBEGIN



class SubproblemCuts
{ 
	private:
		int numResCuts;
		vector <double> resFixCoeff;
		vector< vector<double > >  resArcCoeff, scenResArcCoeff, scenResFixCoeff;
		vector< vector<double > >  resConCoeff;
		FlowCuts *flowCut;
	public:
		SubproblemCuts() 
		{	
		  numResCuts=0;
		  flowCut = new FlowCuts();
		}
//////////////////////////////////
		void genResidualCuts(double* y_SOL, double **X_Value,  Data_S* data_S, int actual_sc)
		{
		  vector <int> commSet;
		  vector <double> dumyVec, dumyVec1;
		  numResCuts=0;
		  double aux, viol;
		  for(int a=0; a <  data_S->getN_arcs()-data_S->getN_od(); a++) 
		      if (data_S->getArcs(a,1) != data_S->getArcs(a,0) && (y_SOL[a]<1 && y_SOL[a]>0) ){
			  //finding the commodity set
			  commSet.clear();
			  for(int k = 0; k < data_S->getN_od(); k++)
			      if( X_Value[a][k] > data_S->getD(k,actual_sc)* y_SOL[a]  )
				commSet.push_back(k);
			  //checking for violation
			  if(commSet.size() > 0){
			      //condition
			      aux=0;
			      for(int p=0; p< commSet.size(); p++)
				  aux += data_S->getD(commSet[p], actual_sc);	
			      double ratio =  aux/data_S->getU(a);
			      if(1e-2 < ratio && ratio < 1-1e-2){				  
				  viol=0;				  
				  viol += -ratio + (ratio-floor(ratio))*ceil(ratio);
				  //ceoff of the y variables
				  dumyVec1.clear();
				  for(int aa=0; aa <  data_S->getN_arcs(); aa++){
				      if (aa != a)
					    dumyVec1.push_back(0);
				      else{
					    dumyVec1.push_back(-ratio + floor(ratio));
					    viol += (-ratio + floor(ratio))*y_SOL[a];
				      }
				  }				  
				  //constraint coeff
				  dumyVec.clear();
				  for ( int aaa = 0; aaa < data_S->getN_arcs() ; aaa++)
				    if (data_S->getArcs(aaa,1) != data_S->getArcs(aaa,0))//if it is not looping arc
					for ( int k = 0; k < data_S->getN_od(); k++)
					    if(aaa == a && (find(commSet.begin(), commSet.end(), k) != commSet.end()) ){
						dumyVec.push_back(-1/data_S->getU(a));
						viol -= (-1/data_S->getU(a)) * X_Value[a][k];
					    }
					    else
						dumyVec.push_back(0);
				  //record the cut if enough violated
				  if(viol>1e-2){
				      numResCuts++;
				      resFixCoeff.push_back(-ratio + (ratio-floor(ratio))*ceil(ratio) );
				      resArcCoeff.push_back(dumyVec1);
				      resConCoeff.push_back(dumyVec);
				      cout << "amount of violation for RCI is : " << viol << endl;
				  }
			      }
			  }
		      }
		}
//////////////////////////////////
		void genFPI(int n_sc, Data_S* data_S, double* y_SOL, double **X_Value, int actual_sc)
		{ 
		    flowCut->flowPack(n_sc, data_S, y_SOL, X_Value, actual_sc);
			// if(flowCut->getNumFPI()>0)
				// cout << "num FPI: " << flowCut->getNumFPI() << " for scenario: " << actual_sc << endl;
 		    for(int i=0; i < flowCut->getNumFPI(); i++){
 			  numResCuts++;
 			  resFixCoeff.push_back(flowCut->getFixPart(i));
 			  resArcCoeff.push_back(flowCut->getArcCoeff(i));
 			  resConCoeff.push_back(flowCut->getConCoeff(i));
 			  scenResFixCoeff.push_back(flowCut->getScenFixPart(i));
			  for(int s=0; s<n_sc; s++)
				scenResArcCoeff.push_back(flowCut->getScenArcCoeff(i, s, n_sc));
 		    }
		}
/////////////////////////////////		
		int getNumNewCuts() 										{return numResCuts;}
		double getFixCoeff(int i) 									{return resFixCoeff[i];}		
		double getScenFixPart(int i, int s) 						{return scenResFixCoeff[i][s];}
		double getScenArcCoeff(int i, int s, int a, int n_sc) 		{return scenResArcCoeff[i*n_sc + s][a];}	
		double getArcCoeff(int i, int a) 							{return resArcCoeff[i][a];}
		double getConCoeff(int i, int t) 							{return resConCoeff[i][t];}
		void   cleanCutMemory() 									{numResCuts=0;}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////						
		~SubproblemCuts()
		{			
		}

};
#endif
