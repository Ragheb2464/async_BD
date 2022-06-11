#ifndef CLASS_DATA_S1_H
#define CLASS_DATA_S1_H



ILOSTLBEGIN
using namespace std;
#include "./class_global_scenarios.h"

class Data_S 
{
	private: 
	  int n_nodes, n_arcs, n_od, **arcs, **od; 	
	  double *f, *p, *u, *c, **d ;
	  string ins_name;
	  vector<string> *test_data; 		// reading new data into a vector
	  double **ins_scenario;			//used to read data
	  bool DataR, DataS;
	  string InsClass;
	  int maxDemandSceID;
	  double *artificialScenario, *alphaWeights;
	public:		
	  Data_S( string insClass)
	  {
		  InsClass = insClass;
	    test_data     = new vector<string>();
	    DataR=false;
	    DataS=false;
	    if( insClass == "S" )
	      DataS =true;
	    if( insClass == "R" )
	      DataR=true;
	  }
/////////////////////////////////Info Function: to read the data ////////////////////////////////////////////////////////////////////////////////		
	  void info(string instName)
	  {
//            		  const char *instancefilepath ="./Data/instance.dat"; 		// the python file copy the appropirate instance into this folder; so its content changes each time
            		  const char *instancefilepath =instName.c_str(); 		// script
		  test_data->clear();		
		  //reading data into the test_data vector
		  fstream Data_file;
		  Data_file.open(instancefilepath, ios::in);
		  
		  string value;
		  while( Data_file >> value )		//it starts from the begining of the file and starts reading data ... value corresponds to something without space 
		    test_data->push_back( value );		  
		  Data_file.close();
		  
		  if(DataS){
		    n_od       = atoi((*test_data)[17].c_str());
		    ins_name   = (*test_data)[2];		   	// instance name	
		    n_nodes    = atoi((*test_data)[5].c_str()); 	 //node		 
		    n_arcs     = n_nodes * n_nodes + n_od;	//arcs (S is a complete network and we have the dummy arcs)				    
		  }
		  if(DataR){		    
		    ins_name   = (*test_data)[0];		   	// instance name	
		    n_nodes    = atoi((*test_data)[1].c_str()); 	 //node	
		    n_od       = atoi((*test_data)[3].c_str());	//num commodities		  	 
		    n_arcs     = atoi((*test_data)[2].c_str()) + n_od;	//arcs 	+ dummies
		  }
	  }		
/////////////////////////////////Info Function: to read the scenarios //////////////////////////////////////////////////////////////////////////////////
	  void info_scenario(Search_Param *search_param, int n_sc,  string scenarioName)
	  {
	    
//            		      const char *scenariofilepath ="./Data/scenario.dat";  //for the python
            		    const char *scenariofilepath =scenarioName.c_str();    //the script


		    ins_scenario = new double*[n_sc];
		    for(int s = 0; s < n_sc; s++)
			ins_scenario[s]= new double[n_od];
		  
		    int  v_count;
		    /////////////////////////////////////////////
		     p = new double[n_sc];
		   /////////////////////////////////////////
		    vector<string> *scenario_data; 
		    
		    scenario_data = new vector<string>();	
		    scenario_data->clear();
		    
// 		    n_od	   = atoi((*test_data)[17].c_str());   

		    //------------------------------------------------------
		    fstream Data_file;
		    Data_file.open(scenariofilepath, ios::in);
		    
		    string value;
		    while( Data_file >> value )
		      scenario_data->push_back( value );
		    Data_file.close();
		    
		    //------------------------------------------------------
		    double maxD=0, auxMaxD;//to drag out the scenario with largest demand
		    if(DataS){
			v_count = n_od + 1;
			for(int s = 0; s < n_sc; s++){
			    auxMaxD=0;
			    p[s] = fabs(atof((*scenario_data)[v_count].c_str()));
			    for(int k = 0; k < n_od; k++){
					ins_scenario[s][k] = fabs(atof((*scenario_data)[v_count+k+1].c_str()));//sometimes we have commodities in the data file with demand of, e.g., -2e-9					
					if(ins_scenario[s][k] < 1e-6)
						ins_scenario[s][k]=0.0;
					auxMaxD += ins_scenario[s][k];
			    }
			    v_count = v_count + n_od+1;
			    if(auxMaxD > maxD){
					maxD = auxMaxD;
					maxDemandSceID=s;
			    }
			} 
		    }
		    if(DataR){
  			// cout << "number of scenarios " << atoi((*scenario_data)[0].c_str()) << " and " << n_sc << endl;
			n_sc=atoi((*scenario_data)[0].c_str());
			v_count =  1;
			for(int s = 0; s < n_sc; s++){
			    auxMaxD=0;
			    p[s] = fabs(atof((*scenario_data)[v_count].c_str()));
			    for(int k = 0; k < n_od; k++){
			      ins_scenario[s][k] = fabs(atof((*scenario_data)[v_count+k+1].c_str()));//sometimes we have commodities in the data file with demand of, e.g., -2e-9					
				  if(ins_scenario[s][k] < 1e-6)
					ins_scenario[s][k]=0.0;
			      auxMaxD += ins_scenario[s][k];
			    }
			    v_count = v_count + n_od + 1;
			    if(auxMaxD > maxD){
				  maxD = auxMaxD;
				  maxDemandSceID=s;
				  
			    }
			} 			
		   }
		 // ReArrangeSce(n_sc, n_od);
		 SceCreRetention(search_param, n_sc);
		  delete scenario_data;
	  }
/////////////////////////////////Assign Function://////////////////////////////////////////////////////////////////////////////////
	  void read_assign(int n_sc)
	  {	
		 
		      int *origin, *destination, v_count,/* n_od,*/ t;
		      
		      
		     // n_od	   = atoi((*test_data)[17].c_str()); //number O-D

		      //////////////////////////////////////////////////////////////
		     //here I can define the size of these arrays I have defined as pointer in the main function because I have passed them as reference  and also I can fill them and still have them in the main
		      arcs = new int*[n_arcs];//I defined this array here because for the S data set we only need the size of 2 while for the R it is 4
		      for(int i = 0; i < n_arcs; i++)
			  arcs[i] = new int[2];		
		      
		      od = new int*[n_od];
		      for(int i = 0; i < n_od; i++)
			od[i] = new int[2];
		     
		      f = new double[n_arcs];
		      c = new double[n_arcs];//for S data set the cost for commodities deos not depend on the commdity itself

		      d = new double*[n_od];
		      for(int k = 0; k < n_od; k++)
			      d[k] = new double[n_sc];

		      u = new double[n_arcs];
		   //////////////////////////////////////////////////////////////////

		     origin      = new int[n_od];
		     destination = new int[n_od];	
			  
		   //----------------------------------------
		    if(DataS){
			  v_count = 32 + 2*n_nodes*n_nodes + 4; 
			  
			  for(int k = 0; k < n_od;k++)
			      origin[k]=atoi((*test_data)[v_count+k].c_str());			      			  
			  
			  v_count += n_od + 2;
			
			  for(int k =0; k<n_od;k++)
			      destination[k] = atoi((*test_data)[v_count+k].c_str());
		      //------------------------------------------------------------			
			for(int k = 0; k<n_od;k++)
			    {
			      od[k][0]= origin[k]-1;
			      od[k][1]= destination[k]-1;
			    }
		    // --------------------------------------------			 			
			t =0; //a counter for number of arcs			
			for (int i = 0; i < n_nodes; i++)	//we have a complete graph for S data set
			  for (int j = 0; j < n_nodes; j++)
			  {   
			      arcs[t][0] = i; //arc`s origin
			      arcs[t][1] = j;
			      c[t]       = atof((*test_data)[32 + t].c_str());
			      u[t]       = atof((*test_data)[14].c_str());
			      f[t]       = 100;
			      t++;
			  }
			
			for(int k = 0; k < n_od; k++) //for dummy arcs
			{
			  arcs[t][0]= origin[k] - 1;	//od[k][0];
			  arcs[t][1]= destination[k] -1;
			  c[t]      = atof((*test_data)[11].c_str());//outsourcing cost ... 
			  u[t]	    = 1000; //capacity of outsourcing is large enough 
			  f[t]	    = 0;//if we outsource we only pay for the routing cost ... not the fixed cost
			  t++;
			}
			n_arcs = t; //number of arcs
		   //________________________________________________________________________________________-			
		      for(int k = 0; k < n_od; k++)	
			      for(int s = 0; s < n_sc; s++)
				      d[k][s] = ins_scenario[s][k];			  
		    }
		    if(DataR){
		      
			  v_count =   4+7*(n_arcs-n_od); 			  
			  for(int k = 0; k < n_od;k++){
			      origin[k]		= atoi((*test_data)[v_count].c_str());
			      destination[k]    = atoi((*test_data)[v_count+1].c_str());
			      v_count += 3;
// 			      cout << origin[k] << "  " << destination[k]<< endl;
			  }
			  
		       //------------------------------------------------------------			
			for(int k = 0; k<n_od;k++)
			    {
			      od[k][0]= origin[k]-1;
			      od[k][1]= destination[k]-1;
			    }
		       // --------------------------------------------			 			
			t =0; //a counter for number of arcs
			 v_count =  4;
			for (int a = 0; a < n_arcs - n_od; a++)	//we have a complete graph for S data set			 
			  {   
			      arcs[a][0] = atoi((*test_data)[v_count].c_str())-1; //arc`s origin			      
			      arcs[a][1] = atoi((*test_data)[v_count+1].c_str())-1;
			      c[a]       = atoi((*test_data)[v_count+2].c_str());
			      u[a]       = atoi((*test_data)[v_count+3].c_str());
			      f[a]       = atoi((*test_data)[v_count+4].c_str());	
			      v_count +=7 ;
			  }
			t = n_arcs - n_od;
			for(int k = 0; k < n_od; k++) //for dummy arcs
			{
			  arcs[t][0]= origin[k] - 1;	//od[k][0];
			  arcs[t][1]= destination[k] -1;
			  c[t]      = 1e5;//outsourcing cost ... 
			  u[t]	    = 1e5; //capacity of outsourcing is large enough 
			  f[t]	    = 0;//if we outsource we only pay for the routing cost ... not the fixed cost
			  t++;
			}									
		//________________________________________________________________________________________-			
		      for(int k = 0; k < n_od; k++)	
			      for(int s = 0; s < n_sc; s++)
				      d[k][s] = ins_scenario[s][k];			  			
		    }
		    
			  

		      delete [] origin;
		      delete [] destination;		      
	   }
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	  void clean_data(int n_sc)
	  {   			
	    if(DataS){
			vector<vector<int> > New_od(2, vector<int>(n_od,0));
			      
			vector<vector<double> > New_d(n_od, vector<double>(n_sc,0));
			
		        int new_n_od; 
			
			
			for (int k =0; k <n_od; k++)
			  {
			    New_od[0][k] = od[k][0];		     
			    New_od[1][k] = od[k][1];
			  }
			    
			  for(int k = 0; k < n_od; k++)	
				for(int s = 0; s < n_sc; s++)
				      New_d[k][s] = d[k][s];
				
			 
			  bool Found_Duplication = false;
			  bool While_Loop =false;
				
			  while( !While_Loop )
			  {
			      for (int k =0; k <(New_od)[0].size(); k++)
			      {
				  for (int kk =k+1; kk <(New_od[0]).size(); kk++)
				    {
				      if ((New_od)[0][k] == (New_od)[0][kk])
					if ((New_od)[1][k] == (New_od)[1][kk])
					{	
					  for(int s = 0; s < n_sc; s++)
					    (New_d)[k][s] += (New_d)[kk][s];
					  
					    New_d.erase( New_d.begin() + kk );
					    New_od[0].erase (New_od[0].begin()+kk);
					    New_od[1].erase (New_od[1].begin()+kk);
					    Found_Duplication =true; 
		  // 			 
					}
				      if(Found_Duplication)
					break;
				  }
				  
				  if(!Found_Duplication)
				    While_Loop = true;
			      }
			  }
			  
			  //remove the commodities with the same O-D node
			  for (int k =0; k <(New_od)[0].size(); k++)
			    if((New_od)[0][k] == (New_od)[1][k]){
// 			      for(int s = 0; s < n_sc; s++)
// 				(New_d)[k][s] -= (New_d)[k][s];
			      cout << (New_od)[0][k] << " - " << (New_od)[1][k] << "  " << (New_d)[k][0] << endl;
 			      New_d.erase( New_d.begin() + k );
 			      New_od[0].erase (New_od[0].begin()+k);
 			      New_od[1].erase (New_od[1].begin()+k);
			    }
			      
// 			 for (int k =0; k <(New_od)[0].size(); k++)
// 			   cout << (New_od)[0][k] << " - " << (New_od)[1][k] << "  " << (New_d)[k][0] << endl;
		    /*
		    cout << endl << "I found duplication1 " <<  New_od[0].size() << endl;
		    cout << endl << "I found duplication1 " <<  New_od[1].size() << endl;
				      
				      
			for (int k =0; k <New_od[0].size(); k++)
			  cout << endl << "Origin of " << k << " -> " << New_od[0][k] << " + " <<  New_od[1][k] << " with demand of: " << New_d[k][0];					
		           */				      
			 if ( New_od[0].size() < n_od  )
			   update_data( n_sc, New_od[0].size(), New_od, New_d);
			 else
			 {
			    New_d.clear();
			    vector<vector<double> >().swap( New_d );	//this is how to delete a vector!!!	      
			    New_od.clear();
			    vector<vector<int> >().swap( New_od );	
			 }   
	    }
	   }
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	  //it will delete the data vectors and fill them with cleaned data
	  void update_data(int n_sc, int new_n_od, vector<vector<int> >& New_od, vector<vector<double> >& New_d)
	  {
		
		if(DataS){
			cout << endl << endl << "$$I found " << n_od - new_n_od << " duplicated commodity$$ " << endl; 
			/////////////demand///////////////////
			  for(int k = 0; k < n_od; k++)
			      delete [] d[k] ;  
			  delete [] d ;
			  
			  d = new double*[new_n_od];
			  for(int k = 0; k < new_n_od; k++)
			      d[k] = new double[n_sc+1];//+1 is to create space for the artificial scenario demands
			    
			for(int k = 0; k < new_n_od; k++)	
			      for(int s = 0; s < n_sc; s++)
				    d[k][s] = New_d[k][s] ;
		      
			New_d.clear();
			vector<vector<double> >().swap( New_d );//this is to free the memory
			
			
					
			///////////OD///////////////////
			for(int i = 0; i < n_od; i++)
			  delete [] od[i];
			delete [] od;
			
			
			od = new int*[new_n_od];
			for(int i = 0; i < new_n_od; i++)
			  od[i] = new int[2];
			
			
		      for (int k =0; k <new_n_od; k++)
			{
			    od[k][0] = New_od[0][k];		     
			    od[k][1] = New_od[1][k];
			}
			  
		      New_od.clear();
		      vector<vector<int> >().swap( New_od );
		      
		      ///////////////////////////////////////////////////
		      for(int i = 0; i < n_arcs; i++)
			  delete [] arcs[i];
		      delete [] arcs;
		      delete [] f; //f= new double[n_arcs];
		      delete [] c;//c = new double[n_arcs];//for S data set the cost for commodities deos not depend on the commdity itself	
		      delete [] u;//
	
		      
		      n_arcs     = n_nodes * n_nodes + new_n_od;
			
		      arcs = new int*[n_arcs];//I defined this array here because for the S data set we only need the size of 2 while for the R it is 4
		      for(int i = 0; i < n_arcs; i++)
			  arcs[i] = new int[2];		
		      f = new double[n_arcs];
		      c = new double[n_arcs];
		      u = new double[n_arcs];
		      
			    
		      int t =0; //a counter for number of arcs
		      
		      for (int i = 0; i < n_nodes; i++)	//we have a complete graph for S data set
			for (int j = 0; j < n_nodes; j++)
			{   
			    arcs[t][0] = i; //arc`s origin
			    arcs[t][1] = j;
			    c[t]       = atof((*test_data)[32 + t].c_str());
			    u[t]       = atof((*test_data)[14].c_str());
			    f[t]       = 100;
			    t++;
			}
		      
		      for(int k = 0; k < new_n_od; k++) //for dummy arcs
		      {
			arcs[t][0]= od[k][0];	//od[k][0];
			arcs[t][1]= od[k][1];
			c[t]      = atof((*test_data)[11].c_str());//outsourcing cost ... 
			u[t]	  = 1000; //capacity of outsourcing is large enough 
			f[t]	  = 0;//if we outsource we only pay for the routing cost ... not the fixed cost
			t++;
		      }
		      
		
		      n_od = new_n_od;
		}
	   }
////////////////////////putting the scenario with largest sum at the end 
	  void  ReArrangeSce(int n_sc, int scID, int i)
	  {
	     double * swap;	     
	     swap 			= ins_scenario[n_sc - 1 - i];
	     ins_scenario[n_sc- 1 - i] 	= ins_scenario[scID];
	     ins_scenario[scID] 	= swap;	
	     double swap1;
	     swap1			= alphaWeights[n_sc - 1 - i];
	     alphaWeights[n_sc - 1 - i]	= alphaWeights[scID];
	     alphaWeights[scID]		=swap1;
	     
	  }
////////////////////
	  void SceCreRetention(Search_Param *search_param, int n_sc)
	  {	      
	  
	      vector <int> gScID;//ID of the global scenarios
	      artificialScenario = new double[n_od];
	      alphaWeights	 = new double[n_sc];
	      GlobalScenarios *sceRetenCreat;
	      sceRetenCreat = new GlobalScenarios();
	      sceRetenCreat->creatModel(search_param, n_od, ins_scenario, n_sc);
	      sceRetenCreat->solveModel(artificialScenario, gScID, alphaWeights);
		  // cout << "RRRRRRRRRR" << gScID.size()<< endl;
	      for(int i=0; i< gScID.size(); i++){
			ReArrangeSce(n_sc, gScID[i], i);
			cout << "real ID of the global scenarios: " <<  gScID[i] << endl;
			if(maxDemandSceID == gScID[i])
			  cout << "*****one of them was max demand scenario******: " << maxDemandSceID << endl;
	      }
	      delete sceRetenCreat;
	  }
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	  void clean_memory(int n_sc)
	  {
	    
		  for(int s = 0; s < n_sc; s++)
		    delete [] ins_scenario[s];
		  delete [] ins_scenario;		  		
       	  		  				  
		  delete test_data;
	  }

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	  int getN_nodes()          {return n_nodes;}
	  
	  int getN_arcs()           {return n_arcs;}
	  
	  int getN_od()             {return n_od;}
	  
	  int getArcs(int i, int j) {return arcs[i][j];}
	  
	  int getOd(int k, int i)   {return od[k][i];}
	  
	  double getF(int a) 	    {return f[a];}
	  
	  double getP(int s)        {return p[s];}
	  
	  double getU(int a)        {return u[a];}
	  
	  double getC(int a)        {return c[a];}
	  
	  double getD(int k, int s) {return d[k][s];}
	  void setD(int k, int s, double value) { d[k][s] = value;}
	  
	  double getArtificialScenario(int k) {return artificialScenario[k];}
	  double getAlphaWeights(int s)	{return alphaWeights[s];}
	  
	  string getInsClass() {return  InsClass;}
	  
	
	  ~Data_S()
	   {
	   }
};

#endif
