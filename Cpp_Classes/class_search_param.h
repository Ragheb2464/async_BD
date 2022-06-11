#ifndef CLASS_Search_Param_H
#define CLASS_Search_Param_H

#include <iostream>
#include <fstream>
#include <stdio.h>
#include <string>


using namespace std;


class Search_Param 
{
	private:	  
		  
		  int n_interations;
		  double tolerance;
		  bool RegSpDuals; 
		  ////~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~		
		  //parallelism related params
		  int 	Tagysol, Tagconti_run, Tagcuts, Tagscnumber, Taglistsc, Tagms_send, Tagms_rec, TagMconti_run, Tagms_ncsend, Tagms_csend, Tagms_run;//these are tag of the messages
		  int 	Tagms_phase, TagNumYsol, Tag_Useless_Sol ;//tag for the message which tells the master to change to the second phase of the McDaniel-Devien	
		  int   solsend_strategy, comm_strategy, cut_strategy, iteration;
		  int 	num_copy_cut_add;
		  double n_cut_wait, initial_core_point;
		  int    num_cut_send; 
		  bool   Data_Clean_UP_Moode, continue_running;		  
		  //global scenarios params
		  int 	master_sc, numArtificialScenarios;
		  //MP valid inequalities
		  int  	LBF_clusters;
		  bool 	LestCapIneq, addCoverIneq, addCardinalityIneq, addRawCardinalityCuts, addCombinatorialCuts, Theta_UB, LB_Lifting, y_Lifting, network_connectivity;
		  //MP features
		  bool 	reducedCostFixing, cleanMaster, looping_arcs, start_from_core_sol;
		  double Reg_Slack, Copy_Slack;			  
		  //Warm start params
		  int 	warmStartIterations;
		  bool 	warmStart;
		  //SP cuts
		  bool 	Strong_Ineq, addResidualCuts, flowPackIneq;
		  //heuristic related params
		  bool 	runHeuristic, runHardHeur, runSoftHeur,  runWholeAlgAsHeuristic;
		  int 	numHeurIter;	  
		  //cut generation related params
		  double mio;
		  bool 	 runModifiedMW, runMW, update_core_point;
		  double timeLim, timeLimI;
		  bool create_art_SPs;
		  bool propagatCut;
		  int numAggrCluster;
		  bool cutAggr;
	public:		
	  
	  Search_Param(int n_sc, int n_workers)
	   {		
			start_from_core_sol	= false;//if true the first solution in the pool will be the initial core point		
			Theta_UB 			= false;//upper bounding on theta if true		
			LB_Lifting 			= false;//true if we want to add the LBF cut
			y_Lifting 			= false;//to add the cut obtained from recourse upper and lower bounding			
			LBF_clusters 		= 1;  //number of clusters in the LBF ineqaulity
			//
			addResidualCuts		= false;//true if add residual cuts to the subproblem, note that when the strongs are added these cuts are not required
			flowPackIneq		= false; //true if your want to add flow pack inequalities
			//
			Strong_Ineq 			= true;//adding strong ineq to the sub-problems and master if true
			network_connectivity 	= true;//to add network connectivity valid inequaliteis  
			LestCapIneq 			= true;//true if we want to add the inequality indicating least capacity that should be opened to cover maximim demand; sum(u_ay_a) >= max_demand
			addCoverIneq 			= true;//true if we want to derive and add cover inequalities
			addCardinalityIneq  	= true; //true if we want to derive and add cardinality inequalities
			addRawCardinalityCuts	= true;//true if we want to add some intial cardinality cuts
			//
			addCombinatorialCuts	= true;//true if we wan to add combinatorial cuts
			reducedCostFixing		= false;// true if you want to add the reduced cost based inequalities (RBI)
			//
			runHeuristic			= false;//true if you want to run the fix and optimize (intermediary) phase
			numHeurIter				= 5; //number of iterations for the heuristics 
			runHardHeur				= true;//true if you want the fix all eligiable variables
			runSoftHeur				= false;//true if you want keep the option to have some variables not fixed and some of them can change their value. 
			runWholeAlgAsHeuristic  = false;//true if you want instead of the second phase only solve the restricted verion of the original DEF formulation
			//			
			update_core_point 	=true; //do not update the core point 
			//**for Papa approach** // the update_core_point=true, runModifiedMW=false, runMW=false, mio=0
			runModifiedMW		=true;//true if you want to run the modified MW subproblem
			runMW				=false;//true if you wan to run the regular Magnanti-Wong subproblem
								//if the first one true and the last two false, then algorithm is based on Papadakos						
			mio = 0.0;// if it is larger than zero, i.e., 1e-7, we generate maximal nondominated cuts in which only one SP will be solved.  //supposedly updating the core point in this method is also beneficial
			//
			warmStartIterations	= 15;
			/*if(!runModifiedMW && !runMW)//if we are using papadakos, no warm start
				warmStart			= false;
			else				
				*/warmStart			= true;//true if we want to deflect the warmStartIterations initiation solutions
			//			
			master_sc 				=0; //number of global scenarios in the partial decomposition
			numArtificialScenarios	=0;//number of artificial scenarios to create 
			create_art_SPs			=false;//true if we create artificial scenarios for each processor 
			propagatCut				=false;//true if you want to propagate the cuts
			numAggrCluster			=min(n_sc,n_sc);//number of clusters to aggregate cuts. 
			cutAggr					=false;
			if(!cutAggr)
				numAggrCluster=n_sc;
			//
			cleanMaster		=false;//true if the clean up procedure is applied
			looping_arcs 	=true;	//to remove arcs that start and end at the same node					
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
			Data_Clean_UP_Moode =true;
			continue_running    =true;		    
			RegSpDuals 			=false;	//if we get the dual solution of regular sub problem, true
			//
			num_cut_send 		=1;		 
			n_interations  		=500;//max number of iterations
			tolerance  			=0.0003;//considered tolerance for stopping		
			initial_core_point 	=0.51; //initial core point (all variables except the dummies take this)			 
			timeLim				=7200;//total run time limit
			timeLimI			=3600;//the time limit of the first phase
					
			solsend_strategy = 1;  // 0 = FIFO,     1 = LIFO, 	2 = Random 	3 = sc ranked 	
								// 4 = last sol then sc ranked
								// 5 = checking last 20 sol and change st. between 0 and 1
								// 6 = arg max sumtheta per fs-solution
								// 7 = sending a convex conbination of those sol that didnt send already
								
								
			comm_strategy    = 1;     //0 = sync (wait for all cuts, solsend_strategy = 1),     1 = async (waits for n_cut_wait)
			
			cut_strategy     = 0;	   //0 = all cuts, 1 =  last sol  2 = sc ranked   3 = random	  
			
			num_copy_cut_add = 0; //we add 2 most violated cuts for now
			
			
			n_cut_wait = ceil(0.4*(n_sc - master_sc));//wait for how many scenarios to return cuts
			if(comm_strategy==0)
			  n_cut_wait = n_sc-master_sc;
					  
			Reg_Slack 	= 150; 
			Copy_Slack 	= 50;
			//communication tags
			Tagysol		=200;  
			Tagysol		=200;  
			Tagconti_run=201;	     
			Tagcuts		=202;	     
			Tagscnumber =203;	     
			Taglistsc 	=204;	     
			Tagms_send 	=300;
			Tagms_rec	=301;	     
			TagMconti_run =302;	     
			Tagms_ncsend =303;	     
			Tagms_csend =304;	     
			Tagms_run	= 305;	
			Tagms_phase =400;	     
			TagNumYsol	=404; 
			Tag_Useless_Sol = 405;
			 
	   }
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	   int getMaster_sc()				{return master_sc;}
	   int setMaster_sc(int s)			{master_sc = s;}
	   
	   int getNum_cut_send()			{return num_cut_send;}
	   void setNum_cut_send(int i)		{ num_cut_send=i;}
	   
	   bool getData_Clean_UP_Moode()	{return  Data_Clean_UP_Moode;}
	   
	   bool getContinue_running()		{return  continue_running;}
	   
	   int getTagysol()				{return Tagysol;}
	   
	   int getTagconti_run()		{return Tagconti_run;}
	   
	   int getTagcuts()				{return Tagcuts;}
	   
	   int getTagscnumber()			{return Tagscnumber;}
	   
	   int getTaglistsc()			{return Taglistsc;}
	   
	   int getTagms_send()			{return Tagms_send;}
	   
	   int getTagms_rec()			{return Tagms_rec;}	   
	   
	   int getTagMconti_run()		{return TagMconti_run;}
	   
	   int getTagms_ncsend()		{return Tagms_ncsend;}
	   
	   int getTagms_csend()			{return Tagms_csend;}
	   
	   int getTagms_run()			{return Tagms_run;}
	   
	   int getTagms_phase()			{return Tagms_phase;}
	   
	   int getTagNumYsol()			{return TagNumYsol;}
	   
	   int getTag_Useless_Sol()		{return Tag_Useless_Sol;}
	   
	   int getN_interations() 		{return n_interations;}
	   
	   double getTolerance() 		{return tolerance;}
	   
	   int getSolsend_strategy()	{return solsend_strategy;}
	   
	   double getN_cut_wait()		{return n_cut_wait;}
	   
	   void setN_cut_wait(double n) {n_cut_wait =n;}
	   
	   const int getNum_copy_cut_add() 	{return num_copy_cut_add;}
	   
	   int getComm_strategy()			{return comm_strategy;}
	   
	   double getInitial_core_point() 	{return initial_core_point;}
	   
	   bool getWroker_start_from_core_sol() {return start_from_core_sol;}
	   
	   double getReg_Slack () 		{return Reg_Slack;}
	   double getCopy_Slack() 		{return Copy_Slack;}
	   	   
	   bool getStrong_Ineq() 		{return Strong_Ineq;}	
	   
	   bool getAddResidualCuts()	{return addResidualCuts;}
	   
	   bool getTeta_UB() 			{return Theta_UB;}
	   
	   bool getLB_Lifting() 		{return LB_Lifting;}
	   
	   bool getY_Lifting() 			{return y_Lifting;}
	   
	   bool getNetwork_Connectivity(){return network_connectivity;}
	   
	   bool getLooping_arcs() 		{return looping_arcs;}
	   
	   bool getUpdateCorePoint() 	{return update_core_point;}
	   void setCorePoint(double c) 	{initial_core_point = c;}
	   
	   bool getRegSpDuals() 		{return RegSpDuals;}
	   
	   int getLBF_clusters() 		{return LBF_clusters;}
	   
	   bool getLestCapIneq()  		{return LestCapIneq;}
	   
	   bool getAddCoverIneq() 		{return addCoverIneq;}
	   
	   bool getAddCardinalityIneq() {return addCardinalityIneq;}
	   
	   bool getAddCombinatorialCuts(){return addCombinatorialCuts;}
	   
	   bool getRunHeuristic() 		 {return runHeuristic;}
	   
	   bool getReducedCostFixing() 	{return reducedCostFixing;}
	   
	   bool getAddRawCardinalityCuts() 	{return addRawCardinalityCuts;}
	   
	   bool getCleanMaster() 		{return cleanMaster;}
	   
	   bool getRunWholeAlgAsHeuristic() {return runWholeAlgAsHeuristic;}
	   bool getRunModifiedMW()			{return runModifiedMW;}
	   bool getRunMW()					{return runMW;}
	   int getWarmStartIterations() 	{return warmStartIterations;}
	   bool getWarmStart() 				{return warmStart;}
	   bool getFlowPackIneq() 			{return flowPackIneq;}
	   int getNumHeurIter() 			{return numHeurIter;}
	   int getNumArtificialScenarios() 	{return numArtificialScenarios;}
	   int setNumArtificialScenarios(int s) 	{ numArtificialScenarios=s;}
	   bool getRunHardHeur() 			{return runHardHeur;}
	   bool getRunSoftHeur() 			{return runSoftHeur;}
	   double getMio()					{return mio;}
	   double getTimeLim()				{return timeLim;}
	   double getTimeLimI()				{return timeLimI;}
	   bool getCreate_art_SPs()			{return create_art_SPs;}
	   bool getPropagatCut()			{return propagatCut;}
	   int getNumAggrCluster() 			{return numAggrCluster;}
	   bool getCutAggr()				{return cutAggr;}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	  ~Search_Param()  {  }
};

#endif
