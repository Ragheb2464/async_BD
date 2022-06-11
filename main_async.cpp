#include <cstring>
#include <sstream>
#include <vector>
#include <iostream>
#include <string>
#include <fstream>
#include <cstdlib>
#include <cmath>
#include <algorithm>
#include <time.h>
#include <ctime>
#include <sys/time.h>
#include <stdio.h>
#include <omp.h>
#include "mpi.h"
#include <pthread.h>

#include "./Cpp_Classes/class_main.h"
#include "./Cpp_Classes/class_sub_problems.h"




#define Comm COMM_WORLD		//means "all the processes in the MPI application.: it name them from 0 to size()-1 I guess
using namespace std;



int main(int argc, char *argv[]) 
{
	
/////////////////////////step 0: parallel env initialization ///////////////////////////////////////////////////////

    
	int n_p, version, subversion, id, namelen, n_workers;
	
	char 	procname[MPI::MAX_PROCESSOR_NAME+1];
	MPI::Init(argc,argv); 		//to enter the MPI world
	n_p = MPI::Comm.Get_size();	//how many processors we have  ... I think in the python file, we have reserved a specific number of processors
	id  = MPI::Comm.Get_rank();	//A process finds out its own rank by calling MPI_Comm_rank(): 
	
	if( !MPI::Is_initialized() ) 	// if MPI is not initialized ... get out ... we have nothing to do
	{
	  cout <<"Error starting MPI program. Terminating.\n";
	  MPI::Comm.Abort(911);
	  return -5;
	}

	MPI::Get_processor_name(procname,namelen);
	procname[namelen] = '\0';
	MPI::Get_version (version, subversion);

	
	n_workers = n_p - 1 ;			//we consider 2 processors for the master (i.e., 0 and n_p-1) and the rest for the workers



//////////////////////////	step 1: Data 	//////////////////////////////////////////////////////////////////////////////	
	if(strcmp(argv[1],"S") != 0 && strcmp(argv[1],"R") != 0)
	{
	    cout<<"\nExit because of No right argaument S or R was given"<<endl;
	    return 0;
	}	


	  int n_sc;
	  std::stringstream charvalue;
	    charvalue << argv[4];  // number of scenario 
	  charvalue  >> n_sc ;				//scenario		  
	  //activate these two whenever you are using the script file
	  string instName = argv[5];	  
	  string scenarioName = argv[6];
	  string insClass = argv[1]; 
//=====================================================================================================================================
		//this is master problem on the processor id=0
//=====================================================================================================================================
	if ( id == 0 ) 	//this is the GSC processor 
	{	
		cout <<"\n\n---------------General Information------------------------------\n";
		cout <<"\nMPI initialized 			"<< MPI::Is_initialized() << "\n";
		cout <<"\nMPI version:                   	"<< version << "." << subversion << endl;
		cout <<"\nNumber of processes:           	"<< n_p<<endl;
		cout <<"\nNumber of Workers: 		        " << n_workers  << "\n";
		cout <<"\nProcesses running on processor: 	"<< procname << endl;
		cout <<"\nFirst instance name:                  " << instName << "\n";
	        cout <<"\nProgram name: 			" << argv[0]  << "\n";
		
		class_main *Main  = new class_main(n_workers, n_sc, instName, scenarioName, insClass);
		Main->main_loop();
		delete Main;
		cout << endl << "end of GSC at processor-" << id << endl;	
      } //end of if (id == 0)
    
//===============================================================================================================================
//the code that we will run on the processor 1 to n_p-2 which solves the sub-problems problem: the workers main loop
//===============================================================================================================================
	if ( id != 0 ) 	//these are the worker ...since they do the same task we have them all here in one if
	{	
		class_sub_problems *SP  = new class_sub_problems( id, n_sc, n_workers, instName, scenarioName, insClass);		
		delete SP;		
		cout<<"\n End of the worker "<< id << endl ;//" with " << iteration << " iterations" << endl;
	}


	MPI::Finalize();
	
	return 0;

}//end main		
