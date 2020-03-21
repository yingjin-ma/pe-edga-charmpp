#include <cmath>
#include <iterator>
#include <iostream>
#include <iomanip>
#include <string>
#include <time.h>
#include <sys/stat.h>
#include <chrono>
#include <vector>
#include <boost/lexical_cast.hpp>
#include <boost/tokenizer.hpp>
#include <boost/tuple/tuple.hpp>
#include <ctime>

using std::cerr;
using std::cout;
using std::endl;

#ifdef USE_AMBIENT
    #include "dmrg/block_matrix/detail/ambient.hpp"
    typedef ambient::numeric::tiles<ambient::numeric::matrix<double> > matrix;
#else
    #include "dmrg/block_matrix/detail/alps.hpp"
    typedef alps::numeric::matrix<double> matrix;
#endif

#include "dmrg/utils/DmrgParameters.h"
#include "dmrg/block_matrix/indexing.h"
#include "dmrg/mp_tensors/mps.h"

#include "ci_encode.hpp"
#include "sampling.hpp"

#include "DetCollection.h"

// The Charm++ part
#include "EDGA_charmc.h"
#include "EDGA_charmc.def.h"
#include <vector>
#include <cassert>
#include "searchEngine.h"

template<class Matrix, class SymmGroup>
    std::vector<std::vector<typename SymmGroup::charge> >
    parse_config(std::string file, std::vector<Index<SymmGroup> > const & site_dims, std::vector<int> order)
   {
    std::ifstream config_file;
    config_file.open(file.c_str());

    typename SymmGroup::charge A(1), B(0), C(0), D(0);
    B[0]=1; C[1]=1;

    std::vector<std::vector<typename SymmGroup::charge> > configs;

    for (std::string line; std::getline(config_file, line); ) {
        std::vector<std::string> det_coeff;
        boost::split(det_coeff, line, boost::is_any_of(" "));

        std::string det = det_coeff[0];

        if (det.size() != site_dims.size())
          throw std::runtime_error("The determinant length doesn't match the mps length\n");

        std::vector<typename SymmGroup::charge> tmp;
        for (std::size_t i = 0; i < det.size(); ++i) {
            int occ = boost::lexical_cast<int>(det[order[i]-1]);
            switch(occ) {
                case 4:
                    tmp.push_back(site_dims[i][0].first); // doubly occ
                    break;
                case 3:
                    tmp.push_back(site_dims[i][1].first); // up
                    break;
                case 2:
                    tmp.push_back(site_dims[i][2].first); // down 
                    break;
                case 1:
                    tmp.push_back(site_dims[i][3].first); // empty
                    break;
            }
        }
        configs.push_back(tmp);
    }
    return configs;
}

Main::Main(CkArgMsg* msg)
	: numRecv(0)
    {
        numQueens = 1;
        initial_grainsize = 1;
	CI_completeness = 0.950;
	CI_threshold    = 0.001;
	mutation_ratio  = 0.350;
        tmp_folder      = "tmp"; 
//        iflag_txt       = 1;	
        int update_interval = 3;
	int iflag_entanglement = 0;
        std::string QuanEnt_file;  
	nloops_main=1;

        clock_t startTime, mid1Time, mid2Time, mid3Time, endTime;
        startTime = clock();

        std::string command;
	command="date";
	int isys1=system(command.c_str());

        if(msg->argc < 2)
        {
	 ckout << " ======================================================================================================"
	       << "\n              The Charm++ version of Determinants Analysis tool for a MPS State "
	       << "\n                         (MPS-to-CI, Monte-Carlo version for testing)"
	       << "\n ----------------------------------------------------------------------------------------------------"
	       << "\n  Usage: charmrun  +p  N  ./EDGA_charmc  <mps.h5>  <determinants_file>  Nchares  "
	       << "\n                          MAX_iter  CI_threshold  Q-Entanglement-file  tmp_Folder  "
	       << "\n ===================================================================================================="
               << endl;
         CkExit();
        }


	// load the MPS
        load(msg->argv[1], mps);

        int L = mps.length();
	nsite_mps=L;

        //loading the archive:
        storage::archive ar_in(msg->argv[1]+std::string("/props.h5"));
        //loading the parameters
        DmrgParameters parms;
        ar_in["/parameters"] >> parms;
        //read ordering from parms:
        std::vector<int> order(mps.length());
        order = parms["orbital_order"].as<std::vector<int> >();

        maquis::cout << " orbital_order " << std::endl;
        for (int p = 0; p < L; ++p)
            {
            //order[p] = p + 1;
            maquis::cout << " " << order[p] << " "; // Here should be mps[order[p]]
            }
        std::cout << std::endl;

        // extract physical basis for every site from MPS
        std::vector<Index<TwoU1PG> > phys_dims;
        for (int p = 0; p < L; ++p)
            phys_dims.push_back(mps[p].site_dim());

        // load the determinants
        std::vector<std::vector<typename TwoU1PG::charge> > determinants = parse_config<matrix, TwoU1PG>(std::string(msg->argv[2]), phys_dims, order);

	// printout few determinants to check
        for (int q = 0; q < 1; ++q){
        for (int p = 0; p < L; ++p){
            std::cout << determinants[q][p];
            std::cout << std::endl;
            }
          std::cout << std::endl;
        }

        // Get the number of electrons & holes (alpha, beta, total)
        nele_alpha=0;
        nele_beta=0;
        nhole_alpha=0;
        nhole_beta=0;
        nele_total=0;
        igroup_sym=0;
        for(int p = L-1; p >= 0; --p)
           {
            if(determinants[0][p][0]==1)
              nele_alpha++ ;
            else
              nhole_alpha++;
            if(determinants[0][p][1]==1)
              nele_beta++  ;
            else
              nhole_beta++ ;
            }
        nele_total=nele_alpha+nele_beta;
        maquis::cout << " nele_alpha:" << nele_alpha << " nele_beta:" << nele_beta << " nhole_alpha:" << nhole_alpha << " nhole_beta:" << nhole_beta << std::endl;

	// Nums. of dets
        int ndets=determinants.size();
	numDets=ndets; 
        std::cout << " ndets : " << numDets << std::endl;

        // The number of chares in Charm++
        if(msg->argc > 3){
          initial_grainsize = atof(msg->argv[3]);}
        if(msg->argc > 4){
          nloops_main       = atof(msg->argv[4]);}
        if(msg->argc > 5){
          CI_threshold      = atof(msg->argv[5]);}
        if(msg->argc > 6){	  
	  iflag_entanglement=1;	
          QuanEnt_file      = std::string(msg->argv[6]);}
        if(msg->argc > 7){	  
          tmp_folder        = std::string(msg->argv[7]);}
//        if(msg->argc > 8){
//	  std::string txt_bin = std::string(msg->argv[8]);
//	  if(txt_bin=="bin")
//	    iflag_txt=0;	   
//	}
	 
	 tmp_folder        = tmp_folder+"/";
	 command = "rm -r " + tmp_folder;
	 int isys2=system(command.c_str()); 
	 command = "mkdir -p " + tmp_folder;
	 int isys3=system(command.c_str()); 

	//If use the entanglement
        std::vector<std::vector<double>> EntMat;
	if(iflag_entanglement==1){ 
          std::cout << " QuanEnt_file : " << QuanEnt_file << std::endl;
	  EntMat.resize(nsite_mps);
          for (int p = 0; p < L; p++){  
            EntMat[p].resize(nsite_mps);
          }	  
          get_entanglement(QuanEnt_file, EntMat, order);
          for(int i = 0; i < nsite_mps; i++){
            for(int j = 0; j < nsite_mps; j++){
              cout.setf(std::ios::fixed);
              std::cout  << std::setprecision(2) << EntMat[i][j] << " " ;}
            std::cout << " " << std::endl ;	  
	  }
	}  

        // Extract the needed data in order to rewrite the extract-coefficient function
        MPSCharm chareMPS=mps;
        extract_coefficient_pre<matrix, TwoU1PG>(chareMPS);

        // Result for saving
	std::vector<ResultCharm> doneCharm;
        doneCharm.resize(1); 
        double sum_ci2=0.0;

        mid1Time = clock();
	std::cout << "The mid1time is: " <<(double)(mid1Time - startTime) / CLOCKS_PER_SEC << "s" << endl;

	// MPS2CI functions 
        unsigned int nchares=initial_grainsize;
	unsigned int nlop=numDets/initial_grainsize;
	unsigned int nmod=numDets%initial_grainsize;
        std::cout << " nlop : " << nlop << " nmod : " << nmod << std::endl;
        std::vector<std::vector<ChargeCharm>> determinants_charm; 

	// Scratch 
//        std::ofstream fs;
//        fs.open("DETs.scr");

	// Initial Dets 
        determinants_charm.resize(numDets);
        for (unsigned int q = 0; q < numDets; q++){
          determinants_charm[q].resize(L);
          for (int p = 0; p < L; p++){        
            determinants_charm[q][p][0]=determinants[q][p][0];
            determinants_charm[q][p][1]=determinants[q][p][1];
            determinants_charm[q][p][2]=determinants[q][p][2];
          }
        }

	// Start the cycles
	for(int icycle=0;icycle<nloops_main+1;icycle++){
		
          if(icycle==0){
            std::cout << " cycle - " << icycle << std::endl;
            mps2ci_in_loop(determinants_charm,chareMPS,nchares,doneCharm,nlop,nmod,icycle);
	  }
          else{
            std::cout << " cycle - " << icycle << std::endl;
	    if(iflag_entanglement==1){ 
              EntGenetic_in_loop(determinants_charm,chareMPS,nchares,doneCharm,nlop,nmod,icycle,EntMat);
            }
	    else{
            //MonteCarlo_in_loop(determinants_charm,chareMPS,nchares,doneCharm,nlop,nmod,icycle);
            //Genetic_in_loop(determinants_charm,chareMPS,nchares,doneCharm,nlop,nmod,icycle);
             }
	  } 
	}

//        std::vector<std::vector<ChargeCharm>> determinants_reservoir;
	// recorded dets 
//        determinants_reservoir.resize(1);
//        for (int q = 0; q < ndets; ++q){
//           determinants_reservoir[q].resize(L);
//	 }   


//        mid2Time = clock();
//  	std::cout << "The mid2time (after Pre) is: " <<(double)(mid2Time - startTime) / CLOCKS_PER_SEC << "s" << endl;
//
//        extract_coefficient_det<matrix, TwoU1PG>(mps,determinants_charm[0]);
        CkPrintf("\n >>>>> grep : need make_left_paried ??? <<<<< \n ");   
//

//        maquis::cout << " ==== Start the Monte-Carlo/Genetic iterations ==== " << std::endl;     
//        boost::random::mt19937 generator;
//        generator.seed(123456+time(NULL));


/*	
        mid3Time = clock();
	std::cout << "The mid3time (before charm++) is: " <<(double)(mid3Time - startTime) / CLOCKS_PER_SEC << "s" << endl;

        mps2ci_charm(determinants_charm,chareMPS,doneCharm);

        maquis::cout << " ==== after mps2ci_charm ==== " << std::endl;

        endTime = clock();
	std::cout << "The end-time is: " <<(double)(endTime - startTime) / CLOCKS_PER_SEC << "s" << endl;

//        CkPrintf(" >>>>> DetCharm[%d].coeff : %f \n ", 0, DetCharm[0].coeff);   
//        CkPrintf(" >>>>> DetCharm[%d].coeff : %f \n ", 1, DetCharm[1].coeff);   

//	CkCallback ckgetCIs(CkIndex_DetCharmClass::cktest(), mainProxy);

//        DetCharm.getCIs(determinants_charm,chareMPS,doneCharm);
       
*/

/*	
	unsigned int iidet=0;
//        std::vector< std::vector<ChargeCharm> >::iterator idet;
	std::vector<ResultCharm>::iterator it;
	for(it=doneCharm.begin();it!=doneCharm.end();it++){
          CkPrintf(" %d-th det : %f  \n", iidet++, it->val);
          for(int j=0; j<2; j++){
            for(int i=0; i<L; i++){
//               CkPrintf(" %d ",it->det[i][j]);
              } 
//              CkPrintf("\n");
            }
          CkPrintf("\n");
	}
	

        if(sum_ci2>CI_completeness){
	  maquis::cout << " ==== Satisfy the CI_completeness, job finish ==== " << std::endl;
	}
	else{	  
          maquis::cout << " ==== Start the Monte-Carlo/Genetic iterations ==== " << std::endl;	   
          boost::random::mt19937 generator;
	  generator.seed(123456+time(NULL));

          // Get the number of electrons & holes (alpha, beta, total)
          nele_alpha=0;
          nele_beta=0;
          nhole_alpha=0;
          nhole_beta=0;
          nele_total=0;
          igroup_sym=0;
          for(int p = L-1; p >= 0; --p)
             {
              if(determinants_charm[0][p][0]==1)
                nele_alpha++ ;
              else
                nhole_alpha++;
              if(determinants_charm[0][p][1]==1)
                nele_beta++  ;
              else
                nhole_beta++ ;
              }
          nele_total=nele_alpha+nele_beta;
	  maquis::cout << " nele_alpha:" << nele_alpha << " nele_beta:" << nele_beta << " nhole_alpha:" << nhole_alpha << " nhole_beta:" << nhole_beta << std::endl;

          boost::uniform_int<> real(1,nele_total);
          int nele_excited;

          // The Monte-Carlo iterations
	  int iter=0; 
	  do {

            nele_excited = real(generator);
	    
            iter=iter+1;
            maquis::cout << " iter-" << iter << " with " << nele_excited << " excited electrons" << std::endl; 	    
	    //DetCharm.monte_carlo(determinants_charm, chareMPS, determinants_reservoir, nele_excited);

	  } while(sum_ci2<CI_completeness && iter<5);
	  maquis::cout << " ==== Done the Monte-Carlo/Genetic iterations with final iter=" << iter <<  " and CI_completeness=" << sum_ci2 << std::endl;
	}	
*/	

/*
        maquis::cout << " ==== debugging ==== " << std::endl; 
//        std::cout << " det_charge_charm[0] in C : "  <<  det_charge_charm[0]  << std::endl; 
//        maquis::cout << "   CI coefficient:   " << extract_coefficient2<matrix, TwoU1PG>(mps, det_charge_charm) << std::endl; 
        maquis::cout << " =================== " << std::endl; 
*/
        
// old version for checking    
        if(msg->argc == 100){
           //compute the CI coefficients for all determinants in the input
           for (typename std::vector< std::vector<TwoU1PG::charge> >::iterator it = determinants.begin(); it != determinants.end(); ++it)
              maquis::cout << "CI coefficient: " << extract_coefficient(mps, *it) << std::endl;
            }

        searchEngineProxy.start();

//    cout << "==========================================================================\n"
//        << "------------------->     The Charm++ EDGA Project     <--------------------\n"
//        << "--> Main Authors: Yingjin Ma, Ting Wang\n"
//        << "-->   Co-Authors: Lian Zhao, Jinrong Jiang \n"
//        << "==========================================================================" << endl;
/*
    if(msg->argc < 4)
    {
        ckout << "USAGE: " << msg->argv[0] << " Folder(mps.h5) [options]" << endl;
        ckout << "options:\n"
            << " -i Initialization_Parameters<string, filename>\n"
            << " -s SeedFile<string, filename>\n"
            << " -e Mutualfile<string, filename>\n"
            << " -p PopulationSize<int>, default = 2000\n"
            << " -c Crossover_Rate<double>, default = 0.5\n"
            << " -m Mutation_Rate<double>, default = 0.3\n"
            << " -t Cutoff_Threshold<double>, default = 1e-6\n"
            << " -f Target_Completeness<double>, default = 0.9999\n"
            << " -l #Loop<int>, default = 10000\n"
            << " -n #threads<int>, default = 1 (serial). use 0 for hardware_concurrency()\n"
            << " -b Ncycles_to_check_list<int>, default = 50"
            << endl;
    }
*/

        delete msg;

// finished();
}

//#include "EDGA_charmc.def.h"
