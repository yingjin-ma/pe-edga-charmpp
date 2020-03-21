#include "EDGA_charmc.decl.h"
#include <vector>
#include <cassert>
#include "dmrg/block_matrix/symmetry/nu1pg.h"
#include "CompatibilityWrapper.h"

#include <boost/random.hpp>
#include <iostream>
#include <boost/random.hpp>
#include <math.h>
#include <string.h>
#include <time.h>
#include <chrono>
#include <fstream>

// Used for Charm++ 
/*readonly*/ CProxy_Main mainProxy;
/*readonly*/ int nElements;
int initial_grainsize;
int numQueens;
int numDets;
int nsite_mps;
int nloops_main;

int nele_alpha;
int nele_beta;
int nhole_alpha;
int nhole_beta;
int nele_total;
int igroup_sym;
double CI_completeness;
double CI_threshold;
double mutation_ratio;
//int iflag_txt;
std::string tmp_folder;

MPSCharm mps; 

// Additional parameter in Charm++

// parameters class
//CalcParams Params;


void print_mat(std::vector<std::vector<double>> const & mat, int const & len1, int const & len2){

    for(int i=0; i<len1; i++){
       for(int j=0; j<len2; j++) {
          std::cout << mat[i][j] << " "; 
        }
      std::cout << std::endl; 
     }
};

template <class Matrix, class SymmGroup>
typename Matrix::value_type extract_coefficient2(MPSCharm const & mps, std::vector<ChargeCharm> const & det)
{
    typedef typename SymmGroup::charge charge;

    int idebug=0000;

    block_matrix<Matrix, SymmGroup> const & block0 = mps[0].data();

    if(idebug>0){
      maquis::cout << "\n mps[0].data().basis_.size() : \n " <<  block0.basis_.size() << std::endl;
// check << in dmrg/block_matrix/dual_index.h
      maquis::cout << "\n mps[0].data().basis_ : \n " <<  block0.basis_ << std::endl;
      maquis::cout << "\n mps[0].data().data_.size() : \n " <<  block0.data_.size() << std::endl;
    }
    if(idebug>0){
      for(unsigned p=0; p < block0.data_.size(); p++)  
      {
       maquis::cout << "\n mps[0].data().data_[" << p << "].num_rows() : \n " <<  block0.data_[p].num_rows() << std::endl;
       maquis::cout << "\n mps[0].data().data_[" << p << "].num_cols() : \n " <<  block0.data_[p].num_cols() << std::endl;
       maquis::cout << "\n mps[0].data().data_[" << p << "] : \n " <<  block0.data_[p] << std::endl;
       }
       maquis::cout << "\n mps[0].data() or block0 : \n " <<  mps[0].data() << std::endl;
    }  

      charge sector = det[0];

      std::size_t b = block0.left_basis().position(sector);
  
    if(idebug>0){
// check << in dmrg/block_matrix/indexing_stable.hpp
      maquis::cout << "\n mps[0].data().left_basis() : \n " <<  block0.left_basis() << std::endl;
    }

      if (b == block0.left_basis().size()) return 0.0;

      if(idebug>0){
        std::cout << "\n size_t b : " << b << "   block0.left_basis().size() : " << block0.left_basis().size() << std::endl;
      }

      Matrix coeff = block0[b];

      if(idebug>0){
        maquis::cout << "\n coeff or block0[b] : \n " << coeff << std::endl;
      }
   
      mps[0].make_left_paired();

      if(idebug>0){
        maquis::cout << "\n mps[0].data() after make_left_paired() : \n " <<  mps[0].data() << std::endl;
      }

    for(unsigned p=1; p<mps.length(); ++p)
    {
        Index<SymmGroup> const & phys = mps[p-1].site_dim();

        if(idebug>0){
          maquis::cout << "\n phys : \n " << mps[p-1].site_dim() << std::endl;
        }

        mps[p].make_left_paired();

        if(idebug>0){
          maquis::cout << "\n mps[" << p << "].data() after make_left_paired() : \n " <<  mps[p].data() << std::endl;
        }

        charge site_charge = det[p];

        if(idebug>0){
          maquis::cout << "\n site_charge : \n " << site_charge << std::endl;
        }

        charge left_input = sector;

        if(idebug>0){
          maquis::cout << "\n left_input  : \n " << left_input  << std::endl;
        }

        sector = SymmGroup::fuse(sector, site_charge);

        if(idebug>0){
          maquis::cout << "\n sector \n " << sector  << std::endl;
        }

        block_matrix<Matrix, SymmGroup> const & data = mps[p].data();

        if(idebug>0){
          maquis::cout << "\n data  : \n " << data  << std::endl;
        }

        std::size_t b = data.left_basis().position(sector);
        if(b == data.left_basis().size())
           return 0.0;

        Matrix const & sector_matrix = data[b];

        if(idebug>0){
          maquis::cout << "\n sector_matrix \n " << sector_matrix  << std::endl;
        }

        // determine voffset - row offset for site_charge in sector_matrix
        // if permformance is an issue, do not recalculate the ProductBasis on every call
        ProductBasis<SymmGroup> left_pb(mps[p].site_dim(), mps[p].row_dim());

        if(idebug>0){
          maquis::cout << "\n mps[" << p << "].site_dim() \n " << mps[p].site_dim() << std::endl;
          maquis::cout << "\n mps[" << p << "].row_dim() \n " << mps[p].row_dim() << std::endl;
          maquis::cout << "\n mps[" << p << "].col_dim() \n " << mps[p].col_dim() << std::endl;
        }
//        maquis::cout << "\n left_pb(mps[p].site_dim(), mps[p].row_dim()) \n " << left_pb(mps[p].site_dim(), mps[p].row_dim())  << std::endl;

        std::size_t voffset = left_pb(site_charge, left_input);

        if(idebug>0){
          maquis::cout << "\n voffset \n " << voffset  << std::endl;

          // determine vlen and set hlen, the size of the matrix subblock
          maquis::cout << "\n sector after PB \n " << sector  << std::endl;
          maquis::cout << "\n site_charge after PB \n " << site_charge  << std::endl;
        }
        unsigned vlen = mps[p-1].col_dim().size_of_block(SymmGroup::fuse(sector, -site_charge));
        unsigned hlen = data.right_basis().size_of_block(sector);
        if(idebug>0){
          maquis::cout << "\n vlen & hlen : \n " << vlen << " " <<  hlen  << std::endl;
        }
    
//        CkExit();

        // extract the subblock
        Matrix site_block(vlen, hlen);
        for (unsigned rowcnt = 0; rowcnt < vlen; ++rowcnt){
            std::copy(sector_matrix.row(voffset+rowcnt).first,
                      sector_matrix.row(voffset+rowcnt).second, site_block.row(rowcnt).first);
//            maquis::cout << "\n rowcnt " << rowcnt << "\n "  << std::endl;
//            maquis::cout << "\n site_block in for loop " << site_block  << std::endl;
           }

        if (coeff.num_cols() != site_block.num_rows())
            throw std::runtime_error("coeff.num_cols() != site_block.num_rows()\n");

        if(idebug>0){
          maquis::cout << "\n site_block in loop [" << p << "] : \n  " << site_block  << std::endl;
        }
//        CkExit();

        // multiply subblock with previous subblock
        Matrix tmp(coeff.num_rows(), site_block.num_cols());

        gemm(coeff, site_block, tmp);
        coeff=tmp;
        if(idebug>0){
          maquis::cout << "\n coeff in loop "<< p << "\n " << coeff << " ||  0,0 elem is " << coeff(0,0)  << std::endl;
        }
    }
    maquis::cout << "\n coeff 0,0 elem is " << coeff(0,0)  << std::endl;

    return coeff(0,0);
}

template <class Matrix, class SymmGroup>
typename Matrix::value_type extract_coefficient_det(MPSCharm const & mps, std::vector<ChargeCharm> const & det) {

//   CkPrintf("\n XXXXXX mps.charm_mps_size[0] : %d \n", mps.charm_mps_size[0] );   

   ChargeCharm sector = det[0]; 
   int b=-1; 
   for(int i=0; i< mps.charm_mps_size[0] ; i++){
//      CkPrintf("\n mps.charm_mps[0].basis[%d].lc : ",i);   
//      CkPrintf(" %d %d %d", mps.charm_mps[0].basis[i].lc[0], mps.charm_mps[0].basis[i].lc[1], mps.charm_mps[0].basis[i].lc[2] ); 
//      std::cout << "\n mps.charm_mps[0].basis[i].lc  : <" << mps.charm_mps[0].basis[i].lc[0] << "," << mps.charm_mps[0].basis[i].lc[1] << "," << mps.charm_mps[0].basis[i].lc[2] << ">" << std::endl;
//      std::cout << "\n                       sector  : <" << sector[0] << "," << sector[1] << "," << sector[2] << ">" << std::endl;  
      if(mps.charm_mps[0].basis[i].lc[0]==sector[0])
        if(mps.charm_mps[0].basis[i].lc[1]==sector[1])
          if(mps.charm_mps[0].basis[i].lc[2]==sector[2])
            b=i;
   }
   if (b == -1) return 0.0; 

//   std::cout << "\n after obtain the b value "<< std::endl; 
//   CkPrintf("\n after obtain the b value ");   

   std::vector<std::vector<double>> coeff;
//   CkPrintf("\n ===>  mps.charm_mps[0].data[b].row : %d \n", mps.charm_mps[0].data[b].row);
   coeff.resize(mps.charm_mps[0].data[b].row);  

   for(int i=0; i<mps.charm_mps[0].data[b].row; i++){ 
      coeff[i].resize(mps.charm_mps[0].data[b].col);  
   }     
   coeff=mps.charm_mps[0].data[b].mat;


//   std::cout << "\n coeff : \n "; 
//   print_mat(coeff, mps.charm_mps[0].data[b].row, mps.charm_mps[0].data[b].col);

//   CkPrintf("\n need make_left_paried ???  \n ");   

   for(int i=1; i< mps.charm_nmps ; i++){
//      CkPrintf("\n start the [%d]-th loop \n ", i);   
      std::vector<ChargeCharm> charm_phys;
      charm_phys.resize(mps.charm_mps_size[i-1]);
      for(int j=0; j< mps.charm_mps_size[i-1]; j++){
         charm_phys[j][0]=mps.charm_mps[i-1].basis[j].lc[0];
         charm_phys[j][1]=mps.charm_mps[i-1].basis[j].lc[1];
         charm_phys[j][2]=mps.charm_mps[i-1].basis[j].lc[2];
//         maquis::cout << "\n charm_phys[ " << j << "] : " << charm_phys[j] << " " << std::endl;
         }

//      CkPrintf("\n after charm_phys");   

      ChargeCharm site_charge = det[i]; 
      ChargeCharm left_input = sector;      
      sector[0]=sector[0]+site_charge[0];
      sector[1]=sector[1]+site_charge[1];
      sector[2]=sector[2]+site_charge[2];
      b=-1; 
      for(int ii=0; ii< mps.charm_mps_size[i] ; ii++){
         if(mps.charm_mps[i].basis[ii].lc[0]==sector[0])
           if(mps.charm_mps[i].basis[ii].lc[1]==sector[1])
             if(mps.charm_mps[i].basis[ii].lc[2]==sector[2])
               b=ii;
      }
      if (b == -1) return 0.0; 

//      CkPrintf("\n after site_charge, left_input, sector");          

      std::vector<std::vector<double>> sector_matrix;
      sector_matrix.resize(mps.charm_mps[i].data[b].row);   
      for(int ii=0; ii<mps.charm_mps[i].data[b].row; ii++){ 
         sector_matrix[ii].resize(mps.charm_mps[i].data[b].col);  
      }     
      sector_matrix=mps.charm_mps[i].data[b].mat;

//      CkPrintf("\n after sector_matrix ");          

//      std::cout << "\n sector_matrix : \n "; 
//      print_mat(sector_matrix, mps.charm_mps[i].data[b].row, mps.charm_mps[i].data[b].col);

      int voffset=0;
      for(int i1=0; i1<mps.charm_mps[i].site_dim.size(); i1++){
//         if(i==5){
//            maquis::cout << "\n site_charge " << site_charge; 
//            maquis::cout << "\n sector " << sector;
//            maquis::cout << "\n mps.charm_mps[" << i << "].site_dim " ; 
//            for(int kk=0; kk<mps.charm_mps[i].site_dim.size(); kk++){
//              maquis::cout << "\n  : " << mps.charm_mps[i].site_dim[kk].charge << " " << mps.charm_mps[i].site_dim[kk].ndim  ; 
//            }           
//            maquis::cout << "\n mps.charm_mps[" << i << "].row_dim " ; 
//            for(int kk=0; kk<mps.charm_mps[i].row_dim.size(); kk++){
//              maquis::cout << "\n  : " << mps.charm_mps[i].row_dim[kk].charge << " " << mps.charm_mps[i].row_dim[kk].ndim  ; 
//            }           
//         }
//
         int ii=-1;
         if(mps.charm_mps[i].site_dim[i1].charge[0]==site_charge[0])
           if(mps.charm_mps[i].site_dim[i1].charge[1]==site_charge[1])
             if(mps.charm_mps[i].site_dim[i1].charge[2]==site_charge[2])
               ii=i1;
         if(ii==-1){ 
           for(int i2=0; i2<mps.charm_mps[i].row_dim.size(); i2++){
//
//              if(i==5)
//                maquis::cout << "\n mps.charm_mps[" << i << "].row_dim[" << i2 << "].charge : " << mps.charm_mps[i].row_dim[i2].charge;
//
              if(mps.charm_mps[i].site_dim[i1].charge[0]+mps.charm_mps[i].row_dim[i2].charge[0]==sector[0])                            
                if(mps.charm_mps[i].site_dim[i1].charge[1]+mps.charm_mps[i].row_dim[i2].charge[1]==sector[1])                            
                  if(mps.charm_mps[i].site_dim[i1].charge[2]+mps.charm_mps[i].row_dim[i2].charge[2]==sector[2])
                    voffset=voffset+mps.charm_mps[i].site_dim[i1].ndim*mps.charm_mps[i].row_dim[i2].ndim; 
           }
         } 
         else{
           break;
         }
      }

//      maquis::cout << "\n The voffset : " << voffset << std::endl; 
//      CkPrintf("\n after obtaining voffset ");
       
      int vlen=-1;
      for(int i1=0; i1<mps.charm_mps[i].row_dim.size(); i1++){
         if(mps.charm_mps[i].row_dim[i1].charge[0]==sector[0]-site_charge[0])
           if(mps.charm_mps[i].row_dim[i1].charge[1]==sector[1]-site_charge[1])
             if(mps.charm_mps[i].row_dim[i1].charge[2]==sector[2]-site_charge[2])
               vlen=mps.charm_mps[i].row_dim[i1].ndim;
         if(vlen>0)
           break;            
      }
//      CkPrintf("\n after obtaining vlen ");

      int hlen=-1;
      for(int i1=0; i1<mps.charm_mps[i].col_dim.size(); i1++){
         if(mps.charm_mps[i].col_dim[i1].charge[0]==sector[0])
           if(mps.charm_mps[i].col_dim[i1].charge[1]==sector[1])
             if(mps.charm_mps[i].col_dim[i1].charge[2]==sector[2])
               hlen=mps.charm_mps[i].col_dim[i1].ndim;
         if(hlen>0)
           break;            
      }
//      CkPrintf("\n after obtaining hlen ");
//      maquis::cout << "\n The vlen & hlen : " << vlen << " " << hlen << std::endl; 

      std::vector<std::vector<double>> site_block;
      site_block.resize(vlen);
      for(int i1=0; i1<vlen; i1++){
         site_block[i1].resize(hlen);
      }
//      CkPrintf("\n before std::copy ; the  vlen & hlen : %d  %d with voffset %d : ", vlen, hlen, voffset);      

      for(int rowcnt = 0; rowcnt < vlen; rowcnt++){
          std::copy(sector_matrix[voffset+rowcnt].begin(),sector_matrix[voffset+rowcnt].end(),site_block[rowcnt].begin());
//
//         CkPrintf("\n rowcnt : %d ", rowcnt);      
//         for(int colcnt=0; colcnt < hlen; colcnt++){
//            CkPrintf("\n colcnt : %d ", colcnt);                 
//            CkPrintf("\n sector_matrix[%d][%d] : %f", voffset+rowcnt, colcnt, sector_matrix[voffset+rowcnt][colcnt]);      
////          std::copy(sector_matrix[voffset+rowcnt].begin(),sector_matrix[voffset+rowcnt].end(),site_block[rowcnt].begin());
//            site_block[rowcnt][colcnt]=sector_matrix[voffset+rowcnt][colcnt];
//            CkPrintf("\n site_block[%d][%d] : %f", rowcnt, colcnt, site_block[rowcnt][colcnt]);      
//         }          
//
      }
//      CkPrintf("\n  after std::copy ");      

//      std::cout << "\n site_block : \n "; 
//      print_mat(site_block, vlen, hlen);
      
      std::vector<std::vector<double>> tmp;
      tmp.resize(coeff.size());
      for(int i1=0; i1<coeff.size(); i1++){
        tmp[i1].resize(hlen);
      }

//      CkPrintf("\n after tmp ");      

//
//      cblas_dgemm("CblasColMajor",
//                  "CblasNoTrans",
//                  "CblasNoTrans",
//                  mps.charm_mps[0].data[b].row,
//                  hlen,
//                  mps.charm_mps[0].data[b].col,
//                  1.0,
//                  coeff,
//                  mps.charm_mps[0].data[b].row,
//                  site_block,
//                  mps.charm_mps[0].data[b].col,
//                  0.0,
//                  tmp,
//                  hlen
//                 );
//
//      gemm(coeff,site_block,tmp);
//      dgemm();      
//
             
//      std::cout << "\n DIMENSION[M][N] of      coeff : " << "[" << coeff.size()      << "][" << coeff[0].size()      << "] \n "; 
//      std::cout << "\n DIMENSION[M][N] of site_block : " << "[" << site_block.size() << "][" << site_block[0].size() << "] \n "; 
//      std::cout << "\n DIMENSION[M][N] of        tmp : " << "[" << tmp.size()        << "][" << tmp[0].size()        << "] \n "; 
      
      for(int i1=0; i1<coeff.size(); i1++){ 
//         std::cout << "\n i1 : " << i1 << std::endl; 
         for(int i2=0; i2<hlen; i2++){
//            std::cout << "\n i2 : " << i2 << std::endl; 
            tmp[i1][i2]=0.0;
            for(int i3=0; i3<site_block.size(); i3++){
//               std::cout << "\n i3 : " << i3 << std::endl; 
               tmp[i1][i2]=tmp[i1][i2]+coeff[i1][i3]*site_block[i3][i2];
            } 
         }   
      }
//      CkPrintf("\n after tmp :: second time ");      

//      std::cout << "\n tmp in loop[" << i << "]: \n "; 
//      print_mat(tmp, coeff.size(), hlen);

      for(int i1=0; i1<coeff.size(); i1++){
        coeff[i1].resize(tmp[0].size());
      }
      coeff=tmp;
     
//      CkPrintf("\n Done the loop[%d] with coeff[0][0] is %f ", i, coeff[0][0]);      
//      std::cout << "\n Done the loop[" << i << "]" << "with coeff[0][0] is " << coeff[0][0] <<  std::endl;
   } 

   return coeff[0][0];
}


template <class Matrix, class SymmGroup>
void extract_coefficient_pre(MPSCharm & mps)
{

//   maquis::cout << "\n In extract_coefficient_pre function "<< std::endl; 
    mps.length_charm(mps.length());
//    maquis::cout << "\n in pre : mps.length() : " << mps.length() << std::endl;

    for(int i=0; i<mps.length(); i++)
      mps.size_charm(i,mps[i].data().basis_.size());

//    maquis::cout << "\n in pre : mps.charm_nmps "<< mps.charm_nmps << std::endl; 

    mps.charm_mps_pre();
    for(int i=0; i<mps.length();i++){

//       maquis::cout << "\n mps[" << i << "].site_dim() :\n" << mps[i].site_dim() << std::endl; 
//       maquis::cout << "\n mps[" << i << "].row_dim()  :\n" << mps[i].row_dim() << std::endl; 
//       maquis::cout << "\n mps.charm_mps_size[" << i << "]" << mps.charm_mps_size[i] << std::endl; 
       //.resize(mps.charm_mps_size[i]);

       mps.charm_mps[i].data.resize(mps[i].data().data_.size());
       mps.charm_mps[i].basis.resize(mps[i].data().basis_.size());
       mps.charm_mps[i].site_dim.resize(mps[i].site_dim().size());
       mps.charm_mps[i].row_dim.resize(mps[i].row_dim().size());
       mps.charm_mps[i].col_dim.resize(mps[i].col_dim().size());


//       maquis::cout << "\n mps["<< i << "].data().basis_ : \n " << mps[i].data().basis_ << std::endl;
//       maquis::cout << "\n mps["<< i << "].site_dim().size() : \n " << mps[i].site_dim().size() << std::endl;

       int ii=0; 
       for (typename Index<SymmGroup>::const_iterator it = mps[i].site_dim().begin(); it != mps[i].site_dim().end();  ++it){
//          maquis::cout << " " << "( " << it->first << ": " << it->second << " )" << std::endl;
          mps.charm_mps[i].site_dim[ii].ndim = it->second;
          mps.charm_mps[i].site_dim[ii].charge[0] = it->first[0];
          mps.charm_mps[i].site_dim[ii].charge[1] = it->first[1];
          mps.charm_mps[i].site_dim[ii].charge[2] = it->first[2];
//          maquis::cout << " shadow : " << "( " << mps.charm_mps[i].site_dim[ii].charge << ": " << mps.charm_mps[i].site_dim[ii].ndim << " )" << std::endl;
          ii=ii+1;
       }

       ii=0;
       for (typename Index<SymmGroup>::const_iterator it = mps[i].row_dim().begin(); it != mps[i].row_dim().end();  ++it){        
//          maquis::cout << " " << "( " << it->first << ": " << it->second << " )" << std::endl;
          mps.charm_mps[i].row_dim[ii].ndim = it->second;
          mps.charm_mps[i].row_dim[ii].charge[0] = it->first[0];
          mps.charm_mps[i].row_dim[ii].charge[1] = it->first[1];
          mps.charm_mps[i].row_dim[ii].charge[2] = it->first[2];
//          maquis::cout << " shadow : " << "( " << mps.charm_mps[i].row_dim[ii].charge << ": " << mps.charm_mps[i].row_dim[ii].ndim << " )" << std::endl;
          ii=ii+1;          
       }

       ii=0;
       for (typename Index<SymmGroup>::const_iterator it = mps[i].col_dim().begin(); it != mps[i].col_dim().end();  ++it){        
//          maquis::cout << " " << "( " << it->first << ": " << it->second << " )" << std::endl;
          mps.charm_mps[i].col_dim[ii].ndim = it->second;
          mps.charm_mps[i].col_dim[ii].charge[0] = it->first[0];
          mps.charm_mps[i].col_dim[ii].charge[1] = it->first[1];
          mps.charm_mps[i].col_dim[ii].charge[2] = it->first[2];
//          maquis::cout << " shadow : " << "( " << mps.charm_mps[i].row_dim[ii].charge << ": " << mps.charm_mps[i].row_dim[ii].ndim << " )" << std::endl;
          ii=ii+1;          
       }

       for(int j=0; j<mps[i].data().data_.size(); j++){
          mps.charm_mps[i].data[j].get_rc(mps[i].data().data_[j].num_rows(),mps[i].data().data_[j].num_cols());
//          maquis::cout << "\n mps.charm_mps[" << i << "].site[" << j << "].row : " << mps.charm_mps[i].site[j].row << std::endl;
//          maquis::cout << "\n mps.charm_mps[" << i << "].site[" << j << "].col : " << mps.charm_mps[i].site[j].col << std::endl;
//          Matrix TM1=mps[i].data().data_[j];  
          mps.charm_mps[i].data[j].mat.resize(mps.charm_mps[i].data[j].row);
//          maquis::cout << "\n after mat_resize() row " << std::endl;

          mps.charm_mps[i].basis[j].ls=mps[i].data().basis_[j].ls;
          mps.charm_mps[i].basis[j].rs=mps[i].data().basis_[j].rs;
          mps.charm_mps[i].basis[j].lc[0]=mps[i].data().basis_[j].lc[0];
          mps.charm_mps[i].basis[j].lc[1]=mps[i].data().basis_[j].lc[1];
          mps.charm_mps[i].basis[j].lc[2]=mps[i].data().basis_[j].lc[2];
          mps.charm_mps[i].basis[j].rc[0]=mps[i].data().basis_[j].rc[0];
          mps.charm_mps[i].basis[j].rc[1]=mps[i].data().basis_[j].rc[1];
          mps.charm_mps[i].basis[j].rc[2]=mps[i].data().basis_[j].rc[2];

//          maquis::cout << "\n mps["<< i << "].data().basis_[" << j <<"].lc : " << mps[i].data().basis_[j].lc << std::endl;
//          maquis::cout << "\n mps["<< i << "].data().basis_[" << j <<"].rc : " << mps[i].data().basis_[j].rc << std::endl;
//          maquis::cout << "\n mps.charm_mps[" << i << "].basis[" << j << "].lc : " << mps[i].data().basis_[j].lc << std::endl;          
//          mps[i].data().basis_ << std::endl;
//
          for(int k=0; k<mps.charm_mps[i].data[j].row; k++){
             mps.charm_mps[i].data[j].mat[k].resize(mps.charm_mps[i].data[j].col);
//             maquis::cout << "\n after mat_resize() col " << std::endl;
             for(int l=0; l<mps.charm_mps[i].data[j].col; l++){
//                maquis::cout << " " << mps[i].data().data_[j](k,l);
                mps.charm_mps[i].data[j].mat[k][l]=mps[i].data().data_[j](k,l);
//                maquis::cout << " " << mps.charm_mps[i].data[j].mat[k][l];
               }
//             maquis::cout << std::endl;
            }
         }
      }

    maquis::cout << "\n done the pre   " << std::endl;

//    finished();
//    CkExit();

//    typedef typename SymmGroup::charge charge;

//    block_matrix<Matrix, SymmGroup> const & block0 = mps[0].data();
//    ChargeCharm sector = det[0];

//    BlockMatCharm blockcharm0;


//
//    int L = det.size();
//    for(int i=0; i < L; i++) {
//      CkPrintf(" <<%d", det[i][0]);
//      CkPrintf(",%d",  det[i][1]);
//      CkPrintf(",%d>> ",det[i][2]);
//      }

////    Index<SymmGroup> bindex = block0.left_basis();
//    IndexCharm bindex;// = block0.left_basis();
////    bindex = blockcharm0.left_basis(); 
//    bindex.size(); 
////    std::size_t b = bindex.position(sector);
//
//    Matrix coeff = block0[0];
////    if (b == block0.left_basis().size()) return 0.0;
//
////    Matrix coeff = block0[b];
//
//    return 0.0;
//

}

inline double mps2ci_wrapper(std::vector<ChargeCharm> const & localdet, MPSCharm const & inMPS){
//  return extract_coefficient2<matrix, TwoU1PG>(mps, localdet);
//  return extract_coefficient_test<matrix, TwoU1PG>(mps, localdet);
  return extract_coefficient_det<matrix, TwoU1PG>(inMPS,localdet);
}

class SyncClass : public CBase_SyncClass{

  public:

    SyncClass(CkMigrateMessage *m) {}
    SyncClass(){}

    ~SyncClass(){}    

    void getCIs(double *TL1){
       CkPrintf(" getCIs -- [%d] Proxy = %f \n", thisIndex, TL1[thisIndex]);
     }  

};	

class DetCharmClass : public CBase_DetCharmClass{
  
  public: 

    DetCharmClass(CkMigrateMessage *m) {}
    DetCharmClass(){}

    ~DetCharmClass(){}    

    int Qinf_ito(std::vector<double> EntVec, double dto){

      double entV=0;	     
      for(int i=0; i<EntVec.size(); i++)
	entV=entV+EntVec[i];
      for(int i=0; i<EntVec.size(); i++)
	EntVec[i]=EntVec[i]/entV;

      double cycle[EntVec.size()];
      entV=0;
      for(int i=0; i<EntVec.size(); i++){
	entV=entV+EntVec[i];      
        cycle[i]=entV;
//        std::cout << cycle[i] << " ";	   
      }

//      std::cout << dto << " ";	   
      int ito=0;
      for(int i=0; i<EntVec.size(); i++){
	if(dto<cycle[i]){
	  ito=i;	
	  break;
	}
      }
//      std::cout << ito << std::endl;	   
      return ito; 
    }

    void monte_carlo(const std::vector<std::vector<ChargeCharm>> &inDets, const MPSCharm &inMPS, int nchares, int nlop, int nmod, int icycle){

      int nexcited=0;

      for(int ilop=0; ilop<nlop+1; ilop++)  {

	int ipos=ilop*nchares;

        boost::random::mt19937 generator;
        generator.seed(123456+time(NULL)+this->thisIndex+ilop+icycle);
        boost::uniform_int<> real(1,nele_total);
        nexcited = real(generator);

        boost::uniform_01<> u01;
        boost::uniform_int<> int_alpha_ele(1,nele_alpha);
        boost::uniform_int<> int_alpha_hole(1,nhole_alpha);
        boost::uniform_int<> int_beta_ele(1,nele_beta);
        boost::uniform_int<> int_beta_hole(1,nhole_beta);

// Initial vectors
        std::vector<int>ele_alpha(nele_alpha);
        std::vector<int>ele_beta(nele_beta);
        std::vector<int>hole_alpha(nhole_alpha);
        std::vector<int>hole_beta(nhole_beta);

        if(ilop==nlop){
	  if(thisIndex<nmod) {
            localdet=inDets[thisIndex+ipos];

//            CkPrintf(" thisIndex+ipos -- [%d] \n", thisIndex+ipos);

            int iexcited=0;
            do {

              int i=0,j=0,k=0,l=0;
              for(int p=0; p<nsite_mps; p++ ){
                 if(localdet[p][0]==1)
                    ele_alpha[i++]=p;
                 else
                    hole_alpha[j++]=p;

                 if(localdet[p][1]==1)
                    ele_beta[k++]=p;
                 else
                    hole_beta[l++]=p;
                 }

              iexcited = iexcited + 1;
              if(u01(generator)<0.5){
                int iele_alpha  = int_alpha_ele(generator)-1;
                int ihole_alpha = int_alpha_hole(generator)-1;
                int ifr = ele_alpha[iele_alpha];
                int ito = hole_alpha[ihole_alpha];
                localdet[ifr][0] = 0;
                localdet[ito][0] = 1;
              }
              else{
                int iele_beta   = int_beta_ele(generator)-1;
                int ihole_beta  = int_beta_hole(generator)-1;
                int ifr = ele_beta[iele_beta];
                int ito = hole_beta[ihole_beta];
                localdet[ifr][1] = 0;
                localdet[ito][1] = 1;
              }
            } while (iexcited < nexcited);

            coeff=mps2ci_wrapper(localdet,inMPS);

            if(fabs(coeff)>CI_threshold)
              dets_scratch(ipos,thisIndex,localdet);

	  }
	}
	else { 
          localdet=inDets[thisIndex+ipos];

//          CkPrintf(" thisIndex+ipos -- [%d] \n", thisIndex+ipos);

          int iexcited=0;
          do {

            int i=0,j=0,k=0,l=0;
            for(int p=0; p<nsite_mps; p++ ){
               if(localdet[p][0]==1)
                  ele_alpha[i++]=p;
               else
                  hole_alpha[j++]=p;
  
               if(localdet[p][1]==1)
                  ele_beta[k++]=p;
               else
                  hole_beta[l++]=p;
               }

            iexcited = iexcited + 1;
            if(u01(generator)<0.5){
              int iele_alpha  = int_alpha_ele(generator)-1;
              int ihole_alpha = int_alpha_hole(generator)-1;
              int ifr = ele_alpha[iele_alpha];
              int ito = hole_alpha[ihole_alpha];
              localdet[ifr][0] = 0;
              localdet[ito][0] = 1;
            }
            else{
              int iele_beta   = int_beta_ele(generator)-1;
              int ihole_beta  = int_beta_hole(generator)-1;
              int ifr = ele_beta[iele_beta];
              int ito = hole_beta[ihole_beta];
              localdet[ifr][1] = 0;
              localdet[ito][1] = 1;
            }
          } while (iexcited < nexcited);

          coeff=mps2ci_wrapper(localdet,inMPS);

          if(fabs(coeff)>CI_threshold) 
            dets_scratch(ipos,thisIndex,localdet);
        }
      }
      contribute(CkCallback(CkReductionTarget(Main, finished), mainProxy));
//      CkCallback cb2(CkReductionTarget(Main,report_coeff2),mainProxy);
//      contribute(sizeof(double), &coeff2, CkReduction::sum_double,cb2);

//      monte_carlo2(inDets, inMPS, out); 
    }

    void mps2ci(const std::vector<std::vector<ChargeCharm>> &inDets, const MPSCharm &inMPS, int nchares, int nlop, int nmod)  {

      for(int i=0; i<nlop+1; i++)  {

	int ipos=i*nchares;

        if(i==nlop){
	  if(thisIndex<nmod){ 
            localdet=inDets[thisIndex+ipos]; 
            coeff=mps2ci_wrapper(localdet,inMPS);
            if(fabs(coeff)>CI_threshold) 
              dets_scratch(ipos,thisIndex,localdet);
	  }
	}
        else { 
          localdet=inDets[thisIndex+ipos]; 
          coeff=mps2ci_wrapper(localdet,inMPS);
          if(fabs(coeff)>CI_threshold) 
            dets_scratch(ipos,thisIndex,localdet);
	}

      }
//          CkPrintf(" thisIndex+ipos -- [%d] \n", thisIndex+ipos);
//            CkPrintf(" thisIndex+ipos -- [%d] \n", thisIndex+ipos);
      contribute(CkCallback(CkReductionTarget(Main, finished), mainProxy));
//      CkCallback cb2(CkReductionTarget(Main,report_coeff2),mainProxy);
//      contribute(sizeof(double), &coeff2, CkReduction::sum_double,cb2); 
    }

    void genetic(const std::vector<std::vector<ChargeCharm>> &inDets, const MPSCharm &inMPS, int nchares, int nlop, int nmod, int icycle, unsigned int nitvl)  {

      for(int i=0; i<nlop+1; i++)  {

	int ipos=i*nchares;

        if(i==nlop){
	  if(thisIndex<nmod){ 

            localdet=inDets[thisIndex+ipos]; 

            unsigned int partnerIndex;
            partnerIndex=thisIndex+nitvl+ipos;
            if(partnerIndex>numDets-1)
              partnerIndex=partnerIndex-numDets;

//            CkPrintf(" partnerIndex -- [%d] [%d] \n", thisIndex, partnerIndex);
            partnerdet=inDets[partnerIndex];

            for(int i=0; i<nsite_mps; i++){
              localdet[i][1]=partnerdet[i][1];
            }

            coeff=mps2ci_wrapper(localdet,inMPS);
            if(fabs(coeff)>CI_threshold) 
              dets_scratch(ipos,thisIndex,localdet);
          }
        }
        else { 
          localdet=inDets[thisIndex+ipos]; 

          unsigned int partnerIndex;
          partnerIndex=thisIndex+nitvl+ipos;
          if(partnerIndex>numDets-1)
            partnerIndex=partnerIndex-numDets;

//          CkPrintf(" partnerIndex -- [%d] [%d] \n", thisIndex, partnerIndex);
          partnerdet=inDets[partnerIndex];

          for(int i=0; i<nsite_mps; i++){
            localdet[i][1]=partnerdet[i][1];
          }

          coeff=mps2ci_wrapper(localdet,inMPS);
          if(fabs(coeff)>CI_threshold) 
            dets_scratch(ipos,thisIndex,localdet);
        }
      }
      contribute(CkCallback(CkReductionTarget(Main, finished), mainProxy));
    }

//     
    void Ent_genetic(const std::vector<std::vector<ChargeCharm>> &inDets, const MPSCharm &inMPS, int nchares, int nlop, int nmod, int icycle, unsigned int nitvl, std::vector<std::vector<double>> EntMat)  {
      int nexcited=0;

      for(int ilop=0; ilop<nlop+1; ilop++)  {

	int ipos=ilop*nchares;
	unsigned int seed1;
	unsigned int seed2;

	seed1=123456+time(NULL)+this->thisIndex+ilop+icycle;
	seed2=654321+time(NULL)+this->thisIndex+ilop+icycle;
//            The mutation steps 1/2
        boost::random::mt19937 generator,generator2,gener1,gener2;
        generator.seed(seed1);
        //generator2.seed(seed2);
        //boost::uniform_int<> real(1,nele_total);
        boost::uniform_int<> real(1,nele_total/2);
        nexcited = real(generator);

        boost::uniform_01<> u01;
        boost::uniform_int<> int_alpha_ele(1,nele_alpha);
//        boost::uniform_int<> int_alpha_hole(1,nhole_alpha);
        boost::uniform_int<> int_beta_ele(1,nele_beta);
//        boost::uniform_int<> int_beta_hole(1,nhole_beta);

// Initial vectors
        std::vector<int>ele_alpha(nele_alpha);
        std::vector<int>ele_beta(nele_beta);
        std::vector<int>hole_alpha(nhole_alpha);
        std::vector<int>hole_beta(nhole_beta);


        if(ilop==nlop){
	  if(thisIndex<nmod){ 

            localdet=inDets[thisIndex+ipos]; 

            unsigned int partnerIndex;
            partnerIndex=thisIndex+nitvl+ipos;
            if(partnerIndex>numDets-1)
              partnerIndex=partnerIndex-numDets;

//            CkPrintf(" partnerIndex -- [%d] [%d] \n", thisIndex, partnerIndex);
            partnerdet=inDets[partnerIndex];

            for(int i=0; i<nsite_mps; i++){
              localdet[i][1]=partnerdet[i][1];
            }

//            The mutation steps 2/2
            if(u01(generator) < mutation_ratio){            
              int iexcited=0;
              do {

                int i=0,j=0,k=0,l=0;
                for(int p=0; p<nsite_mps; p++ ){
                   if(localdet[p][0]==1)
                      ele_alpha[i++]=p;
                   else
                      hole_alpha[j++]=p;

                   if(localdet[p][1]==1)
                      ele_beta[k++]=p;
                   else
                      hole_beta[l++]=p;
                   }

                std::vector<std::vector<double>> EntMatA;
                std::vector<std::vector<double>> EntMatB;

		EntMatA.resize(i);
		for(int i1=0; i1<i; i1++){
		  EntMatA[i1].resize(j);	
		  for(int j1=0; j1<j; j1++){
		    EntMatA[i1][j1]=EntMat[ele_alpha[i1]][hole_alpha[j1]];
		    //std::cout << EntMatA[i1][j1] << " " ;
		  }
		  //std::cout << std::endl;
		}
		EntMatB.resize(k);
		for(int k1=0; k1<k; k1++){
		  EntMatB[k1].resize(l);	
		  for(int l1=0; l1<l; l1++){
		    EntMatB[k1][l1]=EntMat[ele_alpha[k1]][hole_alpha[l1]];
		    //std::cout << EntMatB[k1][l1] << " " ;
		  }
		  //std::cout << std::endl;
		}

                iexcited = iexcited + 1;
                generator2.seed(seed2+iexcited);
                if(u01(generator2)<0.5){
                  gener1.seed(seed1+i+j+iexcited);
                  gener2.seed(seed2+i+j+iexcited);
                  int iele_alpha  = int_alpha_ele(gener1)-1;
                  int ifr = ele_alpha[iele_alpha];

		  double dto=u01(gener2);
                  int ihole_alpha= Qinf_ito(EntMatA[iele_alpha],dto);
                  //int ihole_alpha = int_alpha_hole(generator)-1;
                  int ito = hole_alpha[ihole_alpha];
                  localdet[ifr][0] = 0;
                  localdet[ito][0] = 1;
                }
                else{
                  gener1.seed(seed1+k+l+iexcited);
                  gener2.seed(seed2+k+l+iexcited);
                  int iele_beta   = int_beta_ele(gener1)-1;
                  int ifr = ele_beta[iele_beta];

		  double dto=u01(gener2);
		  int ihole_beta = Qinf_ito(EntMatB[iele_beta],dto);
                  //int ihole_beta  = int_beta_hole(generator)-1;
                  int ito = hole_beta[ihole_beta];
                  localdet[ifr][1] = 0;
                  localdet[ito][1] = 1;
                }
              } while (iexcited < nexcited);    
            }	      

            coeff=mps2ci_wrapper(localdet,inMPS);
            if(fabs(coeff)>CI_threshold) 
              dets_scratch(ipos,thisIndex,localdet);
          }
        }
        else { 
          localdet=inDets[thisIndex+ipos]; 

          unsigned int partnerIndex;
          partnerIndex=thisIndex+nitvl+ipos;
          if(partnerIndex>numDets-1)
            partnerIndex=partnerIndex-numDets;

//          CkPrintf(" partnerIndex -- [%d] [%d] \n", thisIndex, partnerIndex);
          partnerdet=inDets[partnerIndex];

          for(int i=0; i<nsite_mps; i++){
            localdet[i][1]=partnerdet[i][1];
          }

//            The mutation steps 2/2
          if(u01(generator) < mutation_ratio){            
  	    int iexcited=0;
            do {

              int i=0,j=0,k=0,l=0;
              for(int p=0; p<nsite_mps; p++ ){
                 if(localdet[p][0]==1)
                    ele_alpha[i++]=p;
                 else
                    hole_alpha[j++]=p;

                 if(localdet[p][1]==1)
                    ele_beta[k++]=p;
                 else
                    hole_beta[l++]=p;
                 }

              std::vector<std::vector<double>> EntMatA;
              std::vector<std::vector<double>> EntMatB;

              EntMatA.resize(i);
              for(int i1=0; i1<i; i1++){
                EntMatA[i1].resize(j);
                for(int j1=0; j1<j; j1++){
                  EntMatA[i1][j1]=EntMat[ele_alpha[i1]][hole_alpha[j1]];
                  //std::cout << EntMatA[i1][j1] << " " ;
                }
                //std::cout << std::endl;
              }
              EntMatB.resize(k);
              for(int k1=0; k1<k; k1++){
                EntMatB[k1].resize(l);
                for(int l1=0; l1<l; l1++){
                  EntMatB[k1][l1]=EntMat[ele_alpha[k1]][hole_alpha[l1]];
                  //std::cout << EntMatB[k1][l1] << " " ;
                }
                //std::cout << std::endl;
              }

              iexcited = iexcited + 1;
              generator2.seed(seed2+iexcited);
              if(u01(generator2)<0.5){
                gener1.seed(seed1+i+j+iexcited);
                gener2.seed(seed2+i+j+iexcited);
                int iele_alpha  = int_alpha_ele(gener1)-1;
                int ifr = ele_alpha[iele_alpha];

                double dto=u01(gener2);
                int ihole_alpha= Qinf_ito(EntMatA[iele_alpha],dto);
                //int ihole_alpha = int_alpha_hole(generator)-1;
                int ito = hole_alpha[ihole_alpha];
                localdet[ifr][0] = 0;
                localdet[ito][0] = 1;
              }
              else{
                gener1.seed(seed1+k+l+iexcited);
                gener2.seed(seed2+k+l+iexcited);
                int iele_beta   = int_beta_ele(gener1)-1;
                int ifr = ele_beta[iele_beta];

                double dto=u01(gener2);
                int ihole_beta = Qinf_ito(EntMatB[iele_beta],dto);
                //int ihole_beta  = int_beta_hole(generator)-1;
                int ito = hole_beta[ihole_beta];
                localdet[ifr][1] = 0;
                localdet[ito][1] = 1;
              }
            } while (iexcited < nexcited);
          }

          coeff=mps2ci_wrapper(localdet,inMPS);
          if(fabs(coeff)>CI_threshold) 
            dets_scratch(ipos,thisIndex,localdet);
        }
      }
      contribute(CkCallback(CkReductionTarget(Main, finished), mainProxy));
    }

//     

// need HDF5 ? 
    void dets_scratch(unsigned int ipos,int index,std::vector<ChargeCharm> det_save){

      int alpha[nsite_mps];
      int  beta[nsite_mps];
      std::string filename;

      for(int i=0; i<nsite_mps; i++){
         alpha[i]=det_save[i][0];
          beta[i]=det_save[i][1];
       } 

      std::string sindex=std::to_string(index);
      filename=tmp_folder+"DETs.scr-"+sindex;
      //std::cout << " iflag_txt = " << iflag_txt << std::endl; 
      std::ofstream fs;
//      if(iflag_txt==1)
      fs.open(filename,std::ios::app);
//      else
//        fs.open(filename,std::ios::binary|std::ios::app);

      for(int i=0; i<nsite_mps; i++){
	 switch(alpha[i]){
           case 1:{
	     switch(beta[i]){
	       case 1: 
		       fs << "4"; 
		       break;
	       case 0: 
		       fs << "3"; 
		       break;
	     }
	     break;
           }
           case 0:{
	     switch(beta[i]){
	       case 1: 
		       fs << "2"; 
		       break;
	       case 0: 
		       fs << "1"; 
		       break;
	     }
	     break;
           }
         }		  
       } 
      fs << " "<< coeff << endl;

     }  

    void cktest(double *TL1){
       CkPrintf(" getCIs -- [%d] Proxy = %f \n", thisIndex, TL1[thisIndex]);
     }

    void report_coeff2(double result)
    {
     if(result>CI_completeness)	{
       CkPrintf("SumCi^2 : %f reach the threshold %f \n", result, CI_completeness);
       CkExit();
       }
     else{
//       monte_carlo_wrapper();        
       CkPrintf("SumCi^2 : %f not reach the threshold %f \n", result, CI_completeness);
       }
     }

    void pup(PUP::er &p){
      p|coeff;
      p|coeff2;
    }    

//  private:
    std::vector<ChargeCharm> localdet;
    std::vector<ChargeCharm> partnerdet;
    MPSCharm localMPS;
    std::vector<ResultCharm> localout;
    int L;
    double coeff=0.0;
    double coeff2=0.0;
    double *pcoeff;
};

class Main : public CBase_Main {

public:
  Main(CkMigrateMessage *m) {}
  Main(CkArgMsg *m);

  int numRecv;
  CProxy_DetCharmClass DetCharm;

  void mps2ci_in_loop(std::vector<std::vector<ChargeCharm>> determinants_charm,  MPSCharm chareMPS, int nchares, std::vector<ResultCharm> doneCharm, int nlop, int nmod, int icycle){

    CkPrintf(" The mps2ci_in_loop     : nchares[%d] nlop[%d] nmod[%d] icycle[%d]\n", nchares, nlop, nmod, icycle);

    DetCharm=CProxy_DetCharmClass::ckNew(nchares);  
    DetCharm.mps2ci(determinants_charm,chareMPS,nchares,nlop,nmod);    
  }

  void MonteCarlo_in_loop(std::vector<std::vector<ChargeCharm>> determinants_charm,  MPSCharm chareMPS, int nchares, std::vector<ResultCharm> doneCharm, int nlop, int nmod, int icycle){

    CkPrintf(" The MonteCarlo_in_loop : nchares[%d] nlop[%d] nmod[%d] icycle[%d]\n", nchares, nlop, nmod, icycle);

    DetCharm=CProxy_DetCharmClass::ckNew(nchares);  
    DetCharm.monte_carlo(determinants_charm,chareMPS,nchares,nlop,nmod,icycle);    
  }

  void Genetic_in_loop(std::vector<std::vector<ChargeCharm>> determinants_charm,  MPSCharm chareMPS, int nchares, std::vector<ResultCharm> doneCharm, int nlop, int nmod, int icycle){

    CkPrintf(" The Genetic_in_loop    : nchares[%d] nlop[%d] nmod[%d] icycle[%d]\n", nchares, nlop, nmod, icycle);

    unsigned int ninterval=0;
    boost::random::mt19937 generator;
    generator.seed(123456+time(NULL)+nchares+icycle);
    boost::uniform_int<> interval(1,numDets-1);
    ninterval = interval(generator);

    CkPrintf(" The GA crossover interval is [%d] with mutation ratio is [] \n", ninterval);

    DetCharm=CProxy_DetCharmClass::ckNew(nchares);  
    DetCharm.genetic(determinants_charm,chareMPS,nchares,nlop,nmod,icycle,ninterval);    
  }

  void EntGenetic_in_loop(std::vector<std::vector<ChargeCharm>> determinants_charm,  MPSCharm chareMPS, int nchares, std::vector<ResultCharm> doneCharm, int nlop, int nmod, int icycle, std::vector<std::vector<double>>& EntMat){

    CkPrintf(" The EntGenetic_in_loop : nchares[%d] nlop[%d] nmod[%d] icycle[%d]\n", nchares, nlop, nmod, icycle);

    unsigned int ninterval=0;
    boost::random::mt19937 generator;
    generator.seed(123456+time(NULL)+nchares+icycle);
    boost::uniform_int<> interval(1,numDets-1);
    ninterval = interval(generator);

    CkPrintf(" The GA crossover interval is [%d] with mutation ratio is [%f] \n", ninterval,mutation_ratio);

    DetCharm=CProxy_DetCharmClass::ckNew(nchares);  
    DetCharm.Ent_genetic(determinants_charm,chareMPS,nchares,nlop,nmod,icycle,ninterval,EntMat);    

  }

//  void assemble_scratch(std::vector<int> order)
  void assemble_scratch(){
  
    std::ifstream fs;
    fs.open("DETs.scr");
    unsigned int nDets;
    int detorb[nsite_mps];
    double detcoeff;

//    std::vector<int> order;
    int j=0;    
    for(int i=0; i<nsite_mps; i++){
      fs >>  j;
//      std::cout << " " << j;
    }    
    fs >>  detcoeff;
//    std::cout <<  " : " << detcoeff << std::endl;

//    system("cat DETs.scr");
//    cout << std::endl;

//    string dets_show[nDets]; //dets represent
  
  }	  

  void mps2ci_charm(std::vector<std::vector<ChargeCharm>> determinants_charm,  MPSCharm chareMPS, std::vector<ResultCharm> doneCharm){

    DetCharm=CProxy_DetCharmClass::ckNew(numDets);
    CProxy_SyncClass     DetSync =CProxy_SyncClass::ckNew(numDets);
    double TL1[numDets];

    double coeff;
    int i=0;
    coeff=DetCharm[i].ckLocal()->coeff;

//    DetCharm.mps2ci(determinants_charm,chareMPS);

//    CkCallback cb(CkIndex_DetCharmClass::(NULL), thisProxy);
//    coeff=DetCharm[i]->coeff;
    CkPrintf(" cklocal[%d] ci : %f \n", i, coeff);

//    for (int i=0; i< numDets; i++){
//      coeff=DetCharm[i].ckLocal()->coeff;
//      CkPrintf(" cklocal[%d] ci : %f \n", i, coeff);
//    } 
  }	
  	 
  void report_coeff2(double result)
  {
    ++numRecv;	  
   CkPrintf("Summation of ci^2 : %f \n", result);
   if (numRecv == nloops_main+1) {
      CkPrintf("MPS to CI only");
      CkExit();
    }

//   if(result>CI_completeness){
//     CkExit();}
//   else{        
//     }
   }

  

  void finished() {
    ++numRecv;	  
    CkPrintf("numRecv in CkCallback : %d \n", numRecv);
    std::string command="date";
    int isys=system(command.c_str());
    if(numRecv==1)
      assemble_scratch();
//    if (numRecv == nloops_main+1) {
//      CkPrintf("MPS to CI only"); 
//      CkExit();
//    }
   }

    void report_coeff(double result)
    {
     CkPrintf("====>  report coeff  : %f \n", result);
     }

//    void monte_carlo_wrapper(std::vector<std::vector<ChargeCharm>> &inDets, MPSCharm &inMPS, std::vector<ResultCharm> out)
    void monte_carlo_wrapper(std::vector<std::vector<ChargeCharm>> determinants_charm, MPSCharm chareMPS, std::vector<ResultCharm> doneCharm){
      CkPrintf("====>  monte_carlo_wrapper  :  \n");
//      DetCharm.monte_carlo1(determinants_charm, chareMPS, doneCharm);
    }

    void get_entanglement(std::string infile, std::vector<std::vector<double>>& EntMat, std::vector<int> order){

	std::vector<int> Rorder;
	Rorder.resize(nsite_mps);
        for(int i = 0; i < nsite_mps; i++)
        {
          Rorder[order[i] - 1] = i;
        }

//        for(int i = 0; i < nsite_mps; i++)
//        {
//          std::cout << Rorder[i] << " " ;
//	}
//        std::cout << " " << std::endl;

	std::string tmpStr;
        std::ifstream entfile(infile);
        int idx1, idx2;
        double entValue;
        getline(entfile, tmpStr); // skip the 1st line
        while(entfile >> idx1 >> idx2 >> entValue)
          EntMat[Rorder[idx1 - 1]][Rorder[idx2 - 1]] = (entValue >= 0 ? entValue : 0);
        entfile.close(); // read-done

        for(int i = 0; i < nsite_mps; i++)
        {
          for(int j = 0; j < nsite_mps; j++){
	    if(EntMat[i][j]<0)	  
              EntMat[i][j] = 0.000001;
	  }
        }

/*	
        // normalization 
        for(int i = 0; i < nsite_mps; i++)
        {
//          for(int j = 0; j < nsite_mps; j++)
//            std::cout << EntMat[i][j] << " " ;
//          std::cout << " " << std::endl ;
          entValue=0.0;
          for(int j = 0; j < nsite_mps; j++)
            entValue = entValue + EntMat[i][j]*EntMat[i][j];
          for(int j = 0; j < nsite_mps; j++)
            EntMat[i][j] = EntMat[i][j]/sqrt(entValue);
        }

        // averaged to 1
        for(int i = 0; i < nsite_mps; i++)
        {
          entValue=0.0;
          for(int j = 0; j < nsite_mps; j++)
            entValue = entValue + EntMat[i][j];
          for(int j = 0; j < nsite_mps; j++)
            EntMat[i][j] = EntMat[i][j]/entValue;
        }

        // distribution for roulette
        for(int i = 0; i < nsite_mps; i++)
        {
          entValue=0.0;
          for(int j = 0; j < nsite_mps; j++){
            entValue = EntMat[i][j] + entValue;
	    EntMat[i][j] = entValue;
	  }
        }
*/	

    }

  std::vector<std::vector<ChargeCharm>> determinants_charm; 
  MPSCharm chareMPS;
  std::vector<ResultCharm> doneCharm;

};

#define CK_TEMPLATES_ONLY
#include "EDGA_charmc.def.h"
#undef CK_TEMPLATES_ONLY
