#include <vector>
#include <pup_stl.h>
#include "dmrg/block_matrix/symmetry/nu1pg.h"
#include "CompatibilityWrapper.h"

class ChargeCharm : public NU1ChargePG<2> 
{

 public:
    void pup(PUP::er &p){
       PUParray(p,data_,2); 
       }
    inline ChargeCharm &operator=(const ChargeCharm &inDet){
       data_[0]=inDet.data_[0];
       data_[1]=inDet.data_[1];
       data_[2]=inDet.data_[2];
       return *this;
      }
};

class ResultCharm
{ 
 public:
    std::vector<ChargeCharm> det; 
    double val;

    void pup(PUP::er &p){
      p|val;
      p|det;
    }	    

};	

class MPSblockCharm
{
 public: 
    int isite;   
    int row;
    int col; 
    std::vector<std::vector<double>> mat;

    void pup(PUP::er &p){
      p|isite;
      p|row;
      p|col;
      p|mat;
    }    

  void get_rc(const int &in1, const int &in2){
    row=in1;
    col=in2;
  }
/*
   void mat_resize()  {
     mat.resize(row);
     for(int i=0; i<col; i++)
        mat[i].resize(col);
  }
*/
};

class MPSbasisCharm
{
 public:
   ChargeCharm lc;
   ChargeCharm rc;
   int ls;
   int rs;

   void pup(PUP::er &p){
     p|ls;
     p|rs;
     p|lc;
     p|rc;
   }

};

class MPSdimCharm
{
 public:
   ChargeCharm charge;
   int ndim;

   void pup(PUP::er &p){
     p|ndim;
     p|charge;
   }

};

class MPSsiteCharm
{
 public: 
    std::vector<MPSblockCharm> data;
    std::vector<MPSbasisCharm> basis;
    std::vector<MPSdimCharm> site_dim;
    std::vector<MPSdimCharm>  row_dim;
    std::vector<MPSdimCharm>  col_dim;

    void pup(PUP::er &p){
      p|data;
      p|basis;
      p|site_dim;
      p|row_dim;
      p|col_dim;
    }

};

/*
template<class MPSsiteCharm>
inline void PUParray(PUP::er &p,MPSsiteCharm *array,int length) {
   for (int i=0;i<length;i++) p|array[i];
}
*/

//template<class MPSsiteCharm>
//inline void operator|(PUP::er &p, MPSsiteCharm &c){
//   p|c;
//}
//
/*
class test2
{
  public:
  double YY;
  void pup(PUP::er &p){
    p|YY;
  }
   
};

class MPStest
{
  public :
  test2 XX;
  void pup(PUP::er &p){
    p|XX;
  }
};
*/

class MPSCharm : public MPS<matrix, TwoU1PG> 
{
 public:

    int   charm_nmps;
    int  *charm_mps_size;
    std::vector<MPSsiteCharm> charm_mps;
//    std::vector<MPStest> test_mps;
//    int  *charm_mps_size = new int[charm_nmps];

    MPSCharm(void){} 

    void length_charm(const int &inlen){
      charm_nmps = inlen;
      charm_mps_size = new int[charm_nmps]; 
    }

    void size_charm(const int &isize, const int &insize){
      charm_mps_size[isize]=insize;
    }
   
    void charm_mps_pre();

    void pup(PUP::er &p){
      p|charm_nmps;
      if(p.isUnpacking()) charm_mps_size = new int[charm_nmps];
      p(charm_mps_size, charm_nmps);
      p|charm_mps;
//      PUParray(p, charm_mps, charm_nmps) ; 
//      p|*charm_mps_size;  
    }

    ~MPSCharm(){}

//    void mps_rc_charm(){
//       for(int i=0;i<charm_nmps;i++)
//         charm_mps_rdim[i]= new int [charm_mps_size[i]];
//    }
// private:

};


void MPSCharm::charm_mps_pre()
{
//  block_matrix<Matrix, SymmGroup> t = (*this)[length()-1].normalize_left(DefaultSolver());
     (*this).charm_mps.resize(charm_nmps);
};


class IndexCharm : public Index<TwoU1PG>
{
};


class BlockMatCharm : public block_matrix<matrix, TwoU1PG>
{
 public:
   BlockMatCharm& operator=(BlockMatCharm rhs);
   BlockMatCharm& operator=(block_matrix<matrix, TwoU1PG> rhs);
};


