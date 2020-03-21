#include <cmath>
#include <iterator>
#include <iostream>
#include <fstream>
#include <string>
#include <time.h>
#include <sys/stat.h>
#include <chrono>
#include <vector>
#include <map>
#include <boost/lexical_cast.hpp>
#include <boost/tokenizer.hpp>
#include <boost/tuple/tuple.hpp>
#include <boost/random.hpp>

using std::cerr;
using std::cout;
using std::endl;

int exec_cmd(std::string cmd, std::string &res){
  if (cmd.size() == 0){  //cmd is empty 
    return -1;
  }
 
  char buffer[1024] = {0};
  std::string result = "";
  FILE *pin = popen(cmd.c_str(), "r");
  if (!pin) { //popen failed 
    return -1;
  }
 
  res.clear();
  while(!feof(pin)){
    if(fgets(buffer, sizeof(buffer), pin) != NULL){
      result += buffer;
    }
  }
 
  res = result;
  return pclose(pin); //-1:pclose failed; else shell ret
}

void quicksort(std::string dets[], double b[], int left, int right) {

    double pivot = fabs(b[(left+right)/2]);
    int l = left;
    int r = right;
    while(l <= r) {
     while (fabs(b[l]) < pivot) l++;
     while (fabs(b[r]) > pivot) r--;
     if (l <= r) {

        double tmp;
        tmp = b[l];
        b[l]= b[r];
        b[r]= tmp;

	std::string ctmp;
        ctmp   = dets[l];
        dets[l]= dets[r];
        dets[r]= ctmp;

        l++;
        r--;
     }
    };

    if (left < r)  quicksort(dets, b, left, r);
    if (l < right) quicksort(dets, b, l, right);
}


int main(int argc, char ** argv)
{
    std::string ONV_folder;
    std::string ORD_file;
    int Nsites; 

    if (argc < 4)
    {
        std::cout << " Usage :  ./result_sort  [ONVs_folder]  [N-particles] [Ordering]" << std::endl;
        std::cout << "  " << std::endl;
        exit(1);
    }

    ONV_folder = std::string(argv[1]);
    Nsites     = atof(argv[2]);    
    ORD_file   = std::string(argv[3]);

    // Orbital ordering
    int order[Nsites];
    std::ifstream fs0;
    fs0.open(ORD_file);
    for(int i=0; i<Nsites; i++){
      fs0 >> order[i];
//      std::cout <<  order[i] << " " << std::endl;
    }

    // This is used for ONVs
    typedef std::map<std::string, double> Hash_Map_with_value;
    typename Hash_Map_with_value::iterator iter;
    Hash_Map_with_value hash;

    // ONVs index 
    typedef std::map<std::string, long>   Hash_Map_with_index;
    typename Hash_Map_with_index::iterator iter_idx;
    Hash_Map_with_index hash_index;

    // Read-in to Hash/Map
    std::string res;   
    std::string command;
    command="ls "+ ONV_folder +  " -l | grep \"DETs.scr\" | wc -l ";
    exec_cmd(command,res);
    int nfile=stoi(res,0,10);
    command="ls " + ONV_folder + " -l |awk '/DETs/ {print $NF}' > ONVfile.scr";
    int isys=system(command.c_str());

    std::ifstream fsONV;
    fsONV.open("ONVfile.scr");

    unsigned int icounter=0;
    for(int i=0; i<nfile; i++)
    {
      std::string ONV;
      double Donv;
      std::string temp;

      std::string ONV_file;   
      fsONV >> ONV_file;
//      std::cout << "ONV_file : " << ONV_file << std::endl;
      ONV_file=ONV_folder+"/"+ONV_file;

      std::ifstream fs;
      fs.open(ONV_file);

      fs >> ONV >> Donv ;
      hash[ONV]=Donv;
      hash_index[ONV]=icounter;
      while(getline(fs,temp))
      {
        fs >> ONV >> Donv ;
        iter = hash.find(ONV);      
        if(iter == hash.end()){
          icounter=icounter+1;
          hash[ONV]=Donv;
          hash_index[ONV]=icounter;
        }
      }
      icounter=icounter+1;
    }
    icounter=icounter-1;

    double sum_ci2=0.0;
    for(iter=hash.begin(); iter!=hash.end(); iter++){
      sum_ci2=sum_ci2+pow(iter->second,2.0);    
//      std::cout << " " <<  iter->first << " " << iter->second << std::endl; 
    } 

    std::cout << "The total CI^2 : " << sum_ci2 << ", with " << icounter+1 << " ONVs" << std::endl; 

    double    CIs_show[hash.size()];
    std::string   dets_show[hash.size()];
    unsigned int ipos=0;
    for(iter=hash.begin(); iter!=hash.end(); iter++){
      char chtmp[Nsites+1];
      std::string ctmp;
      CIs_show[ipos] = iter->second;
      dets_show[ipos]= iter->first;   
      for(int i=0; i<Nsites; i++){
	int ii=order[i]-1;
        chtmp[ii]=iter->first[i]; 	
      }
      for(int i = 0; i < Nsites; i++)
        ctmp=ctmp+chtmp[i];
      dets_show[ipos]= ctmp;   
      ipos=ipos+1;
    }

    quicksort(dets_show, CIs_show,0, hash.size()-1);

    for(unsigned int i=0; i<hash.size() ; i++)
      std::cout << " The sorted ONVs: " << dets_show[hash.size()-i-1] << " " << CIs_show[hash.size()-i-1] << std::endl;

    std::ofstream ofs;
    ofs.open("ONV_int.restart");
    for(unsigned int i=0; i<hash.size() ; i++)
      ofs <<  dets_show[hash.size()-i-1] << std::endl;

}


