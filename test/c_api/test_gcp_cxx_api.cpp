/* This file is part of mctc-gcp.
 * SPDX-Identifier: LGPL-3.0-or-later
 *
 * mctc-gcp is free software: you can redistribute it and/or modify it under
 * the terms of the Lesser GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * mctc-gcp is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * Lesser GNU General Public License for more details.
 *
 * You should have received a copy of the Lesser GNU General Public License
 * along with mctc-gcp.  If not, see <https://www.gnu.org/licenses/>.
 */

#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <cstring>
#include <string>
#include <vector>
#include <tuple>
#include <gcp.h>

typedef struct {
  std::vector<int> iz;
  std::vector<double> coords;
  std::vector<double> lattice;
  bool pbc;
} Molecule;

bool test_generic(const std::string& title, Molecule& mol, const std::string& method, const double energy_ref){

  char method_f[20];
  strcpy(method_f,method.c_str());
  for(size_t ii=method.size();ii<20;ii++) method_f[ii] = ' ';

  bool dograd  = false;
  bool dohess  = false;
  bool echo    = false;
  bool parfile = false;
  bool pbc     = mol.pbc;

  auto& lat = mol.lattice;
  auto& xyz = mol.coords;
  auto& iz  = mol.iz;
  int   nat = (int)(xyz.size()/3);
  

  double energy = 0.e0;
  std::vector<double> gradient(3*nat,0.e0);
  std::vector<double> gradlat(9,0.e0);

  c_gcp_call(&nat,&xyz[0],&lat[0],&iz[0],&energy,&gradient[0],&gradlat[0],&dograd,&dohess,&pbc,method_f,&echo,&parfile);

  double de = energy - energy_ref;
  
  std::string stat = "\x1B[32m [OK]\x1B[0m";

  bool is_ok = true;
  if (fabs(de) > 1e-9){
    stat = "\x1B[31m [FAILED]\x1B[0m";
    is_ok = false;
  }
  printf("E: %30.20f %20s     %s\n",energy,stat.c_str(),title.c_str());

  return is_ok;

}

int main(int argc, char* argv[]){

  std::vector<std::tuple<std::string,Molecule,std::string,double> > test_systems = {
    {
      "Alanine@HF3c",
      {
        {1,7,6,6,8,6,1,7,6,6,8,6,8,1,1,1,1,1,1,1,1,1,1},
        {-1.93145285444792,+3.69743396003490,+5.07672663197993,
         -2.38225794234874,+1.87122793226519,+5.47939866438837,
         -2.38363826594372,+0.46945474235908,+3.10066131453095,
         +0.05491295406579,+0.84981606018213,+1.57533328139901,
         +1.15553526565357,+2.88302479323393,+1.49584893725220,
         -4.60525415514495,+1.33546037841295,+1.46093111812950,
         -0.94935967135323,+1.23895283176357,+6.59006017993275,
         +0.90719609031382,-1.21518242048549,+0.28080129333493,
         +2.94491704712175,-0.96790880864505,-1.53990197241268,
         +2.33465737391659,+0.79844258337069,-3.73912178828922,
         +3.92673604793746,+1.81703260886247,-5.01158300414526,
         +3.66032988644352,-3.56018817506558,-2.59948629629987,
         -0.17835566103986,+0.97841613903525,-4.20105316431806,
         -0.35820586681744,+2.10507058522216,-5.63692062300423,
         -0.23557408796486,-2.73553018128516,+0.13945189803826,
         -2.61502607287587,-1.54956844388107,+3.54379863895688,
         -4.40788990832190,+3.34964157223239,+1.02628696163526,
         -6.37737290628970,+1.04790338926862,+2.48347133798594,
         -4.67147573464799,+0.30228459708274,-0.33039253201856,
         +4.58101416587973,-0.11358970578901,-0.60029127792819,
         +2.05572744486357,-4.42625244199752,-3.58274737172425,
         +4.26057887986804,-4.81176728720273,-1.06753470682591,
         +5.21404372180499,-3.36409206395890,-3.94412805049104},
        {0.,0.,0.,0.,0.,0.,0.,0.,0.},
        false
      },
     "hf3c",
      -0.12313247650779590714e0
    },
    {
      "Alanine@R2SCAN-3c",
      {
        {1,7,6,6,8,6,1,7,6,6,8,6,8,1,1,1,1,1,1,1,1,1,1},
        {-1.93145285444792,+3.69743396003490,+5.07672663197993,
         -2.38225794234874,+1.87122793226519,+5.47939866438837,
         -2.38363826594372,+0.46945474235908,+3.10066131453095,
         +0.05491295406579,+0.84981606018213,+1.57533328139901,
         +1.15553526565357,+2.88302479323393,+1.49584893725220,
         -4.60525415514495,+1.33546037841295,+1.46093111812950,
         -0.94935967135323,+1.23895283176357,+6.59006017993275,
         +0.90719609031382,-1.21518242048549,+0.28080129333493,
         +2.94491704712175,-0.96790880864505,-1.53990197241268,
         +2.33465737391659,+0.79844258337069,-3.73912178828922,
         +3.92673604793746,+1.81703260886247,-5.01158300414526,
         +3.66032988644352,-3.56018817506558,-2.59948629629987,
         -0.17835566103986,+0.97841613903525,-4.20105316431806,
         -0.35820586681744,+2.10507058522216,-5.63692062300423,
         -0.23557408796486,-2.73553018128516,+0.13945189803826,
         -2.61502607287587,-1.54956844388107,+3.54379863895688,
         -4.40788990832190,+3.34964157223239,+1.02628696163526,
         -6.37737290628970,+1.04790338926862,+2.48347133798594,
         -4.67147573464799,+0.30228459708274,-0.33039253201856,
         +4.58101416587973,-0.11358970578901,-0.60029127792819,
         +2.05572744486357,-4.42625244199752,-3.58274737172425,
         +4.26057887986804,-4.81176728720273,-1.06753470682591,
         +5.21404372180499,-3.36409206395890,-3.94412805049104},
        {0.,0.,0.,0.,0.,0.,0.,0.,0.},
        false
      },
     "r2scan3c",
      0.01836274475757706734e0
    },
  };

  printf("Run C-API tests...\n");
  size_t nok = 0;
  for(auto& mol : test_systems){
    if (test_generic(std::get<0>(mol),std::get<1>(mol),std::get<2>(mol),std::get<3>(mol))) nok++;
  }

  if (nok == test_systems.size())
    printf("All tests passed!\n");
  else
    printf("%li out of %li tests failed!\n",test_systems.size()-nok,test_systems.size());

  return 0;

}

