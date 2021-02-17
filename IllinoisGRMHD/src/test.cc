#include <iostream>

int main() {

  int cctk_lsh[3] = {10,10,10};
  int i = 5;
  int j = 5;
  int k = 5;

  // Auxiliary variables
  const int Nx0 = cctk_lsh[0]-1;
  const int Nx1 = cctk_lsh[1]-1;
  const int Nx2 = cctk_lsh[2]-1;

  // Set the index limits for the averaging loop
  // Minimum values
  const int imin = (i!=0)*(i-1); // i=0 => iimin=0. Otherwise iimin = i-1
  const int jmin = (j!=0)*(j-1); // j=0 => jjmin=0. Otherwise jjmin = j-1
  const int kmin = (k!=0)*(k-1); // k=0 => kkmin=0. Otherwise kkmin = k-1

  // Maximum values  
  const int imax = i + 1*(i!=Nx0); // i=Nx0 => iimax = Nx0. Otherwise iimax = i+1
  const int jmax = j + 1*(j!=Nx1); // j=Nx1 => jjmax = Nx1. Otherwise jjmax = i+1
  const int kmax = k + 1*(k!=Nx2); // k=Nx2 => kkmax = Nx2. Otherwise kkmax = i+1

  int nn = (imax-imin+1)*(jmax-jmin+1)*(kmax-kmin+1);

  int num_neighbors=0;
  for(int kk=kmin;kk<=kmax;kk++)
    for(int jj=jmin;jj<=jmax;jj++)
      for(int ii=imin;ii<=imax;ii++)
        num_neighbors++;

  printf("%d %d %d | %d %d %d | %d %d %d | num_neighbors = %d (expected %d)\n",imin,i,imax, jmin,j,jmax, kmin,k,kmax,num_neighbors,nn);

  return 0;
  
}
