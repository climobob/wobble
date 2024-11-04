#include "grid_math.h"

#define NX 94
#define NY 192

int main(int argc, char *argv[]) {
  grid2<float> x(NX, NY);
  grid2<double> sx(NX, NY), sxx(NX, NY);
  grid2<float> sxy(x.xpoints()*x.ypoints(), x.xpoints()*x.ypoints() );
//  grid2<int> index(x.xpoints()*x.ypoints(), x.xpoints()*x.ypoints() );

  FILE *fin;
  int i, j, count = 0;
  ijpt loc;

  fin = fopen("temps","r");
  sx.set((double) 0.0);
  sxy.set((double) 0.0);
  count = 0;
  for (loc.i = 0;     loc.i < sxy.xpoints(); loc.i++) {
  for (loc.j = loc.i; loc.j < sxy.ypoints(); loc.j++) {
 //   index[loc] = count;
    count++;
  }
  }

  count = 0;
  while (!feof(fin) && count < 8) {
    x.binin(fin);
    if (!feof(fin)) {

      count += 1;

      for (loc.j = 0; loc.j < sxx.ypoints() ; loc.j++) {
      for (loc.i = 0; loc.i < sxx.xpoints() ; loc.i++) {
        sx[loc] += (double) x[loc];
        sxx[loc] += x[loc]*x[loc];
      }
      }

      for (loc.i = 0; loc.i < sxy.xpoints(); loc.i++) {
      // only do from i to ny, as the grids are symmetric (and this is slow)
      for (loc.j = loc.i; loc.j < sxy.ypoints(); loc.j++) {
        sxy[loc] += x[loc.i]*x[loc.j];
      }
      }

    }
  }
  fclose(fin);
  printf("Found %d grids\n",count); fflush(stdout);

  //return 0;


  sx /= (double) count;
  printf("average of all files %f %f %f\n",sx.average(), sx.gridmax(), sx.gridmin() );
  fflush(stdout);

// also output mean, variance, rms, terms used in divisions
  fin = fopen("stats","w");
  sx.binout(fin);
  sxx.binout(fin);
  sxy.binout(fin);
  fclose(fin);

  fin = fopen("stats2", "w");
  double a, b, r;
  ijpt xyloc;
  printf("sxy.xpoints, ypoints %d %d\n",sxy.xpoints(), sxy.ypoints()); fflush(stdout);
  for (i = 0; i < sxy.xpoints(); i++) {
  for (j = i; j < sxy.ypoints(); j++) {
    xyloc.i = i; xyloc.j = j;

    b = (sxy[xyloc] - count * sx[i]*sx[j])/(sxx[i] - count*sx[i]*sx[i]);
    a = sx[j] - b*sx[i];
    r = (sxy[xyloc] - count * sx[i]*sx[j])/sqrt(sxx[i] - count*sx[i]*sx[i])/
                                           sqrt(sxx[j] - count*sx[j]*sx[j]);

    fwrite(&i, sizeof(int), 1, fin);
    fwrite(&j, sizeof(int), 1, fin);
    fwrite(&a, sizeof(double), 1, fin);
    fwrite(&b, sizeof(double), 1, fin);
    fwrite(&r, sizeof(double), 1, fin);

  }
  printf("i %5d j %5d  a b r %f %f %f\n",i, j, a, b, r);
  }


  return 0;
}  
