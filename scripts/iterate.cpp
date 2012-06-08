#include <math.h>
#include <string.h>
#include <iostream>
#include <fstream>
#include <cstdlib>
using namespace std;

void read_data(int n_intervals, double *divs, double **probs){
  char *pch; // a beolvasott bemeneti stringek darabkai
  int i, ii;
  char s[10000]; // max bemenet hossza
  //------------ a div tomb es probs tombok beolvasasa ------//
  // elso beolvasott sor: a div ertekek
  cin.getline(s,10000);
  pch = strtok(s," ,");
  ii=0;
  while(pch != NULL && ii<n_intervals ){
    divs[ii]=strtod(pch,NULL);
    pch = strtok (NULL, " ,");
    ii++;
  }
  // tovabbi sorok a probs 2d tomb ertekei, soronkent
  for(i = 0; i < n_intervals; i++){
    cin.getline(s,10000);
    pch = strtok(s," ,");
    ii=0;
    while(pch != NULL && ii<n_intervals ){
      probs[i][ii]=atof(pch);
      pch = strtok (NULL, " ,");
      ii++;
    }
  }
  return;
}

void read_data_debug(int n_intervals, double *divs, double **probs){
  //
  //  hamis adatbeolvaso fuggveny, debugolasra
  //
  int i,j;
  divs[0]=0;
  divs[1]=0.33333333333333331;
  divs[2]=0.66666666666666663;
  divs[3]=1.0;

  for(i=0;i<3;i++)
      for(j=0;j<3;j++)
          probs[i][j]=0.11111111;
}


void degdist(int n, int mindeg, int maxdeg, int n_intervals, double *divs, double **probs, double *rho)
{
    int i,j,k;
    int ii=0;

    double *lengths;
    lengths = (double *)calloc(n_intervals,sizeof(double));

    double *log_lengths;
    log_lengths = (double *)calloc(n_intervals,sizeof(double));

    double *avgdeg;
    avgdeg = (double *)calloc(n_intervals,sizeof(double));

    double *log_avgdeg;
    log_avgdeg = (double *)calloc(n_intervals,sizeof(double));

    double *log_factorial;
    log_factorial = (double *)calloc((maxdeg+1), sizeof(double));

    double *log_rho_i_length;
    log_rho_i_length = (double *)calloc(maxdeg+1-mindeg,sizeof(double));

    double *log_rho_i;
    log_rho_i = (double *)calloc((maxdeg+1-mindeg), sizeof(double));

    // atmeneti valtozok
    double avgdeg_tmp;

    //lengths = [divs[i] - divs[i-1] for i in range(1, n_intervals)]
    for(i=0; i<n_intervals;i++){
      lengths[i]=divs[i+1] - divs[i];
      //printf("degdist lengths[%d]=%f\n",i,lengths[i]);
    }
    // ez kulonbozik a pythontol, mert a divs-ben itt benne van a kezdo 0,
    // a python divs tomb elejen pedig nincs ott a kezdo nulla
    // lengths[0]=divs[1];

    for(i=0; i<n_intervals; i++){
        avgdeg_tmp = 0;
        for(j=0; j<n_intervals; j++){
          avgdeg_tmp += probs[i][j] * lengths[j];
        }
        avgdeg[i] = n*avgdeg_tmp;
        //printf("%d.=%16.10E ", i, avgdeg[i]);
        log_avgdeg[i] = log(avgdeg_tmp);
        log_lengths[i] = log(lengths[i]);
    }

    for(k=1; k<(maxdeg+1);k++)
        log_factorial[k]=0.5*log(2*M_PI)+(k+0.5)*log(k)-k;
    log_factorial[0]=0;

    for(i=0;i<n_intervals;i++){
        //  # Eq. 5
        //log_rho_i = [(k * math.log(avgdeg[i]) - log_factorial[k]  -  avgdeg[i])
        //     for k in range(mindeg, maxdeg+1)]
        for(ii=0,k=mindeg; k<maxdeg+1; k++,ii++){
          log_rho_i[ii]=k*log(avgdeg[i]) - log_factorial[k] - avgdeg[i];
          //printf("running degdist, log_rho_i[%d]=%f\n",ii, avgdeg[i]);
        }
        // log_rho_i_length = numpy.array(log_rho_i) + log_lengths[i]
        for(ii=0;ii<maxdeg+1;ii++){
          log_rho_i_length[ii]=log_rho_i[ii] + log_lengths[i];
        }
        // rho[mindeg:] += numpy.exp(log_rho_i_length)
        for(ii=0;ii<maxdeg+1;ii++){
          rho[mindeg+ii] += exp(log_rho_i_length[ii]);
        }
    }

    // ------------- tombok felszabaditasa ----------//
    free( lengths );
    free( avgdeg );
    free( log_avgdeg );
    free( log_factorial );
    free( log_rho_i_length );
    free( log_rho_i );
    return;
}


void iterate_divs(int K, int divs_size, double *divs, double *olddivs, double *newdivs){

   if(K==1){
     memcpy(newdivs,olddivs,divs_size*sizeof(double));
     return;
   }
   //----- rekurziv fv hivas uj valtozokkal -------//
   int i,k, olddivs_size,copypos;
   double start, stop, interval_length;
   double *newdivs_in; // ide kerul az eredmeny
   double *olddivs_in; // belso hasznalatra
   int iterated_divs_size=(int)pow(divs_size-1,K)+1; // aktualis K-hoz a divs meret

   newdivs_in=(double *)calloc(iterated_divs_size,sizeof(double)); // ide kerul az eredmeny
   olddivs_in=(double *)calloc(iterated_divs_size,sizeof(double)); // belso hasznalatra

   // iterate_divs(fokszam-1, eredeti divs tomb pointer, eredmeny(iteralt) divs tomb pointer)
   iterate_divs(K-1, divs_size, divs, olddivs, newdivs_in);
   //----------------------------------------------------------------------------------------------//
   //------- a newdivs_in -bol atmasoljuk a fuggvenyhivas eredmenyet a lokalis olddivs_in -be -----//
   // memcpy(destination, source, size )
   memcpy(olddivs_in,newdivs_in,(pow(divs_size-1, K)+1)*sizeof(double) );

   olddivs_size=pow((divs_size-1),K-1);
   copypos=1;
   newdivs[0]=0;
   for(i=0;i<olddivs_size;i++){

       start=olddivs_in[i];
       stop=olddivs_in[i+1];
       interval_length=stop-start;

       // a k-t 1-tol inditjuk, hogy a divs kezdo 0-jat ne szamolja bele
       for(k=1; k<divs_size;k++){
           newdivs[copypos]=start+interval_length*divs[k];
           copypos++;
       }
   }

   free(newdivs_in);
   free(olddivs_in);

   return;
 }

void iterate_probs(int m, int K, double **self_probs, double **iterated_probs){

   int i, j, ii, it;
   int iteration;
   int indices_size;
   indices_size = pow(m,K);

   int **indices_dict;
   indices_dict = (int **)calloc(indices_size,sizeof(int*));
   for(i = 0; i < indices_size; i++)
     indices_dict[i] = (int*)calloc(K,sizeof(int));

   int *indices2;
   indices2 = (int *)calloc(indices_size,sizeof(int));
   //probs = numpy.ones((m**K, m**K))
   for(i=0;i<indices_size;i++)
       for(j=0;j<indices_size;j++)
           iterated_probs[i][j]=1;

     // for iteration in xrange(K):
     //   indices2 = (indices % (m**(iteration+1))) // (m**iteration)
     for(iteration=0; iteration<K; iteration++){
         for(ii=0; ii<indices_size; ii++){
           indices2[ii]= floor( ( ii % (int)pow(m,iteration+1) ) / pow(m,iteration) );
         }
         for(i=0; i<indices_size; i++){
           indices_dict[i][iteration]=indices2[i];
         }
     }
     // probs szamolas
     //
     //   for i in indices:
     //       ilist = indices_dict[i]
     //       for j in indices:
     //           jlist = indices_dict[j]
     //           for it in xrange(K):
     //               probs[i,j] *= self.probs[ilist[it], jlist[it]]
     for(i=0; i<indices_size; i++){
         for(j=0; j<indices_size; j++){
             for(it=0; it<K; it++){
                 iterated_probs[i][j] *= self_probs[ indices_dict[i][it] ][ indices_dict[j][it] ];
             }
         }
     }
     /*
     for(i=0; i<indices_size; i++){
       for(j=0; j<indices_size; j++){
         printf("%2.10e ",iterated_probs[i][j]);
       }
       printf("\n");
     }
     */
     for(i = 0; i < indices_size; i++)
         free( indices_dict[i]);
     free(indices_dict);
     free( indices2 );

     return;
 }

int main(int argc, char **argv)
{
    // program parameterek
    int n_intervals;
    int maxdeg;
    int n;
    int mindeg;
    //
    int i;
    int K=4;
    int m;//m=len(divs)
    int iterated_probs_size, iterated_divs_size;
    char type='a';

    if (argc < 5){
        printf("usage: %s <type> <n_intervals> <K> <maxdeg> <n> <mindeg>\n", argv[0]);
        exit(-1);
    }
    type = argv[1][0];
    if (type == 'd') {
        m = atoi(argv[2]);
        n_intervals = m + 1;
        K = atoi(argv[3]);
        maxdeg = atoi(argv[4]);
        n = atoi(argv[5]);
        mindeg = atoi(argv[6]);
    }

    //vvvvvvvvvvvvvvvvvvvvvvvvvv debug bemeneti parameterek -----//
    /*
    m=3;
    n_intervals=m+1;
    maxdeg=200;
    n=2000;
    mindeg=0;
    */
    //^^^^^^^^^^^^^^^^^^^^^^^^^^ debug bemeneti parameterek -----//

    //-----------------------------------------------------------//
    //     ezeket a tomboket kapja meg a pythontol               //
    //-----------------------------------------------------------//
    double *old_divs;
    double **old_probs;
    old_divs = (double *)calloc(n_intervals,sizeof(double));
    old_probs = (double **)calloc(n_intervals,sizeof(double*));
    for(i = 0; i < n_intervals; i++)
        old_probs[i] = (double *)calloc(n_intervals,sizeof(double));
    //-----------------------------------------------------------//
    //     az iteralt divs es provs tombok                       //
    //-----------------------------------------------------------//
    iterated_probs_size = (int)pow(m,K);
    // az iteralt divs merete a kezdo 0-val egyutt, pl
    // divs={0,0333,0.666,1} -nel n_intervals=4, K=2 -re ez (4-1)^2+1=10
    iterated_divs_size = (int)pow(m,K)+1;
    // az uj divs es probs tombok
    double *iterated_divs;
    double **iterated_probs;
    iterated_divs = (double *)calloc(iterated_divs_size,sizeof(double));
    iterated_probs = (double **)calloc(iterated_probs_size,sizeof(double*));
    for(i = 0; i < iterated_probs_size; i++){
        iterated_probs[i] = (double *)calloc(iterated_probs_size,sizeof(double));
    }

    //-----------------------------------------------------------//
    //     az eloszlas rho tomb                                  //
    //-----------------------------------------------------------//
    double *rho;
    rho = (double *)calloc((maxdeg+1),sizeof(double));


    //-----------------------------------------------------------//
    //-----------------------------------------------------------//

    // 1. beolvassuk a kezdeti adatokat az old_divs es old_probs tombokbe
    read_data(n_intervals, old_divs, old_probs);

    // 2. iteraljuk a divs tombot
    // iterate_divs(K,
    //              eredeti divs tomb merete,
    //              eredeti divs tomb pointere,
    //              az iteracio elozo lepeseben kapott divs tomb pointere (kezdetben az eredeti divs)
    //              eredmeny(iteralt) divs tomb pointere)
    iterate_divs(K,  n_intervals,  old_divs,  old_divs,  iterated_divs);

    // 3. iteraljuk a probs tombot
    iterate_probs(m,K,old_probs, iterated_probs);
    //for (i=0; i<iterated_divs_size; i++)
    i=70;
            //printf("iterated_divs[%d] %15.10e ", i, iterated_divs[i]);

    // 4. kiszamoljuk az eloszlast
    //void degdist(int n, int mindeg, int maxdeg, int n_intervals, double *divs, double **probs, double *rho)
    degdist(n, mindeg, maxdeg, iterated_divs_size-1, iterated_divs, iterated_probs, rho);

    // 5. kiirjuk a rho -t a std outputra
    for(i=0; i<maxdeg+1;i++)
      printf("%E ", rho[i]);

    //-----------------------------------------------------------//
    free(old_divs);
    for(i = 0; i < n_intervals; i++){
        free( old_probs[i]);
     }
    free(old_probs);
    free(iterated_divs);
    for (i = 0; i < iterated_probs_size; i++){
        free( iterated_probs[i]);
    }
    free(iterated_probs);

    return 0;
}

