#include<math.h>
#include<stdio.h>
#include<stdlib.h>

//****** Parametros do gerador aleatorio ********************//
#define IA 16807
#define IM 2147483647
#define AM (1.0/IM)
#define IQ 127773
#define IR 2836
#define NTAB 32
#define NDIV (1+(IM-1)/NTAB)
#define EPS 1.2e-7
#define RNMX (1.0-EPS)

#define tol1 0.00001   /////// tolerancia para identificar os regimes 
#define tol2 0.01

//******   Parametros do sistema ********************//
#define NN 1          // numero of synapses
#define N 4             // numero de equacoes
#define h 0.01         // passo de integracao
#define nn 300000       //300000
#define n nn/h // total number of steps n*h
//#define g_exc 0.15      // acoplamento

//*******gerador de numeros aleatorios******//
float ran1(long *idum);
#define NR_END 1
#define FREE_ARG char*

//*******ponteiros para alocar memoria*****///
void nrerror(char error_text[]);
int *vector(long nl,long nh);
int **imatrix(long nrl, long nrh, long ncl, long nch);
void free_vector(int *v, long nl, long nh);
void free_imatrix(int **m, long nrl, long nrh, long ncl, long nch);
double *dvector(long nl,long nh);
void free_dvector(double *v, long nl, long nh);
void derivs(double y[],double df[]); 
double tau_rec,tau_fac, tau_ina, U;

FILE *o,*p,*q,*s;

int main(void)
{
  double *df,*y,*a,*b,*c,*x;

  int i,j,k,t,max_spike,cont;
  double tempo,tspike,freq,yampl,ymax,yant,yant2,iay,ux;
  int facil,depre,regime,cont_spike,first_facil,first_depre,y0;

  long idum; //semente para gerar números aleatórios
  
  idum=-123456789;

  y=dvector(1,N*NN+1);
  df=dvector(1,N*NN+1);
  x=dvector(1,N*NN+1);
  a=dvector(1,N*NN+1);
  b=dvector(1,N*NN+1);
  c=dvector(1,N*NN+1);

  o=fopen("Dinamics1.dat","wt"); 
  p=fopen("Dinamics2.dat","wt"); 
  q=fopen("Dinamics3.dat","wt");
  s=fopen("Parameter_space.dat","wt"); 
 
  if(o==NULL)
    {
     puts ("erro no arquivo");
     getchar();
     exit(1);
    }
    if(p==NULL)
    {
     puts ("erro no arquivo");
     getchar();
     exit(1);
    }
    if(q==NULL)
    {
     puts ("erro no arquivo");
     getchar();
     exit(1);
    }
    
   //////////////////////// Parameters ////////////////  
   tau_rec=800.0;  /////////// recuperação 800ms
   tau_fac=1000.0;  ///////////facilitação 1000ms
   tau_ina=3.0;  ////////// inativation 
   U=0.95;
   freq=2.0;
    
 //  for(U=0.0;U<=1.005;U=U+0.005)
 //  for(freq=0.01;freq<=10.01;freq=freq+1.0)
   {
   //////////////////////////////////////////////////
   ////                                          ////
   ////          Spike protocoll                 ////
   ////                                          ////
   //////////////////////////////////////////////////
 //  for(j=1;j<=max_spike;j++) 
//   {
 //  tspike = (1000/(freq))*j ;
 //  printf("%f \n", tspike[j]);   
 //  } 

//****************  Condicoes iniciais  **************************   
	 x[1]=1.0; //// x avaliable     x + y + z = 1  !!!! always !!!
	 x[2]=0.0; //// y active 
	 x[3]=0.0; //// z inactive
	 x[4]=U; //// fraction of avaliable vesicules which become actives   
     k=0;
     ymax=0;  yant=0.0;
     facil=0;
     depre=0; first_facil=0; first_depre=0;
     regime=0;cont_spike=0; y0=0;
     iay=0;
     tspike=0;
      //********************* LOOP DO TEMPO  ***************************
     tempo=0.0; 
   for(t=1;t<=n;t++)  
	{                   
	tempo=tempo+h;      //em milisegundos
	yant2=yant;         //// amplitude 2 vezes anterios a atual
	//////////////////////////// Atualiza disparo ////////////////////
	if(tempo>=tspike+1000/freq && tempo<tspike+1000/freq+h) 
	{
	ux =	x[1]*x[4];          //// adicionado/retirada a mesma quantidade
	x[1]= x[1] - ux;   //// x
	x[2]= x[2] + ux;   //// y
	//x[3]=x[3] ;             //// z
	x[4]=x[4] +U*(1-x[4]);	  //// u	
	
	//  printf("%f %f \n", tempo, tspike[j]); 
	 
	 /////////////////////////////////////////// regime analises
	 //// sala y amplitudes 
	tspike=tempo;
	 
	yampl=x[2];
	
	if(x[2]>0.001) y0=1;
	
	if(ymax<yampl) ymax = yampl; 
	
	printf("%f %f %f %f %f \n", tempo, yampl,yant,iay,ux); 
	
	if(yant>0)
	{     
		iay  = yampl-yant;
		if(iay>tol1) facil=1;
		if(iay<-tol1) depre=1; 
		
	//printf("%f %f %f %f \n", tempo, yampl,yant,iay); 
    
		if(facil==1 & depre==0) 
			{
			first_facil=1;
			}	
			
		if(facil==0 & depre==1) 
			{
			first_depre=1;
		   }	
	}
	 
	yant=x[2];

	 fprintf(o,"%f %f %f %f %f %f\n",tempo,x[1],x[2],x[3],x[4],iay);
	//printf("%f %f %d %d\n",ymax,yampl,facil, depre);
	cont_spike++;
	   }
	   
	   if(cont>100) {
//if(freq>0.501 && freq<1.5)	
//fprintf(o,"%f %f %f %f %f\n",tempo,y[1],y[2],y[3],y[4]);cont=0;
//if(freq<=2 && freq>1 & U>=0.05 && U<0.15)	fprintf(p,"%f %f %f %f %f\n",tempo,y[1],y[2],y[3],y[4]);cont=0;
//if(freq<=4 && freq>3)	fprintf(q,"%f %f %f %f %f\n",tempo,y[1],y[2],y[3],y[4]);cont=0;
	}
cont++;
	   
	   
	  
      for(i=1;i<=N;i++) y[i]=x[i]; 
      
	  // ------------ Integrador numérico Runge-Kutta 4ª ordem--------------- 
	  derivs(y,df);
	  for(i=1;i<=N*NN;i++)
	    {
	     a[i]=h*df[i];
	     y[i]=x[i]+a[i]/2.0;
	    }
	  derivs(y,df);
	  for(i=1;i<=N*NN;i++)
	    {
	     b[i]=h*df[i];
	     y[i]=x[i]+b[i]/2.0;
	    }
	  derivs(y,df);
	  for(i=1;i<=N*NN;i++)
	    {
	     c[i]=h*df[i];
	     y[i]=x[i]+c[i]; 
	    }       
	  derivs(y,df);
	  for(i=1;i<=N*NN;i++)
	    x[i]=x[i]+(a[i]+h*df[i])/6.0+(b[i]+c[i])/3.0; 	
	  //--------------------------------------------------------------------



if(x[1]+x[2]+x[3]>=1.001) { printf("Erro 1, Algo explodiu!!\n"); break;}

if(x[1]+x[2]+x[3]<=0.999) 
	{ 
	printf("Erro 1, não conservou!! Soma= %f\n",x[1]+x[2]+x[3]); break;
	}
	}  // fim loop do tempo	
	
    if(facil==1 && depre==0) regime = 1; ////// Facilitation Biphasico
	if(facil==0 && depre==1) regime = 2; ////// Depression
	if(facil==1 && depre==1) regime = 4; ////// Biphasico
	
	if(facil==1 && depre==1)
	{
        if(sqrt((yampl-ymax)*(yampl-ymax))<tol2) regime = 3; ///// pseudo-linear
	if(sqrt((yant2-ymax)*(yant2-ymax))<tol2) regime = 3; ///// pseudo-linear
	}


	 //~ if(regime==4)	
		 //~ { 
		 //~ if(first_facil==1) regime = 5; ////// Biphasico  facil-depress
		 //~ if(first_depre==1)
			 //~ { 	regime = 6; ////// Biphasico  depress-facil
			 //~ printf("###################################################\n############ Depress-facil###########################\n");
			 //~ }
		 //~ }

if(regime==0)
printf("Frequencia=%f, U=%f, Regime = Nada?\n",freq,U);

if(regime==1)
printf("Frequencia=%f, U=%f, Regime = Facilitation\n",freq,U);

if(regime==2)
printf("Frequencia=%f U=%f, Regime = Depression \n",freq,U);

if(regime==3)
printf("Frequencia=%f U=%f, Regime = Pseudo-linear\n",freq,U);

if(regime==4)
printf("Frequencia=%f U=%f, Regime = Biphasic\n",freq,U);

//~ if(regime==5)
//~ printf("Frequencia=%f U=%f, Regime = BP-facil-depres\n",freq,U);

//~ if(regime==6)
//~ printf("Frequencia=%f U=%f, Regime = BP-depres-facil\n",freq,U);
if(y0==0) regime=-1;
fprintf(s,"%f %f %d\n",freq,U,regime);

//printf("Freq=%f,Tau_rec=%f,Tau_fac=%f, Tau_ina=%f,facil=%d,depre=%d,Freq_medida %f \n",freq,tau_rec,tau_fac,tau_ina,facil, depre, cont_spike*1000.0/nn);

if(y0==0) printf("Alerta, y não sai de zero!!\n");

}
	
  free_dvector(y,1,N*NN+1);
  free_dvector(df,1,N*NN+1);  
  free_dvector(x,1,N*NN+1);
  free_dvector(a,1,N*NN+1);
  free_dvector(b,1,N*NN+1);
  free_dvector(c,1,N*NN+1);

  fclose(o);
  fclose(p);
  fclose(q);
  return 0;
}


void derivs(double y[],double df[]) // Equacoes diferenciais acopladas
{
/////                           x[1] = x ,     x[2] = y,     x[3] = z,     x[4] = u
      df[1]= y[3]/tau_rec;    //// dx/dt   				avaliable
      df[2]= -y[2]/tau_ina;    //// dy/dt  				active
      df[3]= y[2]/tau_ina-y[3]/tau_rec;   //// dz/dt   	recovering
      df[4]= -y[4]/tau_fac;     //// dz/dt 				fraction of avaliable NT using in the next spike
}

float ran1(long *idum)
{
 int j;
 long k;
 static long iy=0;
 static long iv[NTAB];
 float temp;
 
 if(*idum<=0 || !iy)
   {
     if(-(*idum)<1) *idum=1;
     else *idum = -(*idum);
    for(j=NTAB+7;j>=0;j--)
      {
       k=(*idum)/IQ;
       *idum=IA*(*idum-k*IQ)-IR*k;
       if(*idum<0) *idum +=IM;
       if(j<NTAB) iv[j]=*idum;
      }
      iy=iv[0];
   }
   k=(*idum)/IQ;
   *idum=IA*(*idum-k*IQ)-IR*k;
   if(*idum<0) *idum += IM;
   j=iy/NDIV;
   iy=iv[j];
   iv[j]=*idum;
   if((temp=AM*iy)>RNMX) return RNMX;
   else return temp;
}

double *dvector(long nl,long nh)
{
   double *v;
   
   v=(double *)malloc((size_t) ((nh-nl+1+NR_END)*sizeof(double)));
   if (!v) nrerror("allocation failure in dvector()");
   return v-nl+NR_END;
}


void free_dvector(double *v, long nl, long nh)
{
   free((FREE_ARG) (v+nl-NR_END));
}

void nrerror(char error_text[])
{
  fprintf(stderr,"Numerical Recipes run-time error...\n");
  fprintf(stderr,"%s\n",error_text);
  fprintf(stderr,"...now exiting to system...\n");
  exit(1);
}

int *vector(long nl,long nh)
{
   int *v;
   
   v=(int *)malloc((size_t) ((nh-nl+1+NR_END)*sizeof(int)));
   if (!v) nrerror("allocation failure in dvector()");
   return v-nl+NR_END;
}

int **imatrix(long nrl, long nrh, long ncl, long nch)
{
   long i, nrow=nrh-nrl+1,ncol=nch-ncl+1;
   int **m;

   m=(int **) malloc((size_t)((nrow+NR_END)*sizeof(int*)));
   if (!m) nrerror("allocation failure 1 in imatrix()");
   m += NR_END;
   m -= nrl;

   m[nrl]=(int *) malloc((size_t)((nrow*ncol+NR_END)*sizeof(int)));
   if (!m[nrl]) nrerror("allocation failure 2 in imatrix()");
   m[nrl] += NR_END;
   m[nrl] -= ncl;

   for(i=nrl+1;i<=nrh;i++) m[i]=m[i-1]+ncol;

   return m;
}

void free_vector(int *v, long nl, long nh)
{
   free((FREE_ARG) (v+nl-NR_END));
}

void free_imatrix(int **m, long nrl, long nrh, long ncl, long nch)
{
   free((FREE_ARG) (m[nrl]+ncl-NR_END));
   free((FREE_ARG) (m+nrl-NR_END));
}
