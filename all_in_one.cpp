#include <iostream>
#include <iomanip>
#include <math.h>
#include<cstdlib>
#include <ctime>
#include <time.h>
#include <chrono>
using namespace std;

double MAXf(double *c, int n){
	
	double max=c[0];
	for (int i=1;i<n*n; i++)
	{
		if (c[i]>max){
			max=c[i];
		}
	}
	
	return max;
}

double MINf(double *c, int n){
	
	
	double min=c[0];
	for (int i=1;i<n*n; i++)
	{
		if (c[i]<min){
			min=c[i];
		}
	}
	
	
	return min;
}

void dgemm0 (double *a, double *b, double *c, int n) {
	
	
	int i,j,k;
	 for (i=0; i<n; i++) 
	 	for (j=0; j<n; j++){
	 		for (k=0; k<n; k++){

	 			c[i*n+j] += a[i*n+k] * b[k*n+j];
	 			
	 			
	 			
	 			}


	 		
	 	}
				
	
 }

void dgemm1 (double *a, double *b, double *c, int n) {
	
	
	int i,j,k;
	for (i=0; i<n; i++)  
		for (j=0; j<n; j++) { 
			register double r = c[i*n+j] ; 
			for (k=0; k<n; k++) 
				r += a[i*n+k] * b[k*n+j]; 
			c[i*n+j] = r;


	}
	
				
	
 }

void dgemm2 (double *a, double *b, double *c, int n) {
	
	
	int i,j,k;
	for (i=0 ;i<n ;i+=2 )
		for (j=0; j<n; j+=2)
			for (k=0 ; k<n; k+=2){
			
				c[i*n+j]= a[i*n+k] * b[k*n+j] + a[i*n+k+1] * b[(k+1)*n+j] + c[i*n+j];
				
				c[(i+1)*n+j]= a[(i+1)*n+k] * b[k*n+j] + a[(i+1)*n+k+1] * b[(k+1)*n+j] + c[(i+1)*n+j];
				
				c[i*n+(j+1)]= a[i*n+k] * b[k*n+(j+1)] + a[i*n+k+1] * b[(k+1)*n+(j+1)] + c[i*n+(j+1)];
				
				c[(i+1)*n+(j+1)]= a[(i+1)*n+k] * b[k*n+(j+1)] + a[(i+1)*n+k+1] * b[(k+1)*n+(j+1)] + c[(i+1)*n+(j+1)];
			}
				
	
 }
 
 

void dgemm2v2(double *a, double *b, double *c, int n) {
	
	
	int i,j,k;
	for (i=0 ;i<n ;i+=2 )
		for (j=0; j<n; j+=2){
			
			register double rc1=c[i*n+j];
			register double rc2=c[(i+1)*n+j];
			register double rc3=c[i*n+(j+1)];
			register double rc4=c[(i+1)*n+(j+1)];
			
			for (k=0 ; k<n; k+=2){
				
				register double ra1=a[i*n+k];
				register double ra2=a[(i+1)*n+k];
				register double ra3=a[i*n+k+1];
				register double ra4=a[(i+1)*n+k+1];
				
				register double rb1=b[k*n+j];
				register double rb2=b[(k+1)*n+j];
				register double rb3=b[k*n+(j+1)];
				register double rb4=b[(k+1)*n+(j+1)];
			
				rc1= ra1 * rb1 + ra3 * rb2 + rc1;
				
				rc2= ra2 * rb1 + ra4 * rb2 + rc2;
				
				rc3= ra1 * rb3 + ra3 * rb4 + rc3;
				
				rc4= ra2 * rb3 + ra4 * rb4 + rc4;
				
				
			}
			c[i*n+j]=rc1;
			c[(i+1)*n+j]=rc2;
			c[i*n+(j+1)]=rc3;	
			c[(i+1)*n+(j+1)]=rc4;
}
 }



void dgemm3(double *a, double *b, double *c, int n) {
	
	
	int i,j,k;
	for (i=0 ;i<n ;i+=3 )
		for (j=0; j<n; j+=3){
			
			register double rc1= c[(i+0)*n+(j+0)];
			register double rc2= c[(i+1)*n+(j+0)];
			register double rc3= c[(i+2)*n+(j+0)];
			
			register double rc4= c[(i+0)*n+(j+1)];
			register double rc5= c[(i+1)*n+(j+1)];
			register double rc6= c[(i+2)*n+(j+1)];	
			
			register double rc7 =c[(i+0)*n+(j+2)];
			register double rc8 =c[(i+1)*n+(j+2)];
			register double rc9 =c[(i+2)*n+(j+2)];	
			

			
		
			for (k=0 ; k<n; k+=3){

				register double a00=a[(i+0)*n+(k+0)];
				register double a10=a[(i+1)*n+(k+0)];
				register double a20=a[(i+2)*n+(k+0)];
			
				register double b00=b[(k+0)*n+(j+0)];
				register double b01=b[(k+0)*n+(j+1)];
				register double b02=b[(k+0)*n+(j+2)];
			
				rc1 += a00*b00;rc2 += a10*b00;rc3 += a20*b00;
			
				/*.........................................................*/
				
				rc4 += a00*b01; rc5 += a10*b01;	rc6 += a20*b01;
				
				
				/*.........................................................*/
				
				rc7 += a00*b02; rc8 += a10*b02; rc9 += a20*b02;

				a00=a[(i+0)*n+(k+1)];
				a10=a[(i+1)*n+(k+1)];
				a20=a[(i+2)*n+(k+1)];
			
				b00=b[(k+1)*n+(j+0)];
				b01=b[(k+1)*n+(j+1)];
				b02=b[(k+1)*n+(j+2)];

				rc1 += a00*b00;rc2 += a10*b00;rc3 += a20*b00;
			
				
				/*.........................................................*/
				
				rc4 += a00*b01;rc5 += a10*b01;rc6 += a20*b01;
				
				
				/*.........................................................*/
				
				rc7 += a00*b02;rc8 += a10*b02;rc9 += a20*b02;				
				
				a00=a[(i+0)*n+(k+2)];
				a10=a[(i+1)*n+(k+2)];
				a20=a[(i+2)*n+(k+2)];
			
				b00=b[(k+2)*n+(j+0)];
				b01=b[(k+2)*n+(j+1)];
				b02=b[(k+2)*n+(j+2)];
				
				rc1 += a00*b00;rc2 += a10*b00;rc3 += a20*b00;
			
				
				/*.........................................................*/
				
				rc4 += a00*b01;rc5 += a10*b01;rc6 += a20*b01;
				
				
				/*.........................................................*/
				
				rc7 += a00*b02;rc8 += a10*b02;rc9 += a20*b02;		
				
				
				
				
			}
			c[(i+0)*n+(j+0)]=rc1;
			c[(i+1)*n+(j+0)]=rc2;
			c[(i+2)*n+(j+0)]=rc3;
			
			
			c[(i+0)*n+(j+1)]=rc4;
			c[(i+1)*n+(j+1)]=rc5;
		    c[(i+2)*n+(j+1)]=rc6;
			
			
			c[(i+0)*n+(j+2)]=rc7;
			c[(i+1)*n+(j+2)]=rc8;
		    c[(i+2)*n+(j+2)]=rc9;
			
				
        }    
 }




main()
{
	cout<<fixed;
	cout<<setprecision (10);
	
	int i,n,k,j;
	cout<<"please enter the value of n: ";

	cin>>n;
	cout<<"n= "<<n<<endl;

	const int capacity=n*n;

	double *c3;
	double *a;
	double *b;

	double *a2;
	double *b2;
	
	double *a3;
	double *b3;
	
	double *a4;
	double *b4;
	
	double *a5;
	double *b5;
	
	
	
    a=new double[capacity];
	b=new double[capacity];
	c3=new double[capacity];
	
	a2=new double[capacity];
	b2=new double[capacity];

	a3=new double[capacity];
	b3=new double[capacity];

	a4=new double[capacity];
	b4=new double[capacity];
	
	a5=new double[capacity];
	b5=new double[capacity];



	double random1;
	double random2;	
	
	for(int i=0;i<(n*n);i++)
	{	
		random1=(double)rand()/RAND_MAX;
		
		random2=(double)rand()/RAND_MAX;
		a[i]=random1; 
		
		b[i]=random2;
		
		a2[i]=random1;
		
		b2[i]=random2;

		a3[i]=random1;
		
		b3[i]=random2;
		
		a4[i]=random1;
		
		b4[i]=random2;
		
		a5[i]=random1;
		
		b5[i]=random2;
		
		
		

	}
	
	cout<<"dgemm3 results:"<<endl;
	cout<<endl;
	
	
	double t3=clock();
	auto start = std::chrono::steady_clock::now();
	
	dgemm3 (a, b, c3, n);
	
	auto end = std::chrono::steady_clock::now();
	
	double elapsT=double (std::chrono::duration_cast<std::chrono::nanoseconds>(end-start).count());
	
	std::cout<<"Eaplsed time in nanoseconds = "<<elapsT/1000000.0<<std::endl;
	
	
	double tt3=(clock()-t3)/CLOCKS_PER_SEC;
	cout<<"Time = "<<tt3<<" Sec"<<endl; 
	
	
	double max3=MAXf(c3,n);
	double min3=MINf(c3,n);
	
	
	
	cout<<"GFlops = "<<(8.0*16.0*n/4.0*n/4.0*n/4.0)/(elapsT)<<endl;
	cout<<"max-min = "<<max3-min3<<endl;
	
	
	
	cout<<"---------------------------------------------------------------------------------------------------------------------"<<endl;
	cout<<"dgemm2v2 results:"<<endl;

	
	cout<<endl;
	double *c2v2;
	c2v2=new double[capacity];
	
	
	double t2v2=clock();
	auto start2v2 = std::chrono::steady_clock::now();
	
	dgemm2v2 (a2, b2, c2v2, n);
	
	auto end2v2 = std::chrono::steady_clock::now();
	
	double elapsT2v2=double (std::chrono::duration_cast<std::chrono::nanoseconds>(end2v2-start2v2).count());
	std::cout<<"Eaplsed time in milliseconds = "<<elapsT2v2/1000000.0<<std::endl;
	
	
	double tt2v2=(clock()-t2v2)/CLOCKS_PER_SEC;
	cout<<"Time = "<<tt2v2<<" Sec"<<endl; 
	
	double max2v2=MAXf(c2v2,n);
	double min2v2=MINf(c2v2,n);	
	
	cout<<"GFlops = "<<(16.0*n/2.0*n/2.0*n/2.0)/(elapsT2v2)<<endl;
	cout<<"max-min = "<<max2v2-min2v2<<endl;
	
	cout<<"---------------------------------------------------------------------------------------------------------------------"<<endl;
	cout<<"dgemm2 results:"<<endl;
	
	cout<<endl;
	double *c2;
	c2=new double[capacity];
	

	
	double t2=clock();
	
	auto start2 = std::chrono::steady_clock::now();
	dgemm2 (a3, b3, c2, n);
	
	auto end2 = std::chrono::steady_clock::now();
	
	double elapsT2=double (std::chrono::duration_cast<std::chrono::nanoseconds>(end2-start2).count());
	std::cout<<"Eaplsed time in milliseconds = "<<elapsT2/1000000.0<<std::endl;

	
	double tt2=(clock()-t2)/CLOCKS_PER_SEC;
	cout<<"Time = "<<tt2<<" Sec"<<endl; 
	
	double max2=MAXf(c2,n);
	double min2=MINf(c2,n);	
	
	cout<<"GFlops = "<<(16.0*n/2.0*n/2.0*n/2.0)/(elapsT2)<<endl;
	cout<<"max-min = "<<max2-min2<<endl;
	
	cout<<"---------------------------------------------------------------------------------------------------------------------"<<endl;
	cout<<"dgemm1 results:"<<endl;
	
	cout<<endl;
	double *c1;
	c1=new double[capacity];
	

	
	double t1=clock();
	
	auto start1 = std::chrono::steady_clock::now();
	dgemm1 (a4, b4, c1, n);
	
	auto end1 = std::chrono::steady_clock::now();
	
	double elapsT1=double (std::chrono::duration_cast<std::chrono::nanoseconds>(end1-start1).count());
	std::cout<<"Eaplsed time in miliseconds = "<<elapsT1/1000000.0<<std::endl;
	
	double tt1=(clock()-t1)/CLOCKS_PER_SEC;
	cout<<"Time = "<<tt1<<" Sec"<<endl; 
	
	double max1=MAXf(c1,n);
	double min1=MINf(c1,n);	
	
	cout<<"GFlops = "<<(2.0*n*n*n/(elapsT1))<<endl;
	cout<<"max-min = "<<max1-min1<<endl;
	
	cout<<"---------------------------------------------------------------------------------------------------------------------"<<endl;
	cout<<"dgemm0 results:"<<endl;
	
	cout<<endl;
	double *c0;
	c0=new double[capacity];
	

	
	double t0=clock();
	auto start0 = std::chrono::steady_clock::now();
	dgemm0 (a5, b5, c0, n);
	
	auto end0 = std::chrono::steady_clock::now();
	
	double elapsT0=double (std::chrono::duration_cast<std::chrono::nanoseconds>(end0-start0).count());
	std::cout<<"Eaplsed time in miliseconds = "<<elapsT0/1000000.0<<std::endl;

	double tt0=(clock()-t0)/CLOCKS_PER_SEC;
	cout<<"Time = "<<tt0<<" Sec"<<endl; 
	
	double max0=MAXf(c0,n);
	double min0=MINf(c0,n);	
	
	cout<<"GFlops = "<<(2.0*n*n*n/(elapsT0))<<endl;
	cout<<"max-min = "<<max0-min0<<endl;
	cout<<"---------------------------------------------------------------------------------------------------------------------"<<endl;
	


	double *diff1;
	double *diff2;
	double *diff3;
	

    diff1=new double[capacity];
    diff2=new double[capacity];
    diff3=new double[capacity];

	cout<<"Correctness Checking:"<<endl;
	cout<<endl;
	
	for (int u=0;u<n*n;u++){
		diff1[u]=c1[u]-c0[u];
		diff2[u]=c2[u]-c1[u];
		diff3[u]=c3[u]-c2[u];
		
	}

	
	
	cout<<"max difference of all elements = ";
	cout<<setprecision (6);
	cout<<MAXf(diff1,n)<<endl;
	
	cout<<"max difference of all elements = ";
	cout<<setprecision (6);
	cout<<MAXf(diff2,n)<<endl;
	
	cout<<"max difference of all elements = ";
	cout<<setprecision (6);
	cout<<MAXf(diff3,n)<<endl;

	system("pause;");
	return 0;
  	
}





