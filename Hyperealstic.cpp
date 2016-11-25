#include "stdafx.h"		
#include"CONFIG.h"	//CONF.cppはインクルードしなくても、CONFIG.hをインクルードするだけでよい
#include"PART.h"		//class PART定義
#include<vector>
#include"function.h"
#include "Hyperelastic.h"

void calc_hyper(mpsconfig *CON,vector<mpselastic>&PART,vector<hyperelastic>&HYPER,vector<hyperelastic2>&HYPER1,int hyper_number,int t)
{
	double le=CON->get_distancebp();
	double Dt=CON->get_dt()*CON->get_interval();
	double V=get_volume(CON);
	double mi=V*CON->get_hyper_density();
	double r=19*le;
	double c10=CON->get_c10();
	double c01=CON->get_c01();
	int repetation=0;

	if(t==1)
	{
		for(int i=0;i<hyper_number;i++)
		{
			for(int D=0;D<DIMENSION;D++)
			{
				PART[i].q0[D]=0;
				PART[i].q0[D]=PART[i].r[D];
			}
		}
		calc_initial_p(PART,HYPER,hyper_number,mi,le,Dt);
		calc_neighbor(PART,HYPER,hyper_number,r);
	}

	calc_DgDq(PART,HYPER,HYPER1,hyper_number,r,V);
	cout<<"DgDq is calculated."<<endl;

	calc_lambda(PART,HYPER,HYPER1,hyper_number,Dt,mi);
	cout<<"Lambda is calculated."<<endl;

	calc_stress(PART,HYPER,hyper_number,r,c10,c01);
	cout<<"Stress is calculated."<<endl;

	calc_half_p(PART,HYPER,HYPER1,hyper_number,Dt,repetation);
	repetation++;
	cout<<"Half momentum is calculated."<<endl;


	renew_q(PART,HYPER,hyper_number,t,Dt,mi);
	cout<<"Position is renewed."<<endl;

	calc_DgDq(PART,HYPER,HYPER1,hyper_number,r,V);
	cout<<"DgDq is renewed."<<endl;

	calc_stress(PART,HYPER,hyper_number,r,c10,c01);
	cout<<"Stress is renewed."<<endl;

	calc_differential_p(PART,HYPER,HYPER1,hyper_number,Dt);
	cout<<"Differential_p is calculated."<<endl;

	renew_lambda(PART,HYPER,HYPER1,hyper_number,Dt,mi);
	cout<<"Lambda is renewed."<<endl;
	
	calc_half_p(PART,HYPER,HYPER1,hyper_number,Dt,repetation);
	cout<<"Half momentum is renewed."<<endl;
	
	
}

void calc_neighbor(vector<mpselastic>&PART,vector<hyperelastic>&HYPER,int particle_number,double r)
{
	for(int i=0;i<particle_number;i++)
	{
		HYPER[i].N=0;
		for(int j=0;j<200;j++)
		{
			HYPER[i].NEI[j]=0;
		}
	}

	double dis=0;
	int N;
	for(int i=0;i<particle_number;i++)
	{
		N=0;
		for(int j=0;j<particle_number;j++)
		{
			dis=sqrt((PART[i].r[A_X]-PART[j].r[A_X])*(PART[i].r[A_X]-PART[j].r[A_X])+(PART[i].r[A_Y]-PART[j].r[A_Y])*(PART[i].r[A_Y]-PART[j].r[A_Y])+(PART[i].r[A_Z]-PART[j].r[A_Z])*(PART[i].r[A_Z]-PART[j].r[A_Z]));
			if(dis<r && i!=j)
			{
				HYPER[i].NEI[N]=j;
				N++;
			}
		}
		HYPER[i].N=N;
	}
}


void calc_lambda(vector<mpselastic> &PART,vector<hyperelastic>&HYPER,vector<hyperelastic2>&HYPER1,int particle_number,double Dt,double mk)
{
	for(int i=0;i<particle_number;i++)	HYPER[i].lambda=0;

	double DFDlambda;

	vector<vector<double> >DFDlam(particle_number,vector<double>(particle_number));
	vector<vector<double> >DgDq_X(particle_number,vector<double>(particle_number));
	vector<vector<double> >DgDq_Y(particle_number,vector<double>(particle_number));
	vector<vector<double> >DgDq_Z(particle_number,vector<double>(particle_number));

	for(int i=0;i<particle_number;i++)
	{
		for(int j=0;j<particle_number;j++)
		{
			for(int k=0;k<particle_number;k++)
			{
				for(int l=0;l<particle_number;l++)
				{
					DgDq_X[k][l]=0;
					DgDq_Y[k][l]=0;
					DgDq_Z[k][l]=0;
				}
			}

			DFDlambda=0;
			for(int k=0;k<particle_number;k++)
			{
				DgDq_X[i][k]+=HYPER1[j*particle_number+i].DgDq[A_X];
				DgDq_Y[i][k]+=HYPER1[j*particle_number+i].DgDq[A_Y];
				DgDq_Z[i][k]+=HYPER1[j*particle_number+i].DgDq[A_Z];
				DgDq_X[j][k]+=HYPER1[j*particle_number+i].DgDq[A_X];
				DgDq_Y[j][k]+=HYPER1[j*particle_number+i].DgDq[A_Y];
				DgDq_Z[j][k]+=HYPER1[j*particle_number+i].DgDq[A_Z];
				DFDlambda+=(-1)*Dt*Dt/2/mk*(DgDq_X[i][k]*DgDq_X[j][k]+DgDq_Y[i][k]*DgDq_Y[j][k]+DgDq_Z[i][k]*DgDq_Z[j][k]);
			}
			DFDlam[i][j]=DFDlambda;
		}//jに関するfor文の終わり
	}//iに関するfor文の終わり

	//lambdaを求める
	double *B=new double [particle_number];
	double *matrix=new double [particle_number*particle_number];

	for(int i=0;i<particle_number;i++)	for(int j=0;j<particle_number;j++)	matrix[i*particle_number+j]=DFDlam[i][j];
	
	for(int i=0;i<particle_number;i++)	B[i]=0;

	gauss_sidel(matrix,B,particle_number);
	
	for(int i=0;i<particle_number;i++)
	{	
		HYPER[i].lambda=B[i];
		cout<<"HYPER["<<i<<"].lambda"<<HYPER[i].lambda<<endl;
	}



	for(int i=0;i<particle_number;i++)
	{
		DFDlam[i].clear();
		DgDq_X[i].clear();
		DgDq_Y[i].clear();
		DgDq_Z[i].clear();
	}

	DFDlam.clear();
	DgDq_X.clear();
	DgDq_Y.clear();
	DgDq_Z.clear();

	delete[] B;
	delete[] matrix;
}



void calc_initial_p(vector<mpselastic>&PART,vector<hyperelastic>&HYPER,int particle_number,double mi,double le,double Dt)
{
	int t=30;
	int b=2;
	for(int i=0;i<particle_number;i++)	for(int D=0;D<DIMENSION;D++)	HYPER[i].p[D]=0;
	for(int i=0;i<particle_number;i++)	
	{
/*		HYPER[i].p[A_Z]=-9.8*mi*PART[i].r[A_Z]*Dt;
		HYPER[i].p[A_X]=0;
		HYPER[i].p[A_Y]=0;
*/		HYPER[i].p[A_X]=-t*mi*(PART[i].r[A_Z]/(9*le))*(PART[i].r[A_Z]/(9*le))*(PART[i].r[A_Z]/(9*le))*PART[i].r[A_Y]+b*(3*(PART[i].r[A_Z]/(9*le))*(PART[i].r[A_Z]/(9*le))-1)+mi*(0.4*PART[i].r[A_X]+1.0);
		HYPER[i].p[A_Y]=t*mi*(PART[i].r[A_Z]/(9*le))*(PART[i].r[A_Z]/(9*le))*(PART[i].r[A_Z]/(9*le))*PART[i].r[A_X]+mi*0.4*PART[i].r[A_Y];
		HYPER[i].p[A_Z]=0+mi*0.4*PART[i].r[A_Z];
		cout<<"HYPER["<<i<<"].p[D]="<<HYPER[i].p[A_X]<<","<<HYPER[i].p[A_Y]<<","<<HYPER[i].p[A_Z]<<endl;
	}
}



void calc_half_p(vector<mpselastic>&PART,vector<hyperelastic>&HYPER,vector<hyperelastic2>&HYPER1,int particle_number,double Dt,int repetation)
{
	int D;
	int D2;

	double partial_half_p[3];
	double E[3][3];

	for(D=0;D<DIMENSION;D++)
	{
		for(D2=0;D2<DIMENSION;D2++)
		{
			if(D==D2)	E[D][D2]=1;
			else
			{
				E[D][D2]=0;
			}
		}
	}

	for(int i=0;i<particle_number;i++)
	{
		for(D=0;D<DIMENSION;D++)	partial_half_p[D]=0;
		for(int j=0;j<particle_number;j++)
		{
			for(D=0;D<DIMENSION;D++)	for(D2=0;D2<DIMENSION;D2++)	partial_half_p[D]+=Dt/2*(HYPER[j].stress[D][D2]-HYPER[j].lambda*E[D][D2])*HYPER1[j*particle_number+i].DgDq[D];
		}//jに関するfor文の終わり

		if(repetation%2==0)
		{
			for(int k=0;k<particle_number;k++)	for(D=0;D<DIMENSION;D++)	HYPER[k].half_p[D]=0;
			for(D=0;D<DIMENSION;D++)	HYPER[i].half_p[D]=HYPER[i].p[D]+partial_half_p[D];
		}

		else
		{
			for(D=0;D<DIMENSION;D++)	HYPER[i].p[D]=HYPER[i].half_p[D]+partial_half_p[D];
		}
	}//iに関するfor文の終わり
}

void calc_stress(vector<mpselastic>&PART,vector<hyperelastic>&HYPER,int particle_number,double r,double c10,double c01)
{
	int D;
	int D2;

	for(int i=0;i<particle_number;i++)	for(D=0;D<DIMENSION;D++)	for(D2=0;D2<DIMENSION;D2++)	HYPER[i].stress[D][D2]=0;

	double ajjn[3];
	double qjjn[3];
	double dis=0;
	double wjjn=0;
	double fj[3][3];
	double Aj[3][3];
	double inverse_Aj[3][3];
	double Fj[3][3];
	double transposed_Fj[3][3];
	double b[3][3];
	double bb[3][3];
	double E[3][3];
	double trace_b;
	double trace_bb;
	double transposed_M[3][3];
	double inverse_M[3][3];
	double transposed_inverseM[3][3];


	for(D=0;D<DIMENSION;D++)
	{
		ajjn[D]=0;
		qjjn[D]=0;
		
		for(D2=0;D2<DIMENSION;D2++)
		{
			inverse_Aj[D][D2]=0;
			transposed_Fj[D][D2]=0;
			transposed_M[D][D2]=0;
			inverse_M[D][D2]=0;
			transposed_inverseM[D][D2]=0;
			if(D==D2)	E[D][D2]=1;
			else
			{
				E[D][D2]=0;
			}
		}
	}

	for(int j=0;j<particle_number;j++)
	{
		for(D=0;D<DIMENSION;D++)
		{
			for(D2=0;D2<DIMENSION;D2++)
			{
				Aj[D][D2]=0;
				fj[D][D2]=0;
			}
		}
		for(int jn=0;jn<HYPER[j].N;jn++)	//粒子番号jに関するAとFを求める *Aiと同様
		{
			for(D=0;D<DIMENSION;D++)
			{
				ajjn[D]=PART[HYPER[j].NEI[jn]].q0[D]-PART[j].q0[D];
				qjjn[D]=PART[HYPER[j].NEI[jn]].r[D]-PART[j].r[D];

			}
			dis=sqrt(ajjn[A_X]*ajjn[A_X]+ajjn[A_Y]*ajjn[A_Y]+ajjn[A_Z]*ajjn[A_Z]);
			wjjn=kernel4(r,dis);

			for(D=0;D<DIMENSION;D++)	//転置＆逆行列化するためにAjを格納	//Fjの一部の計算
			{
				for(D2=0;D2<DIMENSION;D2++)
				{
					Aj[D][D2]+=wjjn*ajjn[D]*ajjn[D2];
					fj[D][D2]+=wjjn*qjjn[D]*ajjn[D2];
				}
			}
		}			

		calc_transposed_inverse_matrix(Aj,transposed_M,inverse_Aj,transposed_inverseM,true,true);	//転置＆逆行列したAjを求める

		for(D=0;D<DIMENSION;D++)	for(D2=0;D2<DIMENSION;D2++)	Fj[D][D2]=0;
		for(D=0;D<DIMENSION;D++)	for(D2=0;D2<DIMENSION;D2++)	for(int D3=0;D3<DIMENSION;D3++)	Fj[D][D2]+=fj[D][D3]*inverse_Aj[D3][D2];	//Fj, 転置＆逆行列したFjを求める
		calc_transposed_inverse_matrix(Fj,transposed_Fj,inverse_M,transposed_inverseM,true,false);


		for(D=0;D<DIMENSION;D++)
		{
			for(D2=0;D2<DIMENSION;D2++)
			{
				b[D][D2]=0;
				bb[D][D2]=0;
			}
		}

		for(D=0;D<DIMENSION;D++)
		{
			for(D2=0;D2<DIMENSION;D2++)
			{
				for(int D3=0;D3<DIMENSION;D3++)
				{
					b[D][D2]+=Fj[D][D3]*transposed_Fj[D3][D2];
					bb[D][D2]+=b[D][D3]*b[D3][D2];
				}
			}
		}

		trace_b=0;
		trace_bb=0;

		for(D=0;D<DIMENSION;D++)
		{
			trace_b+=b[D][D];
			trace_bb+=bb[D][D];
		}

		for(D=0;D<DIMENSION;D++)
		{
			for(D2=0;D2<DIMENSION;D2++)
			{
				HYPER[j].stress[D][D2]=2*((c10+c01*trace_b)*(b[D][D2]-1/3*trace_b*E[D][D2])-c01*(bb[D][D2]-1/3*trace_bb*E[D][D2]));
			}
		}
	}
}


void renew_q(vector<mpselastic>&PART,vector<hyperelastic>&HYPER,int particle_number,int t,double Dt,double mi)
{
	for(int i=0;i<particle_number;i++)	for(int D=0;D<DIMENSION;D++)	PART[i].r[D]=PART[i].r[D]+Dt*HYPER[i].half_p[D]/mi;
}

void calc_differential_p(vector<mpselastic>&PART,vector<hyperelastic>&HYPER,vector<hyperelastic2>&HYPER1,int particle_number,double Dt)
{

	int D;
	int D2;

	for(int i=0;i<particle_number;i++)	for(D=0;D<DIMENSION;D++)	HYPER[i].differential_p[D]=0;

	double partial_differential_p[3];


	double E[3][3];
	for(D=0;D<DIMENSION;D++)
	{
		for(D2=0;D2<DIMENSION;D2++)
		{
			if(D==D2)E[D][D2]=1;
			else
			{
				E[D][D2]=0;
			}
		}
	}

	for(int i=0;i<particle_number;i++)
	{
		for(D=0;D<DIMENSION;D++)	partial_differential_p[D]=0;
		for(int j=0;j<particle_number;j++)
		{
			for(D=0;D<DIMENSION;D++)	for(D2=0;D2<DIMENSION;D2++)	partial_differential_p[D]+=HYPER[j].stress[D][D2]*HYPER1[j*particle_number+i].DgDq[D];
		}
		for(D=0;D<DIMENSION;D++)	HYPER[i].differential_p[D]=HYPER[i].half_p[D]+Dt/2*partial_differential_p[D];
	}
}



void renew_lambda(vector<mpselastic>&PART,vector<hyperelastic>&HYPER,vector<hyperelastic2>&HYPER1,int particle_number,double Dt,double mk)
{
	for(int i=0;i<particle_number;i++)	HYPER[i].lambda=0;

	double N_left;
	double N_right;
	vector<vector<double> > N_Left(particle_number,vector<double>(particle_number));
	vector<double> N_Right(particle_number);

	vector<vector<double > >DgDq_X(particle_number,vector<double>(particle_number));
	vector<vector<double> >DgDq_Y(particle_number,vector<double>(particle_number));
	vector<vector<double> > DgDq_Z(particle_number,vector<double>(particle_number));

	for(int i=0;i<particle_number;i++)
	{
		for(int j=0;j<particle_number;j++)
		{
			for(int k=0;k<particle_number;k++)
			{
				for(int l=0;l<particle_number;l++)
				{
					DgDq_X[k][l]=0;
					DgDq_Y[k][l]=0;
					DgDq_Z[k][l]=0;
				}
			}

			N_left=0;
			N_right=0;
			for(int k=0;k<particle_number;k++)
			{
				DgDq_X[i][k]+=HYPER1[j*particle_number+i].DgDq[A_X];
				DgDq_Y[i][k]+=HYPER1[j*particle_number+i].DgDq[A_Y];
				DgDq_Z[i][k]+=HYPER1[j*particle_number+i].DgDq[A_Z];
				DgDq_X[j][k]+=HYPER1[j*particle_number+i].DgDq[A_X];
				DgDq_Y[j][k]+=HYPER1[j*particle_number+i].DgDq[A_Y];
				DgDq_Z[j][k]+=HYPER1[j*particle_number+i].DgDq[A_Z];

				N_left+=Dt/2/mk*(DgDq_X[i][k]*DgDq_X[j][k]+DgDq_Y[i][k]*DgDq_Y[j][k]+DgDq_Z[i][k]*DgDq_Z[j][k]);
				N_right+=1/mk*(DgDq_X[i][k]*HYPER[k].differential_p[A_X]+DgDq_Y[i][k]*HYPER[k].differential_p[A_Y]+DgDq_Z[i][k]*HYPER[k].differential_p[A_Z]);
			}
			N_Right[i]=N_right;			
			N_Left[j][i]=N_left;
		}//jに関するfor文の終わり
	}//iに関するfor文の終わり


	//lambdaを求める
	double *B=new double [particle_number];
	double *matrix2=new double [particle_number*particle_number];


	for(int i=0;i<particle_number;i++)
	{
		B[i]=N_Right[i];
		for(int j=0;j<particle_number;j++)
		{
			matrix2[i*particle_number+j]=N_Left[i][j];
		}
	}

	gauss_sidel(matrix2,B,particle_number);

	for(int i=0;i<particle_number;i++)
	{
		HYPER[i].lambda=B[i];
		cout<<"HYPER["<<i<<"].lambda"<<HYPER[i].lambda<<endl;
	}

	for(int i=0;i<particle_number;i++)
	{
		DgDq_X[i].clear();
		DgDq_Y[i].clear();
		DgDq_Z[i].clear();
		N_Left[i].clear();
	}
	N_Left.clear();
	N_Right.clear();
	delete[] matrix2;
	delete[] B;
	DgDq_X.clear();
	DgDq_Y.clear();
	DgDq_Z.clear();

}


void calc_DgDq(vector<mpselastic> &PART,vector<hyperelastic>&HYPER,vector<hyperelastic2>&HYPER1,int particle_number,double r,double V)
{

	int D;
	int D2;
	int D3;
	double aiin[3];
	double ajjn[3];
	double qjjn[3];
	double aij[3];
	double dis=0;

	double Ai[3][3];
	double transposed_inverse_Ai[3][3];
	double Aj[3][3];
	double inverse_Aj[3][3];
	double transposed_inverse_Aj[3][3];

	double fj[3][3];
	double Fj[3][3];
	double transposed_inverse_Fj[3][3];

	double transposed_M[3][3];
	double inverse_M[3][3];

	double partial_n0ij[3];
	double n0ij[3];
	double partial_DgDq[3];

	double wiin=0;
	double wjjn=0;
	double wij=0;

	for(int i=0;i<particle_number*particle_number;i++)	for(D=0;D<DIMENSION;D++)	HYPER1[i].DgDq[D]=0;

	for(D=0;D<DIMENSION;D++)
	{
		aiin[D]=0;
		ajjn[D]=0;
		qjjn[D]=0;
		aij[D]=0;
		n0ij[D]=0;
		for(D2=0;D2<DIMENSION;D2++)
		{
			transposed_inverse_Ai[D][D2]=0;
			inverse_Aj[D][D2]=0;
			transposed_inverse_Aj[D][D2]=0;
			Fj[D][D2]=0;
			transposed_inverse_Fj[D][D2]=0;
			transposed_M[D][D2]=0;
			inverse_M[D][D2]=0;
		}
	}

	for(int j=0;j<particle_number;j++)
	{
		for(int i=0;i<particle_number;i++)
		{
			for(D=0;D<DIMENSION;D++)	for(D2=0;D2<DIMENSION;D2++)	Ai[D][D2]=0;
			for(D=0;D<DIMENSION;D++)	partial_n0ij[D]=0;

			//粒子番号iに関するAを求める
			for(int in=0;in<HYPER[i].N;in++)
			{
				//aiin, qiin, wiinを求める
				for(D=0;D<DIMENSION;D++)	aiin[D]=PART[HYPER[i].NEI[in]].q0[D]-PART[i].q0[D];
				dis=sqrt(aiin[A_X]*aiin[A_X]+aiin[A_Y]*aiin[A_Y]+aiin[A_Z]*aiin[A_Z]);
				wiin=kernel4(r,dis);

				//転置＆逆行列化するためにAiを格納
				for(D=0;D<DIMENSION;D++)
				{
					for(D2=0;D2<DIMENSION;D2++)	Ai[D][D2]+=wiin*aiin[D]*aiin[D2];
					partial_n0ij[D]+=wiin*aiin[D];
				}
			}

			//転置＆逆行列したAiを求める
			calc_transposed_inverse_matrix(Ai,transposed_M,inverse_M,transposed_inverse_Ai,true,true);

			//粒子番号jに関するAとFを求める *Aiと同様
			for(D=0;D<DIMENSION;D++)
			{
				for(D2=0;D2<DIMENSION;D2++)
				{
					Aj[D][D2]=0;
					fj[D][D2]=0;
				}
			}

			for(int jn=0;jn<HYPER[j].N;jn++)
			{
				for(D=0;D<DIMENSION;D++)
				{
					ajjn[D]=PART[HYPER[j].NEI[jn]].q0[D]-PART[j].q0[D];
					qjjn[D]=PART[HYPER[j].NEI[jn]].r[D]-PART[j].r[D];
				}
				dis=sqrt(ajjn[A_X]*ajjn[A_X]+ajjn[A_Y]*ajjn[A_Y]+ajjn[A_Z]*ajjn[A_Z]);
				wjjn=kernel4(r,dis);

				//転置＆逆行列化するためにAjを格納
				//Fjの一部の計算

				for(D=0;D<DIMENSION;D++)
				{
					for(D2=0;D2<DIMENSION;D2++)
					{
						Aj[D][D2]+=wjjn*ajjn[D]*ajjn[D2];
						fj[D][D2]+=wjjn*qjjn[D]*ajjn[D2];
					}
				}
			}

			//転置＆逆行列したAjを求める
			calc_transposed_inverse_matrix(Aj,transposed_M,inverse_Aj,transposed_inverse_Aj,true,true);

			//Fj, 転置＆逆行列したFjを求める
			for(D=0;D<DIMENSION;D++)	for(D2=0;D2<DIMENSION;D2++)	Fj[D][D2]=0;
			for(D=0;D<DIMENSION;D++)	for(D2=0;D2<DIMENSION;D2++)	for(D3=0;D3<DIMENSION;D3++)	Fj[D][D2]+=fj[D][D3]*inverse_Aj[D3][D2];

			calc_transposed_inverse_matrix(Fj,transposed_M,inverse_M,transposed_inverse_Fj,true,true);

			if(i!=j)
			{
				for(D=0;D<DIMENSION;D++)	aij[D]=PART[j].q0[D]-PART[i].q0[D];
					dis=sqrt(aij[A_X]*aij[A_X]+aij[A_Y]*aij[A_Y]+aij[A_Z]*aij[A_Z]);
				wij=kernel4(r,dis);
				for(D=0;D<DIMENSION;D++)	n0ij[D]=0;
				for( D=0;D<DIMENSION;D++)	for(D2=0;D2<DIMENSION;D2++)	n0ij[D]+=V*transposed_inverse_Aj[D][D2]*wij*aij[D];
			}
			
			if(i==j)
			{
				for(D=0;D<DIMENSION;D++)	n0ij[D]=0;
				for(D=0;D<DIMENSION;D++)	for(D2=0;D2<DIMENSION;D2++)	n0ij[D]+=V*transposed_inverse_Ai[D][D2]*partial_n0ij[D];
			}
			for(D=0;D<DIMENSION;D++)	partial_DgDq[D]=0;
			for(D=0;D<DIMENSION;D++)	for(D2=0;D2<DIMENSION;D2++)	partial_DgDq[D]+=transposed_inverse_Fj[D][D2]*n0ij[D];

			for(D=0;D<DIMENSION;D++)	HYPER1[j*particle_number+i].DgDq[D]=partial_DgDq[D];
		}
	}
}



	//逆行列と転置行列を求める関数
void calc_transposed_inverse_matrix(double M[3][3],double transposed_M[3][3],double inverse_M[3][3],double transposed_inverseM[3][3],bool transport,bool inversion)
{

	int D;
	int D2;
	double t_inverse_M[3][3];
	double buf[3];

	for(D=0;D<DIMENSION;D++)
	{
		for(D2=0;D2<DIMENSION;D2++)
		{
			transposed_M[D][D2]=0;
			transposed_inverseM[D][D2]=0;
		}
	}

	if(transport==true)
	{
		for(D=0;D<DIMENSION;D++)	for(D2=0;D2<DIMENSION;D2++)	transposed_M[D][D2]=M[D2][D];
	}

	if(inversion==true)
	{
		for(D=0;D<DIMENSION;D++)
		{
			buf[D]=0;
			for(D2=0;D2<DIMENSION;D2++)
			{
				if(D==D2)	inverse_M[D][D2]=1;
				else inverse_M[D][D2]=0;
			}
		}
		
		//about [0][D]
		buf[0]=1/M[0][0];
		for(D=0;D<DIMENSION;D++)
		{
			M[0][D]=M[0][D]*buf[0];
			inverse_M[0][D]=inverse_M[0][D]*buf[0];
			M[1][D]-=M[0][D]*M[1][0];
			M[2][D]-=M[0][D]*M[2][0];
			inverse_M[1][D]-=inverse_M[0][D]*M[1][0];
			inverse_M[2][D]-=inverse_M[0][D]*M[2][0];
		}

		//about [1][D]
		buf[1]=1/M[1][1];
		for(D=0;D<DIMENSION;D++)
		{
			M[1][D]=M[1][D]*buf[1];
			inverse_M[1][D]=inverse_M[1][D]*buf[1];
			M[0][D]-=M[1][D]*M[0][1];
			M[2][D]-=M[1][D]*M[2][1];
			inverse_M[0][D]-=inverse_M[1][D]*M[0][1];
			inverse_M[2][D]-=inverse_M[1][D]*M[2][1];
		}

		//about [2][D]
		buf[2]=1/M[2][2];
		for(D=0;D<DIMENSION;D++)
		{
			M[2][D]=M[2][D]*buf[2];
			inverse_M[2][D]=inverse_M[2][D]*buf[2];
			M[0][D]-=M[2][D]*M[0][2];
			M[1][D]-=M[2][D]*M[1][2];
			inverse_M[0][D]-=inverse_M[2][D]*M[0][2];
			inverse_M[1][D]-=inverse_M[2][D]*M[1][2];
		}
	}


	if((transport==true)&(inversion==true))
	{
		for(D=0;D<DIMENSION;D++)	for(D2=0;D2<DIMENSION;D2++)	M[D][D2]=transposed_M[D][D2];

		for(D=0;D<DIMENSION;D++)
		{
			buf[D]=0;
			for(D2=0;D2<DIMENSION;D2++)
			{
				if(D==D2)	t_inverse_M[D][D2]=1;
				else t_inverse_M[D][D2]=0;
			}
		}

		//about [0][D]
		buf[0]=1/M[0][0];
		for(D=0;D<DIMENSION;D++)
		{
			M[0][D]=M[0][D]*buf[0];
			t_inverse_M[0][D]=t_inverse_M[0][D]*buf[0];
			M[1][D]-=M[0][D]*M[1][0];
			M[2][D]-=M[0][D]*M[2][0];
			t_inverse_M[1][D]-=t_inverse_M[0][D]*M[1][0];
			t_inverse_M[2][D]-=t_inverse_M[0][D]*M[2][0];
		}

		//about [1][D]
		buf[1]=1/M[1][1];
		for(D=0;D<DIMENSION;D++)
		{
			M[1][D]=M[1][D]*buf[1];
			t_inverse_M[1][D]=t_inverse_M[1][D]*buf[1];
			M[0][D]-=M[1][D]*M[0][1];
			M[2][D]-=M[1][D]*M[2][1];
			t_inverse_M[0][D]-=t_inverse_M[1][D]*M[0][1];
			t_inverse_M[2][D]-=t_inverse_M[1][D]*M[2][1];
		}

		//about [2][D]
		buf[2]=1/M[2][2];
		for(D=0;D<DIMENSION;D++)
		{
			M[2][D]=M[2][D]*buf[2];
			t_inverse_M[2][D]=t_inverse_M[2][D]*buf[2];
			M[0][D]-=M[2][D]*M[0][2];
			M[1][D]-=M[2][D]*M[1][2];
			t_inverse_M[0][D]-=t_inverse_M[2][D]*M[0][2];
			t_inverse_M[1][D]-=t_inverse_M[2][D]*M[1][2];
		}

		for(D=0;D<DIMENSION;D++)	for(D2=0;D2<DIMENSION;D2++)	transposed_inverseM[D][D2]=t_inverse_M[D][D2];
	}	
}

void gauss_sidel(double *matrix,double *B,int N)
{
	double aii=0;
	double *BB=new double [N];
	for(int i=0;i<N;i++)
	{
		BB[i]=0;
		for(int j=0;j<N;j++)
		{
			aii=matrix[i*N+i];
			BB[i]+=matrix[j*N+i]*B[i];
		}
		BB[i]=(B[i]-BB[i])/aii;
	}
	for(int i=0;i<N;i++)	B[i]=BB[i];
	delete[] BB;
}


void calc_trace(double M[3][3],double det)
{
	det=0.0;
	det=M[0][0]+M[1][1]+M[2][2];
	det+=M[0][1]*M[1][2]*M[2][0];
	det+=M[1][0]*M[2][1]*M[0][2];
	det-=M[0][0]*M[1][2]*M[2][1];
	det-=M[0][2]*M[1][1]*M[2][0];
	det-=M[0][1]*M[1][0]*M[2][2];
}

//detを求める関数　※自己流のため自信なし
void calc_trace(double **M,double det,int particle_number)
{
	double det_plus=0;
	double det_minus=0;

	for(int i=0;i<particle_number;i++)
	{
		for(int j=0;j<particle_number;j++)
		{
			det_plus*=M[i%particle_number][j%particle_number];
			i=i+1;
			j=j+1;

			det_minus*=M[i%particle_number][j%particle_number];
			i=i-1;
			j=j-1;
		}
		det+=det_plus-det_minus;
	}
}



void calc_jacobi(vector<vector<double> > matrix,vector<double> answer,int particle_number)
{
	double *mean_x= new double [particle_number];
	double *x=new double [particle_number];

	for(int k=0;k<particle_number;k++)
	{
		for(int i=0;i<particle_number;i++)
		{
			mean_x[i]=0;
			for(int j=0;j<particle_number;j++)	mean_x[i]=mean_x[i]+matrix[i][j]*x[j];
			for(int j=0;j<particle_number;j++)	mean_x[i]=-mean_x[i]/matrix[i][j];
		}
		x[k]=mean_x[k];
	}

	delete[] mean_x;
	delete[] x;

}




/*
//余因子行列を求める関数　※今回は使わない
void calc_cofacter_matrix(double M[3][3],double cofacter_M[3][3])
{
	double cofacter[3][3];
	double transposed_M[3][3];

	for(int D=0;D<DIMENSION;D++)
	{
		for(int D2=0;D2<DIMENSION;D2++)
		{
			if((D!=D)&(D2!=D2))
			{
				if((D+D2)%2==0)	cofacter[D][D2]+=M[(D+1)%3][(D+1)%3]*M[(D+2)%3][(D+2)%3]-M[(D+1)%3][(D+2)%3]*M[(D+2)%3][(D+1)%3];
				else cofacter[D][D2]-=M[(D+1)%3][(D+1)%3]*M[(D+2)%3][(D+2)%3]-M[(D+1)%3][(D+2)%3]*M[(D+2)%3][(D+1)%3];
			}
		}
	}

	calc_transposed_inverse_matrix(cofacter,cofacter_M,0,0,true,false);
}
*/