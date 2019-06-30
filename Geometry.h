
void povorot_Mx_x(double **Matrix, double *&X)
{
	int i, j;
	double newX[3];
	for (i = 0; i < 3; i++){
		newX[i] = 0;
		for (j = 0; j < 3; j++)
			newX[i] = newX[i] + Matrix[i][j] * X[j];
	}

	for (j = 0; j < 3; j++)
		X[j] = newX[j];
	i++;
}




class solid_object 
{
	public:
	double R,
		   r,
		   R1,
		   r1,
		   a,
		   b,
		   h;

	double alpha_angle;

	int N_additional_parameters;

	int solid_or_gas;

	bool empty;

	double vx,
		   vy,
		   vz;

	double rho;
	double U;
	double P;
	double T;

	double *volume_ratio;
	#define NAME_LENGTH 100
	char  type_object[NAME_LENGTH];


	double *additional_parameter; 

	double *a_max;//for max_geom
	double *vnz, *va, *center;
	double *min,*max;
	double **Mz1, **Mz2, **My, **Muz1, **Muz2, **Muy;

	double *point;



	void unpovorot_vector( )
	{
		


		if (Muz2[2][2] != 0)
		{
			povorot_Mx_x(Muz2, point);
		}

		if (Muy[1][1]  != 0 )
		{
			povorot_Mx_x(Muy, point);
		
		}

		if ( Muz1[2][2] != 0)
		{
			povorot_Mx_x(Muz1, point);
		}




	}



	void povorot_vector( double x,double y,double z)
	{
			
		point[0] = x - center[0];
		point[1] = y - center[1];
		point[2] = z - center[2];
		


		if ( Mz1[2][2] != 0)
		{
			povorot_Mx_x(Mz1, point);
		}


		if (My[1][1] != 0 )
		{
			povorot_Mx_x(My, point);
		
		}




			if (Mz2[2][2] != 0)
		{
			povorot_Mx_x(Mz2, point);
		}


	}








	void define()
		{
		 R=0;
		 r=0;
		 R1=0;
		 r1=0;
		 a=0;
		 b=0;
		 h=0;
		 
		alpha_angle=0;



		 N_additional_parameters=0;
		 solid_or_gas=0;
		 

		 //empty=false;

		 vx=0;
		 vy=0;
		 vz=0;

		 rho=0;
		 U=0;
		 P=0;
		 T=0;

		AllocateArray(volume_ratio,NS);
		AllocateArray(va,3);
		AllocateArray(vnz,3);

		va[0]=1;
		va[1]=0;
		va[2]=0;

		vnz[0]=0;
		vnz[1]=0;
		vnz[2]=1;

		if(N_additional_parameters>0)
		AllocateArray(additional_parameter,N_additional_parameters);

		AllocateArray(a_max,3);
		AllocateArray(center,3);
		AllocateArray(min,3);
		AllocateArray(max,3);
		AllocateArray(point,3);

		AllocateArray2(Mz1,3,3);
		AllocateArray2(Mz2,3,3);
		AllocateArray2(My,3,3);
		AllocateArray2(Muz1,3,3);
		AllocateArray2(Muz2,3,3);
		AllocateArray2(Muy,3,3);

	}


	void delete_geom()
	{
		if(N_additional_parameters>0)
			DeleteArray(additional_parameter);

		DeleteArray(volume_ratio);

		DeleteArray(a_max);
		DeleteArray(vnz);
		DeleteArray(center);
		DeleteArray(min);
		DeleteArray(max);
		DeleteArray(va);

		DeleteArray2(Mz1,3);
		DeleteArray2(Mz2,3);
		DeleteArray2(My,3);
		DeleteArray2(Muz1,3);
		DeleteArray2(Muz2,3);
		DeleteArray2(Muy,3);
		DeleteArray(point);


	}


	int get_min_max()
		{






		int n=8,i,j;

		double **chek_points;
		AllocateArray2(chek_points,n,3);


		a_max[0]=2*R;
		a_max[1]=2*R;
		a_max[2]=h;



		if(strcmp("BOX", type_object)==0)
		{
			a_max[0]=a;
			a_max[1]=b;
			a_max[2]=h;
		}

		if(strcmp("SPHERE", type_object)==0)
		{
			a_max[0]=2*R;
			a_max[1]=2*R;
			a_max[2]=2*R;
		}




			

		chek_points[0][0] = -a_max[0] / 2;
		chek_points[0][1] = -a_max[1] / 2;
		chek_points[0][2] = 0;

		chek_points[1][0] = +a_max[0] / 2;
		chek_points[1][1] = -a_max[1] / 2;
		chek_points[1][2] = 0;

		chek_points[2][0] = -a_max[0] / 2;
		chek_points[2][1] = +a_max[1] / 2;
		chek_points[2][2] = 0;

		chek_points[3][0] = a_max[0] / 2;
		chek_points[3][1] = a_max[1] / 2;
		chek_points[3][2] = 0;

		chek_points[4][0] = -a_max[0] / 2;
		chek_points[4][1] = -a_max[1] / 2;
		chek_points[4][2] = a_max[2];

		chek_points[5][0] = +a_max[0] / 2;
		chek_points[5][1] = -a_max[1] / 2;
		chek_points[5][2] = a_max[2];

		chek_points[6][0] = -a_max[0] / 2;
		chek_points[6][1] = +a_max[1] / 2;
		chek_points[6][2] = a_max[2];

		chek_points[7][0] = a_max[0] / 2;
		chek_points[7][1] = a_max[1] / 2;
		chek_points[7][2] = a_max[2];
			////////////////////////////////////////


		for (i = 0; i < 3; i++)
			point[i] = chek_points[0][i];

		unpovorot_vector();


		for (i = 0; i < 3; i++)
		{
			min[i] = point[i];
			max[i] = point[i];
		}




			for (j = 1; j < n; j++)
		{
			for (i = 0; i < 3; i++)
				point[i] = chek_points[j][i];


			unpovorot_vector();


			for (i = 0; i < 3; i++)
			{
				if (point[i] < min[i])
					min[i] = point[i];
				if (point[i] > max[i])
					max[i] = point[i];
			}


			



		}

		for (i = 0; i < 3; i++)
		{
			min[i] = min[i] + center[i];
			max[i] = max[i] + center[i];
		}

		if(strcmp("SPHERE", type_object)==0)
		{
			max[2]=max[2]-R;
			min[2]=min[2]-R;
		}

		DeleteArray2(chek_points,n);
			
		return 1;
	}









	void get_Matrix()
	{

	double coss,sinn;
	int i,j;

	for(i=0;i<3;i++)
	  for(j=0;j<3;j++)
		{
		Mz1[i][j]=0;
		Muz1[i][j]=0;

		My[i][j]=0;
		Muy[i][j]=0;

		Mz2[i][j]=0;
		Muz2[i][j]=0;		
		}



		if ( vnz[1] != 0)
		{
			coss = -vnz[0] / pow((vnz[0] * vnz[0] + vnz[1] * vnz[1]), 0.5);
			sinn = +vnz[1] / pow((vnz[0] * vnz[0] + vnz[1] * vnz[1]), 0.5);

			Mz1[0][0] = coss;	Mz1[0][1] = -sinn;		Mz1[0][2] = 0;
			Mz1[1][0] = sinn;	Mz1[1][1] = coss;		Mz1[1][2] = 0;
			Mz1[2][0] = 0;		Mz1[2][1] = 0;			Mz1[2][2] = 1;



			Muz1[0][0] = coss;	Muz1[1][0] = -sinn;		Muz1[2][0] = 0;
			Muz1[0][1] = sinn;	Muz1[1][1] = coss;		Muz1[2][1] = 0;
			Muz1[0][2] = 0;		Muz1[1][2] = 0;			Muz1[2][2] = 1;

			povorot_Mx_x(Mz1, va);
			povorot_Mx_x(Mz1, vnz);
		}


		if (vnz[0] != 0 )
		{
			coss =  vnz[2] / pow((vnz[2] * vnz[2] + vnz[0] * vnz[0]), 0.5);
			sinn = +vnz[0] / pow((vnz[2] * vnz[2] + vnz[0] * vnz[0]), 0.5);

			My[0][0] = coss;	My[0][1] = 0;	My[0][2] = -sinn;
			My[1][0] = 0;		My[1][1] = 1;	My[1][2] = 0;
			My[2][0] = sinn;	My[2][1] = 0;	My[2][2] = coss;


			Muy[0][0] = coss;	Muy[1][0] = 0;	Muy[2][0] = -sinn;
			Muy[0][1] = 0;		Muy[1][1] = 1;	Muy[2][1] = 0;
			Muy[0][2] = sinn;	Muy[1][2] = 0;	Muy[2][2] = coss;


			povorot_Mx_x(My, va);
			povorot_Mx_x(My, vnz);

		}



			if ( va[1] != 0)
		{
			coss =	-va[0]/ pow((va[0] * va[0] + va[1] * va[1]), 0.5);
			sinn = +va[1]/ pow((va[0] * va[0] + va[1] * va[1]), 0.5);

			Mz2[0][0] = coss;	Mz2[0][1] = -sinn;		Mz2[0][2] = 0;
			Mz2[1][0] = sinn;	Mz2[1][1] = coss;		Mz2[1][2] = 0;
			Mz2[2][0] = 0;		Mz2[2][1] = 0;			Mz2[2][2] = 1;


			Muz2[0][0] = coss;	Muz2[1][0] = -sinn;		Muz2[2][0] = 0;
			Muz2[0][1] = sinn;	Muz2[1][1] = coss;		Muz2[2][1] = 0;
			Muz2[0][2] = 0;		Muz2[1][2] = 0;			Muz2[2][2] = 1;
		
			
		}


	}



	bool in_object()
		{
	if(strcmp("BOX", type_object)==0)
			return Box();
			
			
	if(strcmp("CYLINDER", type_object)==0)
			return Cylinder();

		
	if(strcmp("SPHERE", type_object)==0)
			return Sphere();

		
	if(strcmp("RING", type_object)==0)
			return Ring();


	if(strcmp("SECTOR", type_object)==0)
			return Sector();

	if(strcmp("CONE", type_object)==0)
			return Cone();

			return false;

	}


	bool check_cord_object(double x,double y,double z)
	{

	if(x<min[0] || x>max[0])
				return false;
	if(y<min[1] || y>max[1])
				return false;
	if(z<min[2] || z>max[2])
				return false;

			
	povorot_vector(x,y,z);
		if(in_object())
			{
			return true;
			}

	return false;
	}


	/////////////////////	GEOMETRY START	 ///////////////////

	bool Box()
	{
			if(fabs(point[0])<=a/2 && fabs(point[1])<=b/2 && fabs(point[2]-h/2)<=h/2)
				return	true;
			else 
				return false;
	}

	bool Cylinder()
	{
				if((point[0]*point[0]+point[1]*point[1])<=R*R && point[2]<=h &&  point[2]>=0)
					return	true;
				else 
					return false;
	}

	bool Ring()
	{
				if((point[0]*point[0]+point[1]*point[1])>=r*r && (point[0]*point[0]+point[1]*point[1])<=R*R && point[2]<=h &&  point[2]>=0)
					return	true;
				else 
					return false;
	}

	bool Sphere()
	{
				if ((point[0] * point[0] + point[1] * point[1]+point[2]* point[2]) <= R*R  )
				   return	true;
				else
				   return false;
	}


	bool Sector()
	{
		if((point[0]*point[0]+point[1]*point[1])>=r*r && (point[0]*point[0]+point[1]*point[1])<=R*R && point[2]<=h &&  point[2]>=0)
		{
		if(alpha_angle==180)
		{
			if(point[0]>=0)
				return	true;
			else 
				return false;
		}//if(alpha_angle==180)
		if(alpha_angle < 180)
		{
			if(point[0]>=0 && fabs(point[1]/point[0])<=fabs(tan(alpha_angle/2*Pi/180.)))
				 return	true;
			else 
				return false;
		}//(alpha_angle < 180)

		if(alpha_angle > 180)
		{
			if(point[0]<=0 || fabs(point[0]/point[1])<=fabs(tan((alpha_angle-180)/2*Pi/180.)))
				 return	true;
			else 
				return false;
		}//(alpha_angle > 180)

		  return false;
				
		}
		else 
			return false;
	}


	bool Cone()
	{
	double r_cone, R_cone;


	if(r1==r)
	r_cone=r;

	if(R1==R)
	R_cone=R;

	if(r1<r)
	r_cone=r1+(r-r1)/h*(h-point[2]);

	if(R1<R)
	R_cone=R1+(R-R1)/h*(h-point[2]);

	if(r1>r)
	r_cone=r+(r1-r)/h*(point[2]);

	if(R1>R)
	R_cone=R+(R1-R)/h*(point[2]);



				if((point[0]*point[0]+point[1]*point[1])>=r_cone*r_cone && (point[0]*point[0]+point[1]*point[1])<=R_cone*R_cone && point[2]<=h &&  point[2]>=0)
					return	true;
				else 
					return false;

	}




/////////////////////	GEOMETRY END	 ///////////////////

};









