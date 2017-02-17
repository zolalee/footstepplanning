#include "stdio.h"
#include "stdlib.h"
#include "time.h"
//#include "string.h"
#include "float.h"
#include "math.h"
//#include "process.h" 
#include <string.h>
#include <queue>
#include <stack>
#include <iostream>
#include <fstream>

using namespace std;
FILE *stream_finalsteps;
//FILE *stream_pysteps;
FILE *stream_lowest; 
FILE *stream_allinfo; 
FILE *stream_fp; 
FILE *stream_map; 
FILE *stream_deapest;
FILE *stream_allfp;
FILE *stream_all;

///////////////////////////////////////// 
//const double START[3]={32.72,18.63,-120};
//const double GOAL[3]={176.9723,165.5961,-120};

double START[3]={0.0,0.0,0.0};
double GOAL[3]={0.0,0.0,0.0};

const double Tolerence=0.1;  
const double Max_h_stepupon=5;
const double Max_h_stepover=Max_h_stepupon;
const double COST_OFFSET=0;
const int STEP_FPS_NO=0;         // 1 to delete rest child ft  0- not i=FPS_NO  or i++
double k_d=1;    //distance coefficient
double k_a=0.01;  //0.1 when no obstacles  angle coefficient    
double k_s=0.01;  //self coefficient 
double k_dp=0.01;
double k_h=0.0;
//double k_dp=-0.005;  //  0  0.01 0.015   0.04    0.01~0.04  0.014  plan duration increase dramatically

const double RES=0.01;   // resolution
int Node_num_searched;
const int Max_Node=500;      //*1000
const int Max_Depth=10000;
const double Max_duration=30;

const int GoalReachLevel=1;

#define PI 3.1416
#define FPS_NO 10
#define PIXEL 200     // PIXEL*PIXEL-grid map 
const double Safecoefficent=1;   // resolution


/////////////////////////////////////////
bool GoalReached;
float fps[FPS_NO][7];  //information of potential footprints 
const double Pi_Trans=PI/180;

long stepnum=0;
long stepnum1=0;
long stepnum2=0;
long stepnum3=0;
int kk=0;

const double L1=0.16;
const double L2=0.09;
double L3;
double alfa;

double duration_start;
double duration_finish;
double duration;

enum Supporting{LeftFoot,RightFoot}; //判断支撑腿

struct CCOST
{ double cost_angle;  //评估值  g(x);
double cost_distance;  // assessment value  h(3);
double cost_self;
double cost_sum;
double cost_depth;
};

struct FFPINFO
{ double ave_height;
  double max_height;
  double min_height;
  bool isflat;

};

struct FP//足迹
{   int no; //编号
int fpdepth;
double x;//位置 x
double y;//位置 y
double theta;//位置 theta
struct CCOST COST; //总值
bool leftfoot; // leftfoot=true 为左脚 =false为右脚
struct FFPINFO FPINFO;

struct FP * Parent;//父节点
struct FP * Child[FPS_NO]; //子节点
struct FP * Next; // 链表
};

struct Grid
{	int xx;
	int yy;
	double height;   
	int info;
};

Grid ggrid[PIXEL][PIXEL];  //grid map

struct Point
{	int x;
	int y;
};


struct Quadrangle
{	struct Point P[5];
};  

FP *head,*allfplist;

FP DeapestFP, LowestFP;
stack<FP *>Stack_allfp;
int Deapest;
double Lowestcost;

//double __min(double dataA, double dataB)
//{
//	return ((dataA > dataB) ? dataB : dataA);
//}
//	
//	//template <typename type>
//double __max(double dataA, double dataB)
//{
//	return ((dataA > dataB) ? dataA : dataB);
//}

///////////////  FP initialization  ///////////////////////////////////
struct FP * InitFP(struct FP *NewFP)
{     
	NewFP->no=0;
	NewFP->fpdepth=0;
	NewFP->x=0;
	NewFP->y=0;
	NewFP->theta=0;
	NewFP->COST.cost_angle=0;
	NewFP->COST.cost_self=0;
	NewFP->COST.cost_distance=0;
	NewFP->COST.cost_sum=0;	  
	NewFP->leftfoot=true;

	NewFP->FPINFO.ave_height=0;
	NewFP->FPINFO.isflat=false;
	NewFP->FPINFO.max_height=0;
	NewFP->FPINFO.min_height=100;
	
	NewFP->Parent=NULL;
	
	for(int i=0;i<FPS_NO;i++)       
	{ NewFP->Child[i]=NULL; }
	
	NewFP->Next=NULL;
	return NewFP;
}
//////////////////// Grid Map initialization  //////////////////////////////////
void gridmapinit()
{ 
	int temp;
	for(int m=1; m<=PIXEL;m++)
	{for(int n=1; n<=PIXEL; n++)
    { 
		fscanf(stream_map,"%d ",&temp);
		ggrid[m-1][n-1].height=(double)temp;
		ggrid[m-1][n-1].xx=n; 
		ggrid[m-1][n-1].yy=m;
		ggrid[m-1][n-1].info=0;
	}
	}
}


////////////////// print the Footprint Information     ///////////////////////////////////
void PrintFP(struct FP *TheFP)
{
	printf("%d,%3f,%3f,%3f,%3f,%3f,%3f,%3f \n",TheFP->no,TheFP->x,TheFP->y,TheFP->theta/Pi_Trans,TheFP->COST.cost_angle,TheFP->COST.cost_distance,TheFP->COST.cost_self,TheFP->COST.cost_sum);
	//fprintf(stream_pysteps, "%d \n", TheFP->no);
	fprintf(stream_finalsteps,"%d %3f %3f %3f %3f\n",TheFP->no,TheFP->x/RES,TheFP->y/RES,TheFP->theta/Pi_Trans,TheFP->FPINFO.ave_height);
	fprintf(stream_allinfo,"%d %3f %3f %3f %3f %3f %3f %3f %3f\n",TheFP->no,TheFP->x,TheFP->y,TheFP->theta/Pi_Trans,TheFP->FPINFO.ave_height,TheFP->COST.cost_angle,TheFP->COST.cost_distance,TheFP->COST.cost_self,TheFP->COST.cost_sum);
}

void PrintlowestFP(struct FP *TheFP)
{
	//fprintf(stream_lowest,"%d %3f %3f %3f %3f %3f %3f %3f %3f \n",TheFP->no,TheFP->x,TheFP->y,TheFP->theta/Pi_Trans,TheFP->FPINFO.ave_height,TheFP->COST.cost_angle,TheFP->COST.cost_distance,TheFP->COST.cost_self,TheFP->COST.cost_sum);
}
void PrintdeapestFP(struct FP *TheFP)
{
	fprintf(stream_deapest,"%d %3f %3f %3f %3f %3f %3f %3f %3f \n",TheFP->no,TheFP->x,TheFP->y,TheFP->theta/Pi_Trans,TheFP->FPINFO.ave_height,TheFP->COST.cost_angle,TheFP->COST.cost_distance,TheFP->COST.cost_self,TheFP->COST.cost_sum);
}

void PrintallFP(struct FP *TheFP)
{
	fprintf(stream_allfp,"%d %3f %3f %3f %3f %3f %3f %3f %3f \n",TheFP->no,TheFP->x/RES,TheFP->y/RES,TheFP->theta/Pi_Trans,TheFP->FPINFO.ave_height,TheFP->COST.cost_angle,TheFP->COST.cost_distance,TheFP->COST.cost_self,TheFP->COST.cost_sum);
}


//移动足迹
struct FP * SStep(struct FP * CurFP,int FPNR, Supporting Support)  //此函数有很大问题
{   Node_num_searched++;
	FP *NewFP;
	NewFP=new FP();
	NewFP=InitFP(NewFP);
	
	switch(Support)  //移动足迹
	{
	case LeftFoot:
		{  
			NewFP->no=100+FPNR;
			NewFP->fpdepth=CurFP->fpdepth+1;
			NewFP->x=CurFP->x+fps[FPNR][3]*cos(CurFP->theta+fps[FPNR][4]);
			NewFP->y=CurFP->y+fps[FPNR][3]*sin(CurFP->theta+fps[FPNR][4]);  
			NewFP->theta=CurFP->theta+fps[FPNR][2]*Pi_Trans;
			NewFP->COST.cost_self=fps[FPNR][6];
			NewFP->Parent=CurFP;
			NewFP->leftfoot=false;
			CurFP->Child[0]=NewFP; 	    
			break; 
		}
	case RightFoot:
		{
			NewFP->no=200+FPNR;
			NewFP->fpdepth=CurFP->fpdepth+1;
			NewFP->x=CurFP->x+fps[FPNR][3]*cos(PI+CurFP->theta-fps[FPNR][4]);
			NewFP->y=CurFP->y+fps[FPNR][3]*sin(PI+CurFP->theta-fps[FPNR][4]);  
			NewFP->theta=CurFP->theta-fps[FPNR][2]*Pi_Trans;
			NewFP->COST.cost_self=fps[FPNR][6];
			NewFP->Parent=CurFP;
			NewFP->leftfoot=true;
			CurFP->Child[0]=NewFP;
			break;
		}
	default:
		{	   printf(" default entered...\n");
		break;
		}
	}  
        // copy the scripts from the log file if you want to adjust self-cost value
	
	return NewFP;
}
////////////// get the grid info of a footprint///////////////////////////////////  
struct FFPINFO Get_FPinfo(struct FP * TheFP)
{
	Quadrangle Quad;
	int xmin, xmax, ymin, ymax;
	double P1x, P1y, P2x, P2y, P3x, P3y, P4x, P4y;
	double heighttotal=0; 
	double L3_safe=L3*Safecoefficent;  //膨化 脚掌
	int totalingrid=0;
	double max_height=0;
	double min_height=1000;

	double fline[4];
	double flinep[4][2];
	double flinetime[4][2];
	
	P1x=TheFP->x+L3_safe*cos(TheFP->theta+alfa+PI/2);
	P1y=TheFP->y+L3_safe*sin(TheFP->theta+alfa+PI/2);
	
	P2x=TheFP->x+L3_safe*cos(TheFP->theta-alfa+3*PI/2);
	P2y=TheFP->y+L3_safe*sin(TheFP->theta-alfa+3*PI/2);
	
	P3x=TheFP->x+L3_safe*cos(TheFP->theta+alfa+3*PI/2);
	P3y=TheFP->y+L3_safe*sin(TheFP->theta+alfa+3*PI/2);
	
	P4x=TheFP->x+L3_safe*cos(TheFP->theta-alfa+PI/2);
	P4y=TheFP->y+L3_safe*sin(TheFP->theta-alfa+PI/2);

	if((P1x<=0||P1x>=5||P1y<=0||P1y>=5)||
	   (P2x<=0||P2x>=5||P2y<=0||P2y>=5)||
	   (P3x<=0||P3x>=5||P3y<=0||P3y>=5)||
	   (P4x<=0||P4x>=5||P4y<=0||P4y>=5))
	   {
	   TheFP->FPINFO.ave_height=1000;           // avoid the FP exceed the map area.
	   TheFP->FPINFO.isflat=false; 
        
	   return TheFP->FPINFO;
		}

	Quad.P[0].x=P1x/RES;
	Quad.P[0].y=P1y/RES;
	Quad.P[1].x=P2x/RES;
	Quad.P[1].y=P2y/RES;
	Quad.P[2].x=P3x/RES;
	Quad.P[2].y=P3y/RES;
	Quad.P[3].x=P4x/RES;
	Quad.P[3].y=P4y/RES;
	Quad.P[4].x=Quad.P[0].x;
	Quad.P[4].y=Quad.P[0].y;

	xmin=__min(__min(Quad.P[0].x,Quad.P[1].x),__min(Quad.P[2].x,Quad.P[3].x));
	xmax=__max(__max(Quad.P[0].x,Quad.P[1].x),__max(Quad.P[2].x,Quad.P[3].x));
	ymin=__min(__min(Quad.P[0].y,Quad.P[1].y),__min(Quad.P[2].y,Quad.P[3].y));
	ymax=__max(__max(Quad.P[0].y,Quad.P[1].y),__max(Quad.P[2].y,Quad.P[3].y));
	
	flinep[0][0]=(Quad.P[1].x-Quad.P[0].x)*Quad.P[2].y-(Quad.P[1].y-Quad.P[0].y)*Quad.P[2].x+Quad.P[0].x*Quad.P[1].y-Quad.P[0].y*Quad.P[1].x;
	flinep[0][1]=(Quad.P[1].x-Quad.P[0].x)*Quad.P[3].y-(Quad.P[1].y-Quad.P[0].y)*Quad.P[3].x+Quad.P[0].x*Quad.P[1].y-Quad.P[0].y*Quad.P[1].x;
	flinep[1][0]=(Quad.P[2].x-Quad.P[1].x)*Quad.P[3].y-(Quad.P[2].y-Quad.P[1].y)*Quad.P[3].x+Quad.P[1].x*Quad.P[2].y-Quad.P[1].y*Quad.P[2].x;
	flinep[1][1]=(Quad.P[2].x-Quad.P[1].x)*Quad.P[0].y-(Quad.P[2].y-Quad.P[1].y)*Quad.P[0].x+Quad.P[1].x*Quad.P[2].y-Quad.P[1].y*Quad.P[2].x;
	flinep[2][0]=(Quad.P[3].x-Quad.P[2].x)*Quad.P[0].y-(Quad.P[3].y-Quad.P[2].y)*Quad.P[0].x+Quad.P[2].x*Quad.P[3].y-Quad.P[2].y*Quad.P[3].x;
	flinep[2][1]=(Quad.P[3].x-Quad.P[2].x)*Quad.P[1].y-(Quad.P[3].y-Quad.P[2].y)*Quad.P[1].x+Quad.P[2].x*Quad.P[3].y-Quad.P[2].y*Quad.P[3].x;
	flinep[3][0]=(Quad.P[0].x-Quad.P[3].x)*Quad.P[1].y-(Quad.P[0].y-Quad.P[3].y)*Quad.P[1].x+Quad.P[3].x*Quad.P[0].y-Quad.P[3].y*Quad.P[0].x;
	flinep[3][1]=(Quad.P[0].x-Quad.P[3].x)*Quad.P[2].y-(Quad.P[0].y-Quad.P[3].y)*Quad.P[2].x+Quad.P[3].x*Quad.P[0].y-Quad.P[3].y*Quad.P[0].x;
	
	for (int xx=xmin;xx<=xmax; xx++)
		for (int yy=ymin;yy<=ymax; yy++)
		{
			for (int k=0;k<=3;k++)
			{	
				fline[k]=(Quad.P[k+1].x-Quad.P[k].x)*yy-(Quad.P[k+1].y-Quad.P[k].y)*xx+Quad.P[k].x*Quad.P[k+1].y-Quad.P[k].y*Quad.P[k+1].x;
			}
			flinetime[0][0]=fline[0]*flinep[0][0];
			flinetime[0][1]=fline[0]*flinep[0][1];
			
			flinetime[1][0]=fline[1]*flinep[1][0];
			flinetime[1][1]=fline[1]*flinep[1][1];
			
			flinetime[2][0]=fline[2]*flinep[2][0];
			flinetime[2][1]=fline[2]*flinep[2][1];
			
			flinetime[3][0]=fline[3]*flinep[3][0];
			flinetime[3][1]=fline[3]*flinep[3][1];
			
			if(((flinetime[0][0]>=0)||(flinetime[0][1]>=0))&&
				((flinetime[1][0]>=0)||(flinetime[1][1]>=0))&&
				((flinetime[2][0]>=0)||(flinetime[2][1]>=0))&&
				((flinetime[3][0]>=0)||(flinetime[3][1]>=0)))
			{
				heighttotal=heighttotal+ggrid[yy][xx].height;
				totalingrid++;

				if(ggrid[yy][xx].height>max_height)
					max_height=ggrid[yy][xx].height;
				else 
					if(ggrid[yy][xx].height<min_height)
					min_height=ggrid[yy][xx].height;
					else
					{
					};
			}
			
		}

	   TheFP->FPINFO.ave_height=heighttotal/totalingrid; 
	   TheFP->FPINFO.max_height=max_height; 
	   TheFP->FPINFO.min_height=min_height;
	   
	   if((max_height-min_height)>0.01)
	         TheFP->FPINFO.isflat=false; 
	   else
	         TheFP->FPINFO.isflat=true; 
		
        return TheFP->FPINFO;
}

bool FP_is_feasibile(struct FP * TheFP)
{  	
	double min_h;
	FP *tempFP;
	tempFP=new FP();       // important 
	tempFP=InitFP(tempFP);

	if(TheFP->fpdepth<=1)  //depth start from 0,make sure TheFP->Parent->Parent exists.
		return true;
	else	if(pow(TheFP->x-TheFP->Parent->Parent->x,2)+pow(TheFP->y-TheFP->Parent->Parent->y,2)<(0.005*0.005))
				return false;    // to avoid the FP back to its grapa;
		
			else if(TheFP->FPINFO.isflat!=true)
						  return false;
				 else if (fabs(TheFP->FPINFO.ave_height-TheFP->Parent->FPINFO.ave_height)>(Max_h_stepupon+0.01))
							  return false;
					 else  
					 {   	
						 tempFP->x=(TheFP->x+TheFP->Parent->Parent->x)/2;
						 tempFP->y=(TheFP->y+TheFP->Parent->Parent->y)/2;
						 tempFP->theta=(TheFP->theta+TheFP->Parent->Parent->theta)/2;   // to creat a inter FP;
						 tempFP->FPINFO=Get_FPinfo(tempFP);				 
						 min_h=__min(TheFP->FPINFO.ave_height,TheFP->Parent->Parent->FPINFO.ave_height);
						 
						 if(fabs(tempFP->FPINFO.max_height-min_h)>(Max_h_stepover+0.01))    //to avoid collision when swing the leg   
								 return true;
						else  
							 return true;
					 }
}

//////////////////  Node cost assessing function    ///////////////////////////
struct CCOST  cost_compute(struct FP * TheFP,struct FP* StartFP,struct FP * GoalFP)
{   
	double distance_to_goal; 
	double angle_to_goal;
	double cost_a;
	double cost_d;
	double cost_s;
	double cost_dp;
	double cost_h;   //height information

	
	//distance cost to the goal 
	distance_to_goal=sqrt(pow(TheFP->x-GoalFP->x,2)+pow(TheFP->y-GoalFP->y,2));
	TheFP->COST.cost_distance=distance_to_goal;
	cost_d=distance_to_goal/sqrt(pow(GoalFP->x-StartFP->x,2)+pow(GoalFP->y-StartFP->y,2));
	
	//angle cost to the goal
	angle_to_goal=atan2(GoalFP->y-TheFP->y,GoalFP->x-TheFP->x);
	TheFP->COST.cost_angle=fabs(angle_to_goal-(TheFP->theta+PI/2));
	cost_a=TheFP->COST.cost_angle/PI;
	
    cost_s=TheFP->COST.cost_self;
	TheFP->COST.cost_depth=TheFP->fpdepth;
	cost_dp=TheFP->COST.cost_depth;
    
	if(TheFP->Parent==NULL)     //avoid start FP error;
		cost_h=0;
	else 
		cost_h=fabs((TheFP->FPINFO.ave_height-TheFP->Parent->FPINFO.ave_height))/4;
	
    TheFP->COST.cost_sum=k_d*cost_d
						+k_a*cost_a
						+k_s*cost_s
						+k_dp*cost_dp
						+k_h*cost_h;
	
	return TheFP->COST;
}

bool Goalisreached(struct FP * TheFP)
{   
	bool temp;
	switch(GoalReachLevel)  
	{
	case 4: 
		{
			temp=((TheFP->COST.cost_distance<(Tolerence+0.00001))
				&&(TheFP->Parent->COST.cost_distance<(Tolerence+0.00001))
				&&((TheFP->no==202)||(TheFP->no==102)));
			break;
		}
	case 3: 
		{
			temp=((TheFP->COST.cost_distance<(Tolerence+0.00001))
				&&(TheFP->Parent->COST.cost_distance<(Tolerence+0.00001)));
			break;
		}

	case 2:
		{				
			temp=(TheFP->COST.cost_distance<(Tolerence+0.00001))
				 &&((TheFP->no==201)||(TheFP->no==101)); 
			break;
		}
	case 1:
		{	
			temp=(TheFP->COST.cost_distance<(Tolerence+0.00001));	
				break;
		}
	default:
		{   printf(" goalreaching default entered...\n");
		break;
		}
	}
		
		if(temp)
			GoalReached=true;
		else 
			GoalReached=false;
	
		return GoalReached;
}
//////////////////// core function SEARCH  /////////////////////
struct FP * Search(struct FP* Start,struct FP * Goal)   
{
    FP *fp1,*fp2,*fp;
    fp=NULL;
	Supporting Support=(Supporting)0;
	
    do{                                     //search
        fp1=head;
		head=head->Next;
		
		      // if(Collision(fp2,colcheck)!=false)
		      //   fp1=head->Next;
		//   lazy-evaluation approach  . see lazy
		
		if (fp1->leftfoot==true)
			Support=(Supporting)0;
		else 
			Support=(Supporting)1;
		
		//  if fp1->cost_angle < ....   //to realized alternative potential footprints
		//	for (int j=0;j<=20;j++)
		// {	Move2 Mov=(Move2)j;
		
        for(int FPNR=0;FPNR<=(FPS_NO-1);FPNR++)//FPS_NO个可行足迹获得新节点
        {
			fp2=SStep(fp1,FPNR,Support);
			fp2->FPINFO=Get_FPinfo(fp2);           
			fp2->COST=cost_compute(fp2,Start,Goal);			 
			
			 		
			if(FP_is_feasibile(fp2)!=true)
			{  
			    if(STEP_FPS_NO==1)
				FPNR=FPS_NO;   //delete FP if any child of FP collides.
				else
				{
				};
			}
						
			else 
			{   
				Stack_allfp.push(fp2);
					
             if(fp2->fpdepth>Deapest)
				{	Deapest=fp2->fpdepth;
					DeapestFP=*fp2;
				}

            if(fp2->COST.cost_sum<Lowestcost)
				{	Lowestcost=fp2->COST.cost_sum;
					LowestFP=*fp2;
				}

				if(Goalisreached(fp2)==true)	
				{  fp=fp2;
					FPNR=FPS_NO;
				}

				if(head==NULL||head->COST.cost_sum>=(fp2->COST.cost_sum-COST_OFFSET))
				{  FP *temp=head;
				head=fp2;
				fp2->Next=temp; }

				else
				{ 	FP *temp=head;  
					FP *oldtemp=temp;
				
					while(temp!=NULL&&temp->COST.cost_sum<(fp2->COST.cost_sum-COST_OFFSET))
					{  oldtemp=temp;
					temp=temp->Next;
					}
				
					if(temp==NULL||oldtemp->Next==NULL)
					{	oldtemp->Next=fp2;
						fp2->Next=NULL;		
					}
					else
					{	oldtemp->Next=fp2;
						fp2->Next=temp; 	
					}
				}    		

			} 									
		} 
		
        duration_finish=double(clock());
		duration=(duration_finish-duration_start)/CLOCKS_PER_SEC;

		if(duration>=Max_duration)
			return head;

	//	if(Node_num_searched>=Max_Node*1000-1)
      //     return head;   // not the lowest one //////////////////
		
	}while(fp==NULL&&head->Next!=NULL);
	
    return fp;  
} 

void result_print1()
{
	printf("\n\n\n*********************************************\n");
	printf("Candiate FootPrints are are:\n");
	printf("*********************************************\n\n");
	printf("|   x   |   y   | theta | distance | alfa1  |  alfa2  |\n");
	
	fprintf(stream_allinfo,"\n\n\n*********************************************\n");
	fprintf(stream_allinfo,"Candiate FootPrints are are:\n");
	fprintf(stream_allinfo,"*********************************************\n\n");
	fprintf(stream_allinfo,"|   x   |   y   | theta |  distance | alfa1  |  alfa2  |\n");
}

void result_print2()
{       
	printf("\n\n*********************************************\n");
	printf("Predefined parameters are:\n");
	printf("*********************************************\n");
	
	printf("FP Number:                  %d \n",FPS_NO);
	printf("Tolerence:                  %f \n",Tolerence);
	printf("STEP_FPS_NO:                %d /1-yes/0-no\n",STEP_FPS_NO);
	printf("Safecoefficent:             %f \n",Safecoefficent);
	printf("GridMap:                    %d*%d. Resolution=%f \n",PIXEL,PIXEL,RES);
	printf("Start Postion:              %f, %f \n",START[0]*RES,(500-START[1])*RES);
	printf("Goal Postion:               %f, %f \n",GOAL[0]*RES,(500-GOAL[1])*RES);
	printf("COST_OFFSET:                %f \n",COST_OFFSET);
	printf("Max_h_stepupon:             %f \n",Max_h_stepupon);
	printf("Max_h_stepover:             %f \n",Max_h_stepover);
	printf("Cost(k_d/k_a/k_s/k_dp):     %f   %f  %f  %f \n",k_d,k_a,k_s,k_dp);
	printf("GoalReachLevel:             %d (1-easy restriction)\n",GoalReachLevel);

	
	fprintf(stream_allinfo,"\n\n*********************************************\n");
	fprintf(stream_allinfo,"Predefined parameters are:\n");
	fprintf(stream_allinfo,"*********************************************\n\n");
	fprintf(stream_allinfo,"FP Number:                        %d \n",FPS_NO);
	fprintf(stream_allinfo,"Tolerence:                        %f \n",Tolerence);
	fprintf(stream_allinfo,"STEP_FPS_NO(1-y/0-n):             %d \n",STEP_FPS_NO);
	fprintf(stream_allinfo,"Safecoefficent:                   %f \n",Safecoefficent);
	fprintf(stream_allinfo,"GridMap:                          %d*%d. Resolution=%f \n",PIXEL,PIXEL,RES);
	fprintf(stream_allinfo,"Start Postion:                    %f, %f \n",START[0]*RES,(500-START[1])*RES);
	fprintf(stream_allinfo,"Goal Postion:                     %f, %f \n",GOAL[0]*RES,(500-GOAL[1])*RES);
	fprintf(stream_allinfo,"COST_OFFSET:                      %f \n",COST_OFFSET);
	fprintf(stream_allinfo,"Max_h_stepupon:                   %f \n",Max_h_stepupon);
	fprintf(stream_allinfo,"Max_h_stepover:                   %f \n",Max_h_stepover);
	fprintf(stream_allinfo,"Cost(k_d/k_a/k_s/k_dp):           %f   %f  %f  %f \n",k_d,k_a,k_s,k_dp);
	fprintf(stream_allinfo,"GoalReachLevel:                   %d  (1-easy restriction)\n",GoalReachLevel);

	printf("\n*********************************************\n");
	printf("search finisned:\n");
	printf("*********************************************\n\n");
	printf("Total Steps are:            %d \n\n",stepnum);
	printf("Max depth:                  %d \n",Deapest);
	printf("Total Nodes searched are:   %d \n",Node_num_searched);
	printf("Total duration is :         %2.8f seconds\n", duration);
		
    fprintf(stream_allinfo,"\n*********************************************\n");
	fprintf(stream_allinfo,"search finisned, results are: \n");
	fprintf(stream_allinfo,"*********************************************\n\n");
	fprintf(stream_allinfo,"Total Steps are:                  %d \n",stepnum);
	fprintf(stream_allinfo,"Max depth:                        %d \n",Deapest);
	fprintf(stream_allinfo,"Total Nodes searched are:         %d \n",Node_num_searched);
	fprintf(stream_allinfo, "total duration :                 %2.8f seconds\n", duration);
	fprintf(stream_allinfo,"Nodes in List is:                 %d \n",stepnum1);
	fprintf(stream_allinfo,"Lowest-cost Node Depth is:        %d \n",stepnum2);
	fprintf(stream_allinfo,"The largest Depth is:             %d \n",Deapest);
}

void result_print3()
{
	printf("\n\n\n*********************************************\n");
	printf("Search Results\n\n");
	printf("*********************************************\n");
	printf("|N0.|  x |  y | theta | cost_angle | cost_distance |cost_self| cost_sum |\n");
	fprintf(stream_allinfo,"\n\n\n*********************************************\n");
	fprintf(stream_allinfo,"Search Results\n");
	fprintf(stream_allinfo,"*********************************************\n");
	fprintf(stream_allinfo,"|N0.|  x  |  y  | theta | cost_angle |  cost_distance  |  cost_self  |  cost_sum |\n");
}
long filesize(FILE *stream) 
{ 
   long curpos, length; 

   curpos = ftell(stream); 
   fseek(stream, 0L, SEEK_END); 
   length = ftell(stream); 
   fseek(stream, curpos, SEEK_SET); 
   return length; 
}
/*
void updatefootsteptxt(int num_txt)
{
int buffsize;
int begin;
int end;
char* buffer;
int i=num_txt;
char str[6];
itoa(i,str,10);
char *s2=".txt";
strcat(str,s2);
FILE* inFile = fopen("finalsteps.txt","rb");
FILE* outFile = fopen(str,"wb+");
if(inFile!=NULL)
{
buffsize=filesize(inFile);
 printf("%d",buffsize);
fseek(inFile,0,SEEK_SET);
buffer=(char*)malloc(buffsize*sizeof(char));
fread(buffer,buffsize,1,inFile);
fwrite(buffer,buffsize,1,outFile);
cout<<"job done!";
}
else
{
fprintf(stderr,"error");
}
free(buffer);
}
*/


//***********************************************************************//
//***********************************************************************//
//***********************************************************************//
//void FilesBackup(double num_txt)
//{
//	ifstream fin_m("map.txt", ios::app|ios::binary);
//	if(!fin_m)
//	{
//		cout<<"File open error!"<<endl;
//		//exit(1);
//	}
//	int i=num_txt;
//	char str_m[2];
//	char map1[7]="map_";
//	sprintf(str_m, "%d", i); 
//	strcat(map1,str_m);
//	char map2[5]=".txt";
//	strcat(map1,map2);
//
//	ofstream fout_m(map1,ios::binary);
//	char c_m[1024];
//	while(!fin_m.eof())
//	{
//		fin_m.read(c_m,1024);
//		fout_m.write(c_m,fin_m.gcount());
//	}
//	fin_m.close();
//	fout_m.close();
///////////////////////////////////////////////////////////////////////
//	ifstream fin_f("finalsteps.txt", ios::app|ios::binary);
//	if(!fin_f)
//	{
//		cout<<"File open error!"<<endl;
//		//exit(1);
//	}
//	char str_f[2];
//	char finalsteps_f1[12]="finalsteps_";
//	sprintf(str_f, "%d", i); 
//	strcat(finalsteps_f1,str_f);
//	char finalsteps_f2[5]=".txt";
//	strcat(finalsteps_f1,finalsteps_f2);
//
//	ofstream fout_f(finalsteps_f1,ios::binary);
//	char c_f[1024];
//	while(!fin_f.eof())
//	{
//		fin_f.read(c_f,1024);
//		fout_f.write(c_f,fin_f.gcount());
//	}
//	fin_f.close();
//	fout_f.close();
//
//	cout<<"job done!"<<endl;
//}

//***********************************************************************//
//***********************************************************************//



/////////////////////// main  //////////////////////////////////
void main(double S0, double S1, double S2, double G0, double G1, double G2, double PlanNum)
{
    START[0]=S0;
	START[1]=S1;
	START[2]=S2;
	GOAL[0]=G0;
	GOAL[1]=G1;
	GOAL[2]=G2;
	double planNum=PlanNum;
    stream_finalsteps=fopen( "finalsteps.txt", "w" );
	//stream_pysteps=fopen("pysteps.txt","w");
	stream_lowest=fopen( "finalsteps_lowest.txt", "w" );  //for low cost steps
	stream_allinfo=fopen( "results_all_info.txt", "w" );
	stream_allfp=fopen( "finalsteps_allfp.txt", "w" );
	stream_fp=fopen( "FP_10.txt", "r" );
	stream_map=fopen( "map.txt", "r" );
	stream_deapest=fopen( "finalsteps_deapest.txt", "w" );  //for Deapest cost
	stream_all=fopen( "compondmodel_compare.txt", "a" );

///////////////////////////////////////////////////////////////////////	
//////////////////////////////////////////////////////////////////////	
	result_print1();
	
	for (int fpi=0;fpi<=(FPS_NO-1);fpi++)    // fpi--FPNR
    { fscanf(stream_fp,"%f %f %f %f %f %f %f \n",&fps[fpi][0],&fps[fpi][1],&fps[fpi][2],&fps[fpi][3],&fps[fpi][4],&fps[fpi][5],&fps[fpi][6]);
	printf("%d %f %f %f %f %f %f %f \n",fpi,fps[fpi][0],fps[fpi][1],fps[fpi][2],fps[fpi][3],fps[fpi][4],fps[fpi][5],fps[fpi][6]);
	fprintf(stream_allinfo,"%d %f %f %f %f %f %f %f \n",fpi,fps[fpi][0],fps[fpi][1],fps[fpi][2],fps[fpi][3],fps[fpi][4],fps[fpi][5],fps[fpi][6]);
	}; 
    
	L3=sqrt(L1*L1+L2*L2);
	alfa=atan(L2/L1);
	GoalReached=false;
	
	FP  Start,Goal,* T;

	gridmapinit();
	Deapest=0;
	Lowestcost=1000;
	
	DeapestFP=*InitFP(&DeapestFP);
	Start=*InitFP(&Start);
	Goal=*InitFP(&Goal);
	Goal.no=300;
	
	Goal.theta=START[2]*Pi_Trans;
	Goal.x=GOAL[0]*RES;
	Goal.y=(PIXEL-GOAL[1])*RES; 
	Goal.FPINFO=Get_FPinfo(&Goal);
	
	Start.no=0;
	Start.theta=START[2]*Pi_Trans;
	Start.x=START[0]*RES;
	Start.y=(PIXEL-START[1])*RES;
	Start.COST=cost_compute(&Start,&Start,&Goal);
	Start.COST.cost_sum=100;
	Start.FPINFO=Get_FPinfo(&Start);


    head=&Start;  
    
	duration_start=double(clock());     //to compute the planning duration  
	T=Search(&Start,&Goal);  
	duration_finish=double(clock());
	duration=(duration_finish-duration_start)/CLOCKS_PER_SEC;
	
	if(GoalReached==true)   // print the final prints
    {   	
        FP *fp=T;
        stack<FP *>Stack1;
        while(fp->Parent!=NULL)
        {
            Stack1.push(fp);
            fp=fp->Parent;
        }
		
		result_print3();
		PrintFP(&Start);
	    PrintFP(&Goal);    
		
        while(!Stack1.empty())
        {  
			PrintFP(Stack1.top());
			if(Stack1.top()->FPINFO.isflat)
			printf("1\n");
			Stack1.pop();
			stepnum++;
        }
		  
		result_print2();
	
	}
   
	{         
		FP *fp3;
		FP *fp4;
		fp3=&LowestFP;
		fp4=&DeapestFP;

        stack<FP *>Stack2;
		stack<FP *>Stack3;
		

        while(fp3->Parent!=NULL)
        {
            Stack2.push(fp3);
            fp3=fp3->Parent;
        }

		while(fp4->Parent!=NULL)
        {
            Stack3.push(fp4);
            fp4=fp4->Parent;
        }

		PrintlowestFP(&Start);
		PrintlowestFP(&Goal);  
		PrintdeapestFP(&Start);
		PrintdeapestFP(&Goal);
		PrintallFP(&Start);
		PrintallFP(&Goal);
 
        while(!Stack2.empty())
        {  
			PrintlowestFP(Stack2.top());
			Stack2.pop();
			stepnum1++;
        }
		while(!Stack3.empty())
        {  
			PrintdeapestFP(Stack3.top());
			Stack3.pop();
			stepnum2++;
        }
		while(!Stack_allfp.empty())
        {  
			PrintallFP(Stack_allfp.top());
			Stack_allfp.pop();
			stepnum3++;
		 
        }
 
 
         if(GoalReached==false)
		 {
		  printf("\n\nOooooooooops...\n search failed.....\n\n");
	      printf("Total Nodes searched are:   %d \n",Node_num_searched);
	      printf("Total duration is :         %2.8f seconds\n", duration);
         }

		printf("\nLowest-cost Nodes List is:  %d \n",stepnum1);
		printf("largest Depth is:           %d %d\n",Deapest,stepnum2);
		printf("Max_Node is:                %d * 1000\n",Max_Node);
		printf("All_available_Node is:      %d \n",stepnum3);
	    printf("Max_Depth is:               %d \n\n",Max_Depth);
		printf("GoalReachLevel:             %d (1-easy restriction)\n",GoalReachLevel);
		printf("\nLowest-cost is:             %f \n\n",Lowestcost);

	};
		
		fclose(stream_finalsteps);
		//fclose(stream_pysteps);
		//fclose(stream_lowest);
		fclose(stream_allinfo);  
		fclose(stream_fp);
		fclose(stream_map);
	    fclose(stream_deapest);
		fclose(stream_allfp);

		fprintf(stream_all,"%d %d %d %2.8f %d \n",FPS_NO,Node_num_searched,stepnum3,duration,stepnum);
        fclose(stream_all);
		//FilesBackup(planNum);


		//updatefootsteptxt(5);
		//cin.get();
		//cin.get();

}
//extern "C"
//{
//   void test0(double S0, double S1, double S2, double G0, double G1, double G2,double PlanNum)
//   {
//      main(S0, S1, S2, G0, G1, G2, PlanNum);
//   }
//}
