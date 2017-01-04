// This is an example c++ (really c-style c++) code for analyzing MD results
// The basic components are essentially the same as the perl script lmps.analysis.pl. However, note the differences in implementing and speed of these two languages.

//Tasks
//(1) Setup the neighbors of each cell (Expand function SimulationBox::SetupCell() )
//(2) Throw atoms in Each Cell (Write function SimulationBox::DistributeAtoms() )
//(3) Calculate RDF using LinkCell

//Memory check in CCNI
// valgrind --leak-check=full ./md.cell dump.melt haha 1 3

#include "stdio.h"
#include "stdlib.h"
#include "math.h"
#include "string.h"

#define ONELINEMAX 2000
#define MAXDUMP 10000
#define MaxCellAtoms 1000
#define MaxCellNbr 100
#define PI 3.1415926
#define VERBOSE 1


/////////////////////////////////
// Global variables defination //
/////////////////////////////////

class Atom {  // this is about the information of every single atom
public:
  int id;
  int type;
  float x[3];
  float v[3];
  float f[3];
  void PrintInfo() {
    printf("Atom: id[%d] type[%d] x[%f,%f,%f]\n", id,type,x[0],x[1],x[2]);
  };
};

class Cell {   // this is about the cell formed by dividing the box 
public:
  int member[MaxCellAtoms];
  int nmember;
  int nbr[MaxCellNbr];
  int nnbr;
  float origin[3];
  float size[3];
  Cell() {nmember=0;nnbr=0;};
};

class SimulationBox { 
public:
  Atom **myatom;
  Cell ****mycell;
  int ncell[3]; //dimension of the cells
  float L[3];
  float origin[3]; //origin of the simulation box, default will be (0,0,0)
  int nmols;
  void SetupCell(float cellsize) {    
#if VERBOSE ==1 
    printf("Initiating setupcell\n");
#endif
    int i,j,k,p,q,r,pi,qi,ri,l;
    if(nmols==0) {
      printf("Simulation Box is empty!\n");
      exit(0);
    }
    for(int d=0;d<3;d++) {
      ncell[d] = (int)(L[d]/cellsize)-1; //make it a little larger than necessary
      if(ncell[d]<4) {ncell[d]=1;} //if the cells are too few, there are no advantages of doing link-cell.
    }
#if VERBOSE ==1 
    printf("cell[%d,%d,%d]\n", ncell[0],ncell[1],ncell[2]);
#endif
    if(ncell[0]==1||ncell[1]==1||ncell[2]==1){printf ("Warning!!! Cells are too few to use linkcell!!! Please decrese cellsize.\n"); exit(0);}
    mycell = new Cell ***[ncell[0]];
    for(i=0;i<ncell[0];i++) {
      mycell[i] = new Cell **[ncell[1]];
      for(j=0;j<ncell[1];j++) {
	mycell[i][j] = new Cell *[ncell[2]];
	for(k=0;k<ncell[2];k++) {
	  mycell[i][j][k] = new Cell();
	}
      }
    }

    //Setup NBR cells, Coding here
///////////////////////////////////////////
///////////////////////////////////////////
///////Yongjian code following/////////////
///////////////////////////////////////////
///////////////////////////////////////////
   //loop of cells
  for(i=0;i<ncell[0];i++) {
     for(j=0;j<ncell[1];j++) {
       for(k=0;k<ncell[2];k++) {
//printf("k has a value: %d\n",k);
 
          //loop of neighbour cells 
          mycell[i][j][k]->nnbr=27;l=0; 
        for (pi=i-1;pi<i+2;pi++){ 
           for (qi=j-1;qi<j+2;qi++){
              for (ri=k-1;ri<k+2;ri++){
//printf("ri has a value: %d\n",ri);
         //convert neighbour cells outside the range (periodic boundary condition) 
         p=pi;q=qi;r=ri;l=(pi-i+1)*9+(qi-j+1)*3+(ri-k+1);     
      	   if (l<0 || l>26) {printf( "l:%d is wrong",l); exit(0);}
           if (p==-1){p=ncell[0]-1;}
           if (p==ncell[0]){p=0;}  //if (p==-1||p==ncell[0]){p=abs(abs(p)-ncell[0]);}
           if (q==-1){q=ncell[1]-1;}
           if (q==ncell[1]){q=0;}  //if (q==-1||q==ncell[1]){q=abs(abs(q)-ncell[1]);} 
           if (r==-1){r=ncell[2]-1;} 
           if (r==ncell[2]){r=0;}  //if (r==-1||r==ncell[2]){r=abs(abs(r)-ncell[2]);}
          //convert every neighbour's vector index into a scalor index 
            mycell[i][j][k]->nbr[l]=p*ncell[1]*ncell[2]+q*ncell[2]+r;
//printf("l has a value: %d\nnbr[l] has a value:%d\ni=%d j=%d k=%d\n",l,mycell[i][j][k]->nbr[l],i,j,k);

} } }
//printf("mycell[i][j][k] has %d neighbour\ni j k are:%d %d %d\n", mycell[i][j][k]->nnbr,i,j,k);     
   }
  }
 }
printf("i j k are  %d %d %d\n", i,j,k);
};     
//////////////////////////////////////////
//////////////////////////////////////////
////////Yongjian code ends////////////////
//////////////////////////////////////////
//////////////////////////////////////////


  SimulationBox() {ncell[0]=ncell[1]=ncell[2]=0; nmols=0;mycell=NULL;};
  
  ~SimulationBox() {
    int i,j,k;
    if(mycell!=NULL) {
#if VERBOSE == 1
      printf("Now destroy mycell.\n");
#endif
      for(i=0;i<ncell[0];i++) {
	for(j=0;j<ncell[1];j++) {
	  for(k=0;k<ncell[2];k++) {
	    delete mycell[i][j][k];
	  }
	  delete [] mycell[i][j];
	}
	delete [] mycell[i];
      }      
      delete [] mycell;
    }
  };
  void PrintInfo() {
    printf("SimulationBox: ncell(%d,%d,%d), L(%f,%f,%f), origin(%f,%f,%f),nmols(%d)\n",
	   ncell[0],ncell[1],ncell[2],
	   L[0],L[1],L[2],
	   origin[0],origin[1],origin[2],
	   nmols);
  };
};


/////////////////
// Subroutines //
/////////////////


////////////////////
// Main functions //
////////////////////

int main(int argc, char **argv)
{
  int i=0,j=0,k=0;
  char inputfilename[100], outputfilename[100];
  int inputparameter;
  FILE *ifp, *ofp;
  char templine[ONELINEMAX], str1[100], str2[100];
  long current_timestep, nmols;
  float xlo, xhi, ylo, yhi, zlo, zhi;
  int id, type;
  float ix, iy, iz, ivx, ivy, ivz, ifx, ify, ifz;
  SimulationBox **mysystem;
  Atom **myatom; //convenient pointer
  int iframe=0; //input number of frames
  int totalframes;
  float cellsize;

  /************* Command line parameter processing ************/

  if(argc!=5) {
  // Provide the correct syntax of this analysis code
    printf("Correct syntax: md.analysis.cxx inputfile outputfile inputparameter cellsize\n");
    exit(0);
  }
  else {
    //Note that argv[0] is "md.analysis.cxx"
    sscanf(argv[1], "%s", inputfilename);
    sscanf(argv[2], "%s", outputfilename);
    sscanf(argv[3], "%d", &inputparameter);
    sscanf(argv[4], "%f", &cellsize);
  }

  /************** Import the data file ********************/
  ifp = fopen(inputfilename, "r");
  ofp = fopen(outputfilename, "w");

  if(ifp==NULL) {
    printf("File %s does not exist!\n", inputfilename);
    exit(0);
  }

  //initialize Simulation Boxes
  mysystem = new SimulationBox *[MAXDUMP];
  for(i=0;i<MAXDUMP;i++) {mysystem[i] = new SimulationBox();}
printf("Initialization of Simulation Boxes has been finished! %d\n",i);
  //Read in the entire file
  while(!feof(ifp)) {
    fgets(templine, ONELINEMAX, ifp); // TIMESTEP line
    sscanf(templine, "%s %s", str1, str2);
#if VERBOSE==1
    //    printf("For the first line: (%s) (%s) (%s)\n", templine, str1, str2);
#endif
    if(strcmp(str1, "ITEM:")==0) {
      fgets(templine, ONELINEMAX, ifp); //actual timesteps
      sscanf(templine, "%ld", &current_timestep);
      fgets(templine, ONELINEMAX, ifp); // ITEM: NUMBER OF ATOMS
      fgets(templine, ONELINEMAX, ifp); // nmols
      sscanf(templine, "%ld", &nmols);  //number of atoms


      fgets(templine, ONELINEMAX, ifp); // ITEM: BOX BOUND
      fgets(templine, ONELINEMAX, ifp); // xlo xhi
      sscanf(templine, "%f %f", &xlo, &xhi);
      fgets(templine, ONELINEMAX, ifp); // ylo yhi
      sscanf(templine, "%f %f", &ylo, &yhi);
      fgets(templine, ONELINEMAX, ifp); // zlo zhi
      sscanf(templine, "%f %f", &zlo, &zhi);
      fgets(templine, ONELINEMAX, ifp); //empty line

      mysystem[iframe]->origin[0] = xlo;
      mysystem[iframe]->origin[1] = ylo;
      mysystem[iframe]->origin[2] = zlo;

      mysystem[iframe]->L[0] = xhi - xlo;
      mysystem[iframe]->L[1] = yhi - ylo;
      mysystem[iframe]->L[2] = zhi - zlo;

      //initialize myatom
      if(nmols>0) {
	mysystem[iframe]->nmols = nmols;
	mysystem[iframe]->myatom = new Atom *[nmols];
	myatom = mysystem[iframe]->myatom;
	for(i=0;i<nmols;i++) {
	  myatom[i] = new Atom();
	}
printf("Initialization of myatom  has been successfully finished! %d \n",i);
      }
      else {
	printf("Error in nmols (%d) \n", nmols);
	exit(0);
      }

      //Be aware that id may start from 1
      //id is the TRUE identification of atoms
      //your i index is NOT the indentification of atoms

#if VERBOSE == 1
      printf("Importing %d atoms...\n", nmols);
#endif
      for(i=0;i<nmols;i++) {
	fgets(templine, ONELINEMAX, ifp);
	sscanf(templine, "%d %d %f %f %f %f %f %f %f %f %f",
	       &id, &type, &ix, &iy, &iz, &ivx, &ivy, &ivz, &ifx, &ify, &ifz);
	myatom[i]->id = id; if (id==4000){printf("%d\n",id);}
	myatom[i]->type = type;

	myatom[i]->x[0] = ix;
	myatom[i]->x[1] = iy;
	myatom[i]->x[2] = iz;
      }

      //Setup the LinkCell
      mysystem[iframe]->SetupCell(cellsize);
      printf("%d\n",iframe);
      printf("mycell-111 nmember:%d nnbr:%d\n",mysystem[iframe]->mycell[1][1][1]->nmember,mysystem[iframe]->mycell[1][1][1]->nnbr); 
      //Any analysis routine



//////////////////////////////////////////////////////
//////Yongjian Coded the following part///////////////
//////////////////////////////////////////////////////


      //Throw in atoms, Coding here  
      float x1=0,x2=0,x3=0,originx,originy,originz,cellsizex,cellsizey,cellsizez;
      int p1,p2,p3,memberid;
      originx= mysystem[iframe]->origin[0]; printf("originx=%f\n",originx);
      originy= mysystem[iframe]->origin[1];
      originz= mysystem[iframe]->origin[2];
      cellsizex=(float)(mysystem[iframe]->L[0]/mysystem[iframe]->ncell[0]); printf("cellsizex=%f\n",cellsizex);
      cellsizey=(float)(mysystem[iframe]->L[1]/mysystem[iframe]->ncell[1]);
      cellsizez=(float)(mysystem[iframe]->L[2]/mysystem[iframe]->ncell[2]);
      //loop of every atoms of one frame in dump file 
      for (i=0;i<nmols;i++) {
         x1=myatom[i]->x[0]; // if (i==400){printf("%d\n",i);}
         x2=myatom[i]->x[1];
         x3=myatom[i]->x[2]; 
         //according to the atom's position, allocate them to the right cell. Since cellsizex is float number, to make it safe, when the index from atom's position is calculated from cellsizex, bring it back to within the cell. 
         p1=(int)((x1-originx)/cellsizex);
           if (p1==mysystem[iframe]->ncell[0]){p1=p1-1;} 
         p2=(int)((x2-originy)/cellsizey);
           if (p2==mysystem[iframe]->ncell[1]){p2=p2-1;} 
         p3=(int)((x3-originz)/cellsizez);
           if (p3==mysystem[iframe]->ncell[2]){p3=p3-1;}
         memberid=mysystem[iframe]->mycell[p1][p2][p3]->nmember;
//printf ("p1=%d p2=%d p3=%d\n",p1,p2,p3);
         //allocate the atom's id to the right cell 
         mysystem[iframe]->mycell[p1][p2][p3]->member[memberid]=myatom[i]->id;
//printf("mycell%d%d%d->member[%d]:%d\n",p1,p2,p3,memberid,mysystem[iframe]->mycell[p1][p2][p3]->member[memberid]);
      
//if (p1==0 && p2==0 && p3==0){printf("mycell-000 has atom(id): %d\n",mysystem[iframe]->mycell[p1][p2][p3]->member[mysystem[iframe]->mycell[p1][p2][p3]->nmember]);} 
         //increase the nmember of cell when a new atom is added into it.    
         mysystem[iframe]->mycell[p1][p2][p3]->nmember=mysystem[iframe]->mycell[p1][p2][p3]->nmember+1; 
//if(i==400){printf("nmember of mycell[%d][%d][%d] is %d\n",p1,p2,p3,mysystem[iframe]->mycell[p1][p2][p3]->nmember);}   
         }


//for (i=0;i<7;i++) {printf("mycell[0][%d][0] has nmember:%d\n",i,mysystem[iframe]->mycell[0][i][0]->nmember);}      





/////////////////////////////////////////////
      //Do RDF calculations, Coding here
/////////////////////////////////////////////

      float lx=mysystem[iframe]->L[0], ly=mysystem[iframe]->L[1],lz=mysystem[iframe]->L[2];

      float dr=0.01,dx,dy,dz,dist,area;

      int nbins=int(lx/dr)+1,index;
      int bin[nbins]; 
      float g[nbins];
//          printf("Done 1\n");
        //initialize the array!
        for(i=0;i<nbins;i++) {bin[i]=0;g[i]=0.0;}

                        //for (iframe=0, iframe<MAXDUMP,iframe++) //loop of different atoms in different cell;
			int cellx,celly,cellz,nmember,imember; //cell index
                        int ncellx,ncelly,ncellz; //neighbour cell index
                        int ncell[3],iii,ii,c,inmember;
			int comparei,comparej;

//printf("Done 2\n");
                        
                        for(iii=0;iii<3;iii++) {ncell[iii]=mysystem[iframe]->ncell[iii];} //initialize array ncell[3]
                        //loop among different cell 
                        for (cellx=0;cellx<ncell[0];cellx++){
			for (celly=0;celly<ncell[1];celly++){
			for (cellz=0;cellz<ncell[2];cellz++){

//printf("mysystem[%d]->mycell[%d][%d][%d]->nmember:%d\n",iframe,cellx,celly,cellz,mysystem[iframe]->mycell[cellx][celly][cellz]->nmember); 
                          
 
                        nmember=mysystem[iframe]->mycell[cellx][celly][cellz]->nmember;
                        //loop among all the atoms contained in one cell  
                                        for (imember=0;imember<nmember;imember++){ 
                                        //loop among all neighbour cells
                                        for (ii=0;ii<27;ii++){
					c=mysystem[iframe]->mycell[cellx][celly][cellz]->nbr[ii];
					
				        ncellx=(int)(c/(ncell[1]*ncell[2]));
					ncelly=(int)(fmod((int)(c/ncell[2]),ncell[1]));
					ncellz=(int)(fmod(c,ncell[2]));
//printf("mysystem[%d]->mycell[%d][%d][%d]: ncellx:%d ncelly:%d ncellz:%d\n",iframe,cellx,celly,cellz,ncellx,ncelly,ncellz); 
					inmember=mysystem[iframe]->mycell[ncellx][ncelly][ncellz]->nmember;
					//loop among all the neighbour's atoms
			         	for (j=0;j<inmember;j++){
					comparei=mysystem[iframe]->mycell[cellx][celly][cellz]->member[imember];
					comparej=mysystem[iframe]->mycell[ncellx][ncelly][ncellz]->member[j];
//printf("comparei:%d  comparej:%d\n",comparei,comparej);
//if (comparei==1066 && comparej==1067){printf("atom pair 1066-1067 has been included in the loops\n1066 in cell %d %d %d\n1067 in cell %d %d %d\n",cellx,celly,cellz,ncellx,ncelly,ncellz);}
					//Avoid double counting of pair comparei-comparej
                                        if (comparei<comparej) {
//printf("comparei:%d  comparej:%d\n",comparei,comparej);
					float xi,yi,zi,xj,yj,zj;
                                        xi=myatom[comparei-1]->x[0];xj=myatom[comparej-1]->x[0];
					yi=myatom[comparei-1]->x[1];yj=myatom[comparej-1]->x[1];
					zi=myatom[comparei-1]->x[2];zj=myatom[comparej-1]->x[2];
                                        dx=xi-xj;
//printf("dx=%f\n",dx);
                                        dy=yi-yj;
                                        dz=zi-zj;
                                        while (dx>lx/2.0) {dx=dx-lx;}
                                        while (dx<-lx/2.0){dx=dx+lx;}
                                        while (dy>ly/2.0) {dy=dy-ly;}
                                        while (dy<-ly/2.0){dy=dy+ly;}
                                        while (dz>lz/2.0) {dz=dz-lz;}
                                        while (dz<-lz/2.0){dz=dz+lz;}
                                        dist=sqrt(dx*dx+dy*dy+dz*dz);
//printf("Done 3\n");




                                        index=int(dist/dr);
//printf("index=%d\n",index);
                                        if(index<=nbins){bin[index]++;
//printf("bin%d=%d\n",index,bin[index]);
                                        }
					}
	              			}	
					} 

                            }
                        }}}

                               printf( "This is frame:(%d) \n", iframe);
                         float ndensity=float(nmols/(lx*ly*lz));
                         float norm=(1.0/nmols);
//printf("nmols:%d lx:%f ly:%f ly:%f ndensity:%f\n",nmols,lx,ly,lz,ndensity);
                         for (index=0;index<nbins;index++){
                              dist=index*dr;
                              if(index==0){g[index]=0;}
                              else {area=4*PI*dist*dist*dr;g[index]=2*bin[index]/(area*ndensity)*norm;}
                              printf("gofr %d %f %f %d \n",iframe,dist, g[index],bin[index]);

                              }

      //The simulation box information is in mysystem[iframe]. For instance mysystem[iframe]->L[0] will be X-dimension of the system
      //The atom information is in myatom. For instance myatom[i]->x[0] is the x-coordinate of myatom      
////////////////////////////////////////////////////
////////////////////////////////////////////////////
/////////// End of Yongjian's coding////////////////
////////////////////////////////////////////////////
////////////////////////////////////////////////////
      //Increment
      iframe++;
      if(iframe>=MAXDUMP) {
	printf("Too many dumps (%d), increase MAXDUMP!\n", iframe);
	exit(0);
      }
    
fprintf(ofp, "\nThis is frame: (%d) \nindex  distance  gof\n", iframe-1);
for (index=0;index<nbins;index++)
{
 dist=index*dr;
 fprintf (ofp, "%d %f %f \n",index,dist,g[index]);
  }

}
    else {
	if(!feof(ifp)) 
      printf("Not the right syntax, end of file.\n");
	else
      printf("Finish reading in input file\n");
    }
 } 

  totalframes = iframe;
  printf("Total number of frames is %d\n", totalframes);
  //  printf("Atom %f\n", mysystem[0]->myatom[1599]->x[0]);
  //  printf("atom: %f %f %f\n", mysystem[0]->myatom[1599]->x[0], mysystem[0]->myatom[1599]->x[1], mysystem[0]->myatom[1599]->x[2]);

  //mysystem[0]->PrintInfo();
  //output
  for(j=0;j<totalframes;j++) {
    for(i=0;i<mysystem[j]->nmols;i++) {
      if(mysystem[j]->myatom[i]->id==inputparameter)
      fprintf(ofp, "%d %f %f %f\n", j, 
	      mysystem[j]->myatom[i]->x[0],
	      mysystem[j]->myatom[i]->x[1],
	      mysystem[j]->myatom[i]->x[2]
	      );
    }
  }
 

  fclose(ifp);
  fclose(ofp);

  //deallocate memory
  for(j=0;j<totalframes;j++) {
    myatom = mysystem[j]->myatom;
    for(i=0;i<mysystem[j]->nmols;i++) {
      delete myatom[i];
    }
    delete [] myatom;    
  }
  for(j=0;j<MAXDUMP;j++) {
    delete mysystem[j];
  }
  delete [] mysystem;

}
