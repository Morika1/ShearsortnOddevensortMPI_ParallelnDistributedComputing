// Authors: Sharon Tavdy 315856880, Mor Tavdy 315856898
#include <mpi.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <sys/types.h>
#include <unistd.h>

#define TAG 0
#define NDIMS 2
#define FILE_NAME "Rectangles.dat.txt"
#define FINAL_FILE_NAME "Result.dat.txt"
enum Tags {ROOT, MIN, MAX};


typedef struct {
int id;
double height;
double width;
double area;
}Rectangle;

Rectangle* allRectangles;
int rows, cols;

void readFromFile(const char* fileName, int size);
void create(int rank, Rectangle *my_rectangle, int size, MPI_Comm cart_comm, Rectangle* finalArray, MPI_Datatype RectangleMPIType, const char* fileName);
void shearSort(int* coordinates, Rectangle* my_rectangle, MPI_Comm cart_comm, MPI_Datatype RectangleMPIType);
void oddEven(int loopSize, int neighbour1, int neighbour2, int* coordinates, int coordinate, MPI_Comm cart_comm, Rectangle* my_rectangle, int i, MPI_Datatype RectangleMPIType);
void switchRectangles(int neighbour,Rectangle* my_rectangle, MPI_Comm cart_comm, enum Tags switchOption, MPI_Datatype RectangleMPIType, int* coordinates);
void testAndSwitch(int* coordinates, int i, int neighbour,Rectangle* my_rectangle, MPI_Comm cart_comm,enum Tags switchOption1, enum Tags switchOption2, MPI_Datatype RectangleMPIType);
void printAndSaveArray(Rectangle* finalArray, const char* fileName);


int main(int argc, char* argv[])
{
    int nprocs, rank; 
    Rectangle rectangle; 
    Rectangle finalArray[nprocs];
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm cart_comm;
    MPI_Datatype RectangleMPIType;
    MPI_Datatype type[4] = { MPI_INT, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE};
    int blockLen[4] = {1, 1, 1, 1};
    MPI_Aint disp[4];

    disp[0] = (char*)&rectangle.id - (char*)&rectangle;
    disp[1] = (char*)&rectangle.height - (char*)&rectangle;
    disp[2] = (char*)&rectangle.width - (char*)&rectangle;
    disp[3] = (char*)&rectangle.area - (char*)&rectangle;
    MPI_Type_create_struct(4, blockLen, disp, type, &RectangleMPIType);
    MPI_Type_commit(&RectangleMPIType); // creatring MPIDataType to pass processes struct

    if(rank==ROOT)
        readFromFile(FILE_NAME, nprocs);

    MPI_Scatter(allRectangles,1, RectangleMPIType, &rectangle, 1, RectangleMPIType, TAG, MPI_COMM_WORLD); // each process takes one rectangle
    
    rows=log2(nprocs); //do square root
    cols=log2(nprocs);
    int size=  rows*cols;

    create(rank, &rectangle, size, cart_comm,finalArray, RectangleMPIType, FINAL_FILE_NAME); // to create grid in cartesian communicator and sort

     if(rank==0)
        free(allRectangles);

     MPI_Finalize();

    return 0;
}

void readFromFile(const char* fileName, int size)
{
    FILE* fp;
    // Open file for reading rectangles
    if ((fp = fopen(fileName, "r")) == 0) 
    {
        printf("cannot open file %s for reading\n", fileName);
        exit(0);
    }
    // Allocate array of rectangles end Read data from the file
    allRectangles = (Rectangle*)malloc( size* sizeof(Rectangle));
    if (allRectangles == NULL) 
    {
        printf("Problem allocating memory\n");
        MPI_Abort(MPI_COMM_WORLD, 1);
    
    }
    //can suppose nprocs= numberOfRectangles
    for (int i = 0; i < size; i++) 
    {
        fscanf(fp, "%d %lf %lf", &allRectangles[i].id, &allRectangles[i].height, &allRectangles[i].width);
        allRectangles[i].area=allRectangles[i].height*allRectangles[i].width;
    }
    fclose(fp);
}

void create(int rank, Rectangle *rectangle, int size, MPI_Comm cart_comm, Rectangle* finalArray, MPI_Datatype RectangleMPIType, const char* fileName )
{
    int dims[NDIMS]= {rows, cols};
    int periods[NDIMS]={0,0};
    int reorder=0;

    MPI_Cart_create(MPI_COMM_WORLD, NDIMS, dims, periods, reorder, &cart_comm);
    int coordinates[NDIMS]; 
    MPI_Cart_coords(cart_comm, rank, NDIMS, coordinates);

    Rectangle my_rectangle=*rectangle;
    shearSort(coordinates, &my_rectangle, cart_comm, RectangleMPIType); 

    (*rectangle)=my_rectangle; 
    MPI_Barrier(cart_comm);
    MPI_Gather(rectangle, 1, RectangleMPIType,finalArray, 1, RectangleMPIType, TAG, cart_comm); // each process sends to procs 0 including procs 0 the rectangle he holds after sorting

    if(rank==0)
    {
        printf("Final array of ids sorted by rectangles area:\n");
        printAndSaveArray(finalArray, FINAL_FILE_NAME);
    }
}

void shearSort(int* coordinates, Rectangle* my_rectangle, MPI_Comm cart_comm, MPI_Datatype RectangleMPIType)
{
    int numOfIterations = 2 * ((int)log2(cols)) + 1;
    int left, right, up, down;
    int displacement=1;

    MPI_Cart_shift(cart_comm, 1, displacement, &left, &right); // left and right neighbours of process in grid
    MPI_Cart_shift(cart_comm, 0, displacement, &up, &down); // up and down neighbours of process in grid

    for(int i=0; i<numOfIterations; i++)
    {
        if(i%2==0) // row sort
            oddEven(cols, left, right, coordinates, coordinates[1],  cart_comm, my_rectangle, i, RectangleMPIType);
        
        else // col sort
            oddEven(rows, up, down,coordinates, coordinates[0], cart_comm, my_rectangle, i, RectangleMPIType);
    }
}

void oddEven(int loopSize, int neighbour1, int neighbour2, int* coordinates, int coordinate, MPI_Comm cart_comm, Rectangle* my_rectangle, int i, MPI_Datatype RectangleMPIType)
{
    for(int j=0; j<loopSize; j++)
    {
        int neighbour;
        if(j%2==0) // rows: even cols talk right, cols: even rows talk down
        {
            if(coordinate%2==0) // rows: if col in even , cols: if row is even
            {
                neighbour= neighbour2; // right or down
                testAndSwitch(coordinates, i, neighbour, my_rectangle, cart_comm, MIN, MAX, RectangleMPIType);
            }
            else // rows: if col is odd, cols: if row is odd
            {
                neighbour= neighbour1; // left or up
                testAndSwitch(coordinates, i, neighbour, my_rectangle, cart_comm, MAX, MIN, RectangleMPIType);
            }
        }
        else // rows: even cols talk left, cols: even rows talk up
        {
             if(coordinate%2==0) // rows: if col in even , cols: if row is even
            {
                neighbour= neighbour1; // left or up
                if(neighbour!=-2)
                    testAndSwitch(coordinates, i, neighbour, my_rectangle, cart_comm, MAX, MIN, RectangleMPIType);        
            }
            else // rows: if col is odd, cols: if row is odd
            {
                neighbour= neighbour2; // right or down
                if(neighbour!=-2)
                    testAndSwitch(coordinates, i, neighbour, my_rectangle, cart_comm, MIN, MAX, RectangleMPIType);
            }
        }
    }
}

void switchRectangles(int neighbour,Rectangle* my_rectangle, MPI_Comm cart_comm, enum Tags switchOption, MPI_Datatype RectangleMPIType, int* coordinates)
{
    Rectangle other_rectangle;
    MPI_Sendrecv(my_rectangle,1,RectangleMPIType,neighbour,TAG, &other_rectangle,1,RectangleMPIType,neighbour,TAG,cart_comm,MPI_STATUS_IGNORE);
    double my_area= my_rectangle->area;
    double other_area= other_rectangle.area;
    int my_id = my_rectangle->id;
    int other_id = other_rectangle.id;

    if(switchOption == MIN)
    {
      if( my_area < other_area || (my_area == other_area && my_id< other_id ))
      {
          my_rectangle->id= other_rectangle.id;
          my_rectangle->height = other_rectangle.height;
          my_rectangle->width = other_rectangle.width;
          my_rectangle->area= other_rectangle.area;  
      }
    }

    if(switchOption == MAX)
    {
        if( my_area > other_area || (my_area == other_area && my_id > other_id ))
      {
          my_rectangle->id= other_rectangle.id;
          my_rectangle->height = other_rectangle.height;
          my_rectangle->width = other_rectangle.width;
          my_rectangle->area= other_rectangle.area;  
      }
    }
}

void testAndSwitch(int* coordinates, int i, int neighbour,Rectangle* my_rectangle, MPI_Comm cart_comm,enum Tags switchOption1, enum Tags switchOption2, MPI_Datatype RectangleMPIType)
{
    if(i%2==0 && coordinates[0]%2!=0) // descending row sort
        switchRectangles(neighbour, my_rectangle, cart_comm, switchOption1,RectangleMPIType, coordinates); 

    else
        switchRectangles(neighbour, my_rectangle, cart_comm, switchOption2,RectangleMPIType, coordinates); 
}

void printAndSaveArray(Rectangle* finalArray, const char* fileName)
{
    FILE* fp;
    fp= fopen(FINAL_FILE_NAME, "w+");
    if(!fp)
    {
        printf("Problem with opening file\n");
        exit(1);
    }
    for(int i=0; i<rows; i++)
    {
        if(i%2==0)
        {
            for( int j=0; j< cols; j++ )
            {
                printf("%d ", finalArray[i*cols+j].id);
             fprintf(fp, "%d ",finalArray[i*cols+j].id );
            }
             
        }
        else
        {
            for( int j=cols-1; j>=0; j-- )
            {
                 printf("%d ", finalArray[i*cols+j].id);
                 fprintf(fp, "%d ", finalArray[i*cols+j].id);
            }
               
        } 
    }
    fclose(fp);
}

