/*---------------------------------------------------------------------------------- 
    Name       : Abhinav Bohra
    Roll No.   : 18CS30049
----------------------------------------------------------------------------------*/
#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#define MAX 10000
#define PI 3.141592653589793
/*------------------------------------------------------------------------------------------------------------------------------------------------*/
/*----- User-defined Structure Definitions -----*/

typedef struct point{
    long double x;                   //x co-ordinate of point
    long double y;                   //y co-ordinate of point 
}POINT;

typedef struct arc{
  POINT point;                      //Center point of arc  
  long double startAngle;           //Angle made by sector's left line segment with x-axis 
  long double endAngle;             //Angle made by sector's right line segment with x-axis
}ARC;

typedef struct tangent{
  POINT startpoint;                 //Starting point of tangent
  POINT endpoint;                   //Ending point of tangent
}TANGENT;

typedef struct stack{
  POINT *items;                     //Array storing elements of stack
  int top;                          //Top element in stack
  int size;                         //NUmber of elements currently present in stack
}STACK;
/*------------------------------------------------------------------------------------------------------------------------------------------------*/
/*----- Function Prototyping -----*/

void printcontzone(int u,int l,TANGENT *T,ARC *A);            //Function to print the boundaries (specifies by arcs & tangents) of the containment zone
void printCenters(POINT *S, int n);                           //Function to print the x & y co-ordinate of centers
void sortCenters(POINT *S, int n);                            //Function to sort centers w.r.t their x co-ordinates
void mergeSort(POINT *S, int low, int high);                  //Function to perform 'MERGE SORT'
void merge(POINT * S,int low, int mid, int high);             //Function to merge the sorted partitions
long double side(POINT p,POINT q,POINT r);                    //Function to compute orientation of a point r w.r.t line segment pq
long double computeAngle(POINT p, POINT q);                   //Function to compute angle made by line segment pq with x-axis
int CH(POINT *S,int n,int flag,POINT *H);                     //Function to determine the points in upper/lower hull
void contzone(POINT *UH,int u,POINT *LH,int l,long double r,TANGENT *T,ARC *A);  //Function to compute the boundaries (specifies by arcs & tangents) of the containment zone
/*----- Stack Functions Prototyping -----*/

STACK* createStack();            							  //Function to create an empty stack	
int isfull(STACK *s);                                         //Function to check whether a stack is full or not
int isempty(STACK *s);                                        //Function to check whether a stack is empty or not
POINT pop(STACK *s);                                          //Function to delete & return top element of stack
void push(STACK *s, POINT newitem);                           //Function to insert new element on top of stack
POINT top(STACK *s);                             			  //Function to fetch top most element of stack
POINT previousTop(STACK *s); 								  //Function to fetch 2nd top most element of stack
/*------------------------------------------------------------------------------------------------------------------------------------------------*/
/*----- Main Driver Function----*/

int main(){
    
    int n;             //Number of circles
    long double r;     //Radius of each circle
    POINT *S;          //Array S(of dtype 'POINT') to store center of n circles
    
    /*--------------------------------- Part 1: Input Handling  ------------------------------------------*/

    scanf("%d",&n);
    scanf("%Lf",&r);
    S = (POINT *)malloc((n+1)*(sizeof(POINT)));
    for(int i=1;i<=n;i++) scanf("%Lf %Lf",&S[i].x,&S[i].y);

    /*--------------------------------- Part 2: Sort the given points  ---------------------------------------------*/

    
    sortCenters(S,n);                                      				//Sort centers w.r.t their x co-ordinates
    printf("\n+++ Circles after sorting\n");        
    printCenters(S,n);                                    				//Print centers after sorting 
    
    /*--------------------------------- Part 3:  Compute the convex hull of the centers  --------------------------*/

    
    POINT *UH = (POINT *)malloc((n+1)*(sizeof(POINT)));    				//Array UH to store point in upper hull
    POINT *LH = (POINT *)malloc((n+1)*(sizeof(POINT)));    				//Array LH to store point in lower hull
    int n_UH = CH(S,n,0,UH);                               				//n_UH -> number of edges in upper hull
    int n_LH = CH(S,n,1,LH);                               				//n_LH -> number of edges in lower hull
    printf("+++ Upper Hull\n");
    printCenters(UH,n_UH + 1);                         					//Print points in upper hull
    printf("+++ Lower Hull\n");
    printCenters(LH,n_LH + 1);                            				//Print points in lower hull
    
    /*--------------------------------- Part 4: Compute the segments and the arcs of the boundary -----------------*/

    
    ARC *A = (ARC*)malloc((n_LH+1 + n_UH+1 + 1)*sizeof(ARC));           //Array A (of dtype 'ARC') to store arcs of the containment zone
    TANGENT *T = (TANGENT*)malloc((n_LH + n_UH + 1)*sizeof(TANGENT));   //Array T (of dtype 'TANGENT') to store tangents of the containment zone
    contzone(UH,n_UH,LH,n_LH,r,T,A);                                    //Compute containment zone
    printf("+++ The containment zone\n");      
    printcontzone(n_UH,n_LH,T,A);                                       //Print containment zone 

    return 0;
}
/*------------------------------------------------------------------------------------------------------------------------------------------------*/
/*----- Print Centers of Circle Function----*/

void printCenters(POINT *S, int n){

    /*
     * Arguments      : S-> Pointer to array S of POINTs
     *                  n-> number of points in S (size of S)
     * Task Performed : Prints x & y co-ordinates of points in S (space separated)
     * Returns        : void                
    */

    for (int i=1;i<=n; i++)
    {
        printf("    %0.15Lf %0.15Lf\n",S[i].x,S[i].y);
    }
    printf("\n");
    return; 
}
/*------------------------------------------------------------------------------------------------------------------------------------------------*/
/*-----Compute Convex Hull Function----*/
int CH(POINT *S,int n,int flag,POINT *H){

    /*
     * Arguments      : S-> Pointer to array of POINTs storing centers of circles
     *                  n-> number of points in S (size of S)
     *                  flag -> 0 denotes upper hull, 1 denotes lower hull
     *                  H -> Pointer to array of POINTs storing points in upper/lower hull
     * Task Performed : Determines points in upper/lower hull using 'GRAHAM SCAN' Algorithm
     * Returns        : n_edges -> number of edges in upper/lower hull                
    */
    
    STACK *T = createStack();    //Create an empty stack
    int n_edges=0;               //number of edges in upper/lower hull

    if(n<3) return n_edges;      //Base Case -> Hull Cannot be formed
    
    if(flag==0)   //Find points in Upper Hull
    {
        push(T,S[1]);
        push(T,S[2]);
        for(int i=3;i<=n;i++){
            while(T->size>=2 && side(previousTop(T),top(T),S[i]) > 0) {POINT temp = pop(T);}
            push(T,S[i]);
        }
        n_edges = T->size-1;
    }
    else        //Find points in Lower Hull
    {
        push(T,S[n]);
        push(T,S[n-1]);
        for(int i=n-2;i>=1;i--){
            while(T->size>=2 && side(previousTop(T),top(T),S[i]) > 0) {POINT temp = pop(T);}
            push(T,S[i]);
        }
        n_edges = T->size -1;
    }
    
    for(int i= n_edges+1;i>=1;i--) H[i]=pop(T);   //Populate Hull with points left in stack
    return n_edges;                               //Return number of edges in hull
}
/*------------------------------------------------------------------------------------------------------------------------------------------------*/
/*-----Containment Zone Computation Function----*/
void contzone(POINT *UH,int u,POINT *LH,int l,long double r,TANGENT *T,ARC *A){

    /*
     * Arguments      : UH -> Pointer to array of POINTs storing centers in upper hull
     *                  u  -> Number of edges in upper hull
     *                  LH -> Pointer to array of POINTs storing centers in lower hull
     *                  l  -> Number of edges in lower hull
     *                  r  -> Radius of each circle
     *                  T  -> Pointer to array of TANGENTs storing tangents of containment zone
     *                  A  -> Pointer to array of ARCs storing arcs of containment zone
     * Task Performed : Populates the arrays A and T using the below mentioned logic:-
     *                  Arc Computation:
     *                  Base Case -> startAngle(leftmost point) = PI, endAngle(rightmost point) = 0
     *                  For other points: startAngle of ith arc = endAngle of (i-1)th arc
     *                                    endAngle of ith arc   = Angle made by line segment joining ith & (i+1)th point
     *                                    taken in clockwise direction with the x-axis + PI/2
     *                  Note : In the lower section, all the angles are kept in the range [−PI,0] by subtracting 2*PI from positive angles
     *
     *                  Tangent Computation:
     *                  Base Case -> startAngle(leftmost point) = PI, endAngle(rightmost point) = 0
     *                  For other points: Starting point of ith tangent = Center of ith Arc + (r, endAngle(ith arc)  <- Polar Form
     *                  			      Ending point of ith tangent = Center of (i+1)th Arc + (r, endAngle(ith arc) <- Polar Form
     *                  (In Cartesian form : x_new = x + r*cos(theta), y_new = y + r*sin(theta) )    
     *                 
     * Returns        : void              
    */

    //Populating Arcs
    //Upper Hull -> Leftmost Point
    A[1].point=UH[1];
    A[1].startAngle = PI;
    A[1].endAngle   = computeAngle(UH[1],UH[2])+PI/2;

    //Upper Hull -> Others Points
    for(int i=2;i<=u;i++){
        A[i].point = UH[i];
        A[i].startAngle = A[i-1].endAngle;
        A[i].endAngle = computeAngle(UH[i],UH[i+1])+PI/2 ;
    }

    //Upper Hull -> Rightmost Point
    A[u+1].point=UH[u+1];
    A[u+1].startAngle = A[u].endAngle;
    A[u+1].endAngle   = 0.0;
    
    //Lower Hull -> Others Points
    for(int i=1;i<=l;i++){
        A[u+i+1].point = LH[i];
        A[u+i+1].startAngle = (A[u+i].endAngle <= PI)? A[u+i].endAngle : A[u+i].endAngle - 2*PI;
        A[u+i+1].endAngle = (computeAngle(LH[i],LH[i+1])+PI/2 <= PI)? computeAngle(LH[i],LH[i+1])+PI/2 : computeAngle(LH[i],LH[i+1])+PI/2- 2*PI;
    }

    //Lower Hull -> Rightmost Point
    A[u+l+2].point=LH[l+1];
    A[u+l+2].startAngle = A[u+l+1].endAngle;
    A[u+l+2].endAngle   = -PI;  

    //Populating Tangents
    for(int i=1,j=1;i<=u+l+1;i++,j++)
    {
        if(j==u+1) 
        {
            i--;
            continue;
        }
        
        T[i].startpoint.x = A[j].point.x + r*cos(A[j].endAngle);
        T[i].startpoint.y = A[j].point.y + r*sin(A[j].endAngle);
        T[i].endpoint.x   = A[j+1].point.x + r*cos(A[j].endAngle);
        T[i].endpoint.y   = A[j+1].point.y + r*sin(A[j].endAngle);   
       
     }

    T[u+l+1].endpoint.x   = A[1].point.x + r*cos(A[u+l].endAngle);
    T[u+l+1].endpoint.y   = A[1].point.y + r*sin(A[u+l].endAngle);

    return;
}
/*------------------------------------------------------------------------------------------------------------------------------------------------*/
/*-----Print Containment Zone Function----*/
void printcontzone(int u,int l,TANGENT *T,ARC *A){

    /*
     * Arguments      : u  -> Number of edges in upper hull
     *                  l  -> Number of edges in lower hull
     *                  T  -> Pointer to array of TANGENTs storing tangents of containment zone
     *                  A  -> Pointer to array of ARCs storing arcs of containment zone
     * Task Performed : Prints arcs & tangents of containment zone of upper section & lower section respectively
     * Returns        : void                
    */

    printf("\n--- Upper Section:\n\n");
    for(int i=1;i<=u+1;i++){
        printf("    Arc     : (%0.15Lf, %0.15Lf) From %0.15Lf to %0.15Lf\n",A[i].point.x,A[i].point.y,A[i].startAngle,A[i].endAngle);
        if(i==u+1) continue;
        printf("    Tangent : From (%0.15Lf, %0.15Lf) to (%0.15Lf, %0.15Lf)\n",T[i].startpoint.x,T[i].startpoint.y,T[i].endpoint.x,T[i].endpoint.y);
    }

    printf("\n--- Lower Section:\n\n");
    for(int i=u+2;i<=u+l+2;i++){
        printf("    Arc     : (%0.15Lf, %0.15Lf) From %0.15Lf to %0.15Lf\n",A[i].point.x,A[i].point.y,A[i].startAngle,A[i].endAngle);
        if(i==u+l+2) continue;
        printf("    Tangent : From (%0.15Lf, %0.15Lf) to (%0.15Lf, %0.15Lf)\n",T[i-1].startpoint.x,T[i-1].startpoint.y,T[i-1].endpoint.x,T[i-1].endpoint.y);
    }
    printf("\n\n");
}
/*------------------------------------------------------------------------------------------------------------------------------------------------*/
/*-----Utility Functions for mathematics----*/
long double side(POINT p,POINT q,POINT r){
	
	/*
     * Arguments      : Points p,q,r
     * Task Performed : Determines the orientation of point r w.r.t line segment pq using the following determinant
     *
     *        					|  1    1    1  |
     * 					side =  | p.x  q.x  r.x |
     *         					| p.y  q.y  r.y |
     *
     * Returns        : return 'side'                
    */
    long double side = (q.x*r.y - r.x*q.y) - (p.x*r.y - r.x*p.y) + (p.x*q.y - q.x*p.y);
    return side;
}

long double computeAngle(POINT p, POINT q){
    
    /*
     * Arguments      : Points p,q,r
     * Task Performed : Computes angle made by line segment pq with x-axis
     * Returns        : returns computed angle 'theta'              
    */
    long double X = q.x - p.x;
    long double Y = q.y - p.y;
    long double theta = atan2(Y,X);
    return theta;
}

/*------------------------------------------------------------------------------------------------------------------------------------------------*/
/*----- MergeSort Function----*/

void sortCenters(POINT *S, int n){
	
	/*
     * Arguments      : S-> Pointer to array S of POINTs
     *                  n-> number of points in S (size of S)
     * Task Performed : Sorts the array S w.r.t x co-ordinates of points in O(n*log(n)) time using 'MERGE SORT' algorithm
     * Returns        : void                
    */
    
    mergeSort(S,1,n);
    return ;
}

void mergeSort(POINT *S, int low, int high){

    if(low < high){
        int mid = (low + high)/2;
        mergeSort(S,low,mid);
        mergeSort(S,mid+1, high);
        merge(S,low,mid,high);
    }
}
void merge(POINT * S,int low, int mid, int high){

    int n1 = mid - low + 1;
    int n2 = high - mid;
    POINT *L = (POINT *)malloc((n1)*sizeof(POINT));
    POINT *R = (POINT *)malloc((n2)*sizeof(POINT));
 
    for(int i = 0; i < n1; i++) L[i] = S[low + i];
    for(int j = 0; j < n2; j++) R[j] = S[mid + 1 + j];
 
    int i = 0, j = 0,k = low;
     
    while (i < n1 && j < n2)
    {
        if (L[i].x <= R[j].x) S[k++] = L[i++];
        else                  S[k++] = R[j++];
    }
    while (i < n1) S[k++] = L[i++];
    while (j < n2) S[k++] = R[j++];

}
/*------------------------------------------------------------------------------------------------------------------------------------------------*/
/*-----Stack Functions----*/

STACK* createStack() {
  
  STACK* s =(STACK*)malloc(sizeof(STACK));
  s->top = -1;
  s->size=0;
  s->items=(POINT*)malloc((MAX)*sizeof(POINT));
  return s;
}

int isfull(STACK *s) {
  if (s->top == MAX - 1)  return 1;
  else                    return 0;
}

int isempty(STACK *s) {
  if (s->top == -1)  return 1;
  else               return 0;
}

void push(STACK *s,POINT newitem) {
  if (isfull(s))   printf("STACK FULL");
  else 
  {
    s->top++;
    s->items[s->top] = newitem;
    s->size++;
  }
}

POINT pop(STACK *s) {
  POINT point;
  point.x=-1,point.y=-1;

  if (isempty(s)) return point;
  else            point = s->items[s->top--];

  s->size--;
  return point;
}

POINT top(STACK *s){
    POINT point;
    point.x=-1,point.y=-1;

    if(isempty(s)) return point;
    else           return s->items[s->top];
}

POINT previousTop(STACK *s){
    POINT point;
    point.x=-1,point.y=-1;

    if(s->size<=1) return point;
    else           return s->items[(s->top)-1];
}
/*------------------------------------------------------------------------------------------------------------------------------------------------*/
