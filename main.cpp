#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#define SIZE 7

void find(long int m[SIZE][SIZE], long int *min, long int *max, int index, long int total);
void print(long int m[SIZE][SIZE], long int mark1, long int mark2);
inline void copy(long int source[SIZE][SIZE], long int destiny[SIZE][SIZE]);
inline void swap(int i, int j, int p, int q, long int m[SIZE][SIZE]);
int distance(int i, int j, int p, int q);
int gcd(int m, int n);
inline long int eval(long int m[SIZE][SIZE]);
long int eval3(long int m[SIZE][SIZE], long int swap1, long int swap2);
inline void delete_val(long int m[SIZE][SIZE], long int n);
int check_consistency(long int m[SIZE][SIZE]);
long int eval_single(long int m[SIZE][SIZE], int p, int q);
//long int eval3(long int m[SIZE][SIZE]);
  long int d[SIZE][SIZE], small[SIZE][SIZE], big[SIZE][SIZE];
  
  clock_t init, totaltime;  
  long int total, current, distances[SIZE*SIZE][SIZE*SIZE], gcds[SIZE*SIZE][SIZE*SIZE];
  
  long int fixed[SIZE][SIZE] =  { {1,  0,  0,  0, 0, 0, 0},
                                  {0,  0,  0,  0, 0, 0, 0},
                                  {0,  0,  0,  1, 0, 0, 0},
                                  {0,  0,  1,  1, 1, 0, 0},
                                  {0,  0,  0,  1, 1, 0, 0},
                                  {0,  0,  0,  0, 0, 0, 0},
                                  {0,  0,  0,  0, 0, 0, 1}};
                       
  
int main()
{
  int i, j, k, l;

  for(i=0;i<SIZE*SIZE;i++)
  {
    for(j=0;j<SIZE*SIZE;j++)
    {
      distances[i][j] = distance(i/SIZE, i%SIZE, j/SIZE, j%SIZE);
    }
  }
  
  for(i=0;i<SIZE*SIZE;i++)
  {
    for(j=0;j<SIZE*SIZE;j++)
    {
      gcds[i][j] = gcd(i+1, j+1);
    }
  }
  
  //total = (SIZE*SIZE)*(SIZE*SIZE-1)*(SIZE*SIZE-2)*(SIZE*SIZE-3)*(SIZE*SIZE-4)*(SIZE*SIZE-5);
  current = 0;
  //long int n[2][3] = {{4, 1, 3},{6, 5, 2}};
  //printf("Evaluada, la matriz da: %lu\n", eval_max(n));
  //return 0;
  long int m[SIZE][SIZE], min, max;
/*
  long int n[SIZE][SIZE] = { {16,  2,  3,  4,  5},
                            {6,   7,  11, 9,  10},
                            {8,   1,  1,  1,  22},
                            {12,  25, 1,  1,  21},
                            {20,  18, 14, 15, 24}};
  */
/*  
  long int n[SIZE][SIZE] = {{12,  9,  5,   8, 10},
                            { 6,  7,  2,  11, 21},
                            { 3,  1,  1,   1,  4},
                            { 8, 24,  1,   1, 15},
                            {20, 22, 14,  18,  5}};
*/
  /*
  long int n[SIZE][SIZE] = { {18,  2,  3,  4,  5},
                                      {6,   7, 8, 9,  10},
                                      {11,  1, 1, 1, 22},
                                      {16,  25, 1,  1, 21},
                                      {20,  15, 14, 12, 24}};
    */                                  
    
  long int n[SIZE][SIZE] = {  { 48, 44,  8,  9, 14, 30, 36},
                              { 40, 26, 11, 34, 25,  6, 21},
                              { 15, 19,  1, 29,  2, 23, 28},
                              { 27,  7, 31, 37, 41,  5, 39},
                              { 35,  4,  3, 43, 47, 38, 32},
                              { 12, 46, 13, 17,  49, 33, 20},
                              { 42, 18, 10, 16, 22, 45, 24}};
                                       
  for(i=0;i<SIZE;i++)
  {
    for(j=0;j<SIZE;j++)
    {
      d[i][j]=0;
    }
  }
  for(i=0;i<SIZE*SIZE;i++)
  {
    for(j=0;j<SIZE*SIZE;j++)
    {
      long int dist = distance(i/SIZE, i%SIZE, j/SIZE, j%SIZE);
      //printf("The distance from %d,%d to %d,%d is : %lu\n", i/SIZE, i%SIZE, j/SIZE, j%SIZE, dist);
      d[i/SIZE][i%SIZE]+=dist;
    }
  }
  print(d, 0, 0);
  
  //we initialize the matrix
  for( i=0; i<SIZE ;i++)
  {
    for( j=0; j<SIZE; j++)
    {
      m[i][j] = i*SIZE+j+1;
    }
  }
  
  //eliminate numbers with gcd of 1 (as of others)
  for( i=0; i<SIZE ;i++)
  {
    for( j=0; j<SIZE; j++)
    {
      int flag_ok=0;
      for(k=1;k<SIZE*SIZE;k++)
      {
        if(m[i][j]==(k+1))
          continue;
        if(gcd(m[i][j], k+1)!=1)
        {
          flag_ok=1;
        }
      }
      if(flag_ok==0)
        m[i][j]=1;
    }
  }
  //eliminate unnecesary numbers
  int u_numbers=0;
  for( i=0; i<SIZE ;i++)
  {
    for( j=0; j<SIZE; j++)
    {
      int check_number = m[i][j]-1;
      int flag_out=0;
      //si hay un multiplo del numero, es mejor salirse
      
      for(k=1;k<=SIZE*SIZE;k++)
      {
        if(k%m[i][j]==0 && k!=m[i][j])
          flag_out=1;
      }
      if(flag_out)
        continue;
      
      for(k=(SIZE*SIZE-1);k>0;k--)
      {
        if(k==check_number)
          continue;
        for(l=1, flag_out=0;l<(SIZE*SIZE-2) && flag_out==0; l++)
        {
          //printf("checking %d, k: %d, l: %d, gcds[check_number][l]: %lu, gcds[k][l]: %lu\n", check_number+1, k+1, l, gcds[check_number][l], gcds[k][l]);
          if(l==check_number)
            continue;
          if(gcds[check_number][l] != gcds[k][l])
          {
            flag_out=1;
          }
        }  
        if(flag_out==0 && m[i][j]>(k+1))
        {
          printf("#%d: %ld = %d\n", u_numbers+1, m[i][j], k+1);
          m[i][j]=k+1;
          u_numbers++;
          break;
        }
      }   
    }
  }
  printf("\n%d redundant numbers were modified...\n", u_numbers);
  
  
  
  init = clock();
  total = clock();
  min=eval(n);
  max=eval(n);
  printf("At the start, the min value is %lu, and the max value is %lu.\n", min, max);
  find(n, &min, &max, 0, eval(m));
  printf("\nPrinting biggest: \n");
  print(big, max, 0);
  printf("\nPrinting smallest: \n");
  print(small, 0, min);
  //eval2(small);
  return 0;
}


void find(long int m[SIZE][SIZE], long int *min, long int *max, int index, long int total)
{
  long int swap1, swap2, newtotal, eval1, total_tmp;
  total_tmp=total;
  
  if (total>(*max))
  {
    (*max)=total;
    copy(m, big);
    printf("big found... %lu\n", *max);
    print(big, *max, 0);
    
    clock_t t = clock()-init;
    init = clock();
    float ti  = ((float)t)/CLOCKS_PER_SEC;
    float total_t  = ((float)(clock()-totaltime))/CLOCKS_PER_SEC;
    total_t = total_t/60;
    printf ("It took me %.3f s, total: %.3f m (%d)\n", ti, total_t, SIZE);
  }
  if (total<(*min))
  {
    (*min)=total;
    copy(m, small);
    printf("small found... %lu\n", *min);
    print(small, 0, *min);
    clock_t t = clock()-init;
    init = clock();
    float ti  = ((float)t)/CLOCKS_PER_SEC;
    float total_t  = ((float)(clock()-totaltime))/CLOCKS_PER_SEC;
    total_t = total_t/60;
    printf ("It took me %.3f s, total: %.3f m (%d)\n", ti, total_t, SIZE);
  }
  int i, j;
  i = index;
  //long int mtable[SIZE][SIZE];
  //copy(d, mtable);
  //delete_val(mtable, mtable[i/SIZE][i%SIZE]);
  eval1 = eval_single(m, i/SIZE, i%SIZE);
  for(j=i+1;j<=(SIZE*SIZE-1);j++)
  {
    //printf("Now checking j: %d, i: %d, d[i/SIZE][i%%SIZE]: %lu, d[j/SIZE][j%%SIZE]: %lu\n", j, i, d[i/SIZE][i%SIZE], d[j/SIZE][j%SIZE]);
    //if((mtable[j/SIZE][j%SIZE]==0 && j!=i) && ((d[i/SIZE][i%SIZE]==d[j/SIZE][j%SIZE] || m[i/SIZE][i%SIZE]==m[j/SIZE][j%SIZE]) || d[m[j/SIZE][j%SIZE]/SIZE][m[j/SIZE][j%SIZE]%SIZE]==d[j/SIZE][j%SIZE] || d[m[i/SIZE][i%SIZE]/SIZE][m[i/SIZE][i%SIZE]%SIZE]==d[i/SIZE][i%SIZE]))
    //if(fixed[i/SIZE][i%SIZE])
    //  i++;
    if(m[i/SIZE][i%SIZE]==m[j/SIZE][j%SIZE] || fixed[j/SIZE][j%SIZE] || fixed[i/SIZE][i%SIZE])
      find(m, min, max, i+1, total_tmp);
    else
    {
      swap1 = eval1 + eval_single(m, j/SIZE, j%SIZE);
      swap(i/SIZE, i%SIZE, j/SIZE, j%SIZE, m);
      swap2 = eval_single(m, i/SIZE, i%SIZE) + eval_single(m, j/SIZE, j%SIZE);
      newtotal = total_tmp - swap1 + swap2;
      //if(check_consistency(m))
      find(m, min, max, i+1, newtotal);
      swap(i/SIZE, i%SIZE, j/SIZE, j%SIZE, m);
      //delete_val(mtable, mtable[j/SIZE][j%SIZE]);
    }
  }
  return;
}

void print(long int m[SIZE][SIZE], long int mark1, long int mark2)
{
  int i, j;
  printf("\n");
  for(i=0;i<SIZE;i++)
  {
    for(j=0;j<SIZE;j++)
    {
      if(j==0) printf("(");
      printf("\t%lu", m[i][j]);
      if(j!=SIZE-1) printf(",");
    }
    if(i!=SIZE-1) printf("),");
    else  printf(")");
    if(mark1>0 && i==0)
      printf("\t mark1: %lu", mark1);
    if(mark2>0 && i==1)
      printf("\t mark2: %lu", mark2);
    printf("\n");
  }
}

void copy(long int source[SIZE][SIZE], long int destiny[SIZE][SIZE])
{
  int i, j;
  for(i=0;i<SIZE;i++)
  {
    for(j=0;j<SIZE;j++)
    {
      destiny[i][j]=source[i][j];
    }
  }
}

void swap(int i, int j, int p, int q, long int m[SIZE][SIZE])
{
  long int tmp;
  tmp = m[i][j];
  m[i][j]=m[p][q];
  m[p][q]=tmp;
}

//long int eval_max(long int m[SIZE][SIZE])
long int eval(long int m[SIZE][SIZE])
{
  int i, j;
  long int res;
  //printf("a\tb\tgcd(a, b)\tDa,b\n");
  for(res=0, i=0;i<(SIZE*SIZE);i++)
  {
    for(j=(i+1);j<(SIZE*SIZE);j++)
    {
      //printf("%lu\t%lu\t      %d \t  %d\n", m[i/SIZE][i%SIZE], m[j/SIZE][j%SIZE], gcd(m[i/SIZE][i%SIZE], m[j/SIZE][j%SIZE]), distance(i/SIZE, i%SIZE, j/SIZE, j%SIZE));
      //res += (gcd(m[i/SIZE][i%SIZE], m[j/SIZE][j%SIZE])*distance(i/SIZE, i%SIZE, j/SIZE, j%SIZE));
      res += gcds[m[i/SIZE][i%SIZE]-1][m[j/SIZE][j%SIZE]-1]*distances[i][j];
    }
  }
  return res;
}
/*
long int eval3(long int m[SIZE][SIZE], long int swap1, long int swap2)
{
  int i, j;
  long int res;
  //printf("a\tb\tgcd(a, b)\tDa,b\n");
  for(res=0, i=0;i<(SIZE*SIZE);i++)
  {
    for(j=(i+1);j<(SIZE*SIZE);j++)
    {
      //printf("%lu\t%lu\t      %d \t  %d\n", m[i/SIZE][i%SIZE], m[j/SIZE][j%SIZE], gcd(m[i/SIZE][i%SIZE], m[j/SIZE][j%SIZE]), distance(i/SIZE, i%SIZE, j/SIZE, j%SIZE));
      //res += (gcd(m[i/SIZE][i%SIZE], m[j/SIZE][j%SIZE])*distance(i/SIZE, i%SIZE, j/SIZE, j%SIZE));
      res += gcds[m[i/SIZE][i%SIZE]-1][m[j/SIZE][j%SIZE]-1]*distances[i][j];
    }
  }
  return res;
}
*/


long int eval2(long int m[SIZE][SIZE])
{
  int i, j;
  long int res;
  printf("a\tb\tgcd(a, b)\tDa,b\n");
  for(res=0, i=0;i<(SIZE*SIZE);i++)
  {
    for(j=(i+1);j<(SIZE*SIZE);j++)
    {
      printf("%lu\t%lu\t      %d \t  %d\n", m[i/SIZE][i%SIZE], m[j/SIZE][j%SIZE], gcd(m[i/SIZE][i%SIZE], m[j/SIZE][j%SIZE]), distance(i/SIZE, i%SIZE, j/SIZE, j%SIZE));
      res += (gcd(m[i/SIZE][i%SIZE], m[j/SIZE][j%SIZE])*distance(i/SIZE, i%SIZE, j/SIZE, j%SIZE));
    }
  }
  return res;
}

int gcd(int m, int n)
{
  if(m<2 || n<2)
    return 1;
  
  int cd = 1;
  
  int i;
  for(i=2;i<=m;)
  {
    if(m%i==0 && n%i==0)
    {
      cd = cd*i;
      m = m/i;
      n = n/i;
    }
    else
      i++;
  }
  return cd;
}

int distance(int i, int j, int p, int q)
{
  int tmp1, tmp2;
  tmp1 = (p-i)*(p-i);
  tmp2 = (q-j)*(q-j);
  return tmp1+tmp2;
}

void delete_val(long int m[SIZE][SIZE], long int n)
{
  int i, j;
  for(i=0;i<SIZE;i++)
  {
    for(j=0;j<SIZE;j++)
    {
      if(m[i][j]==n)
        m[i][j]=0;
    }
  }
}

int check_consistency(long int m[SIZE][SIZE])
{
    //if there are 1's bot at the outer edge and at any point inside that edge, this shouldn't be a good path
    int i,j, edge, center;
    for(i=0, edge=0, center=0; i<SIZE; i++)
    {
      for(j=0;j<SIZE;j++)
      {
        if(m[i][j]==1)
        {
          if((i==0 || i==(SIZE-1) || j==0 || j==(SIZE-1)) && m[i][j]==1)
            edge=1;
          else
            center=1;
        } 
      }
    }
    return ((edge+center)==1);
}
long int eval_single(long int m[SIZE][SIZE], int p, int q)
{
  int i;
  long int res;
  //printf("a\tb\tgcd(a, b)\tDa,b\n");
  for(res=0, i=0;i<(SIZE*SIZE);i++)
  {
    if((p*SIZE+q)==i)
      continue;
    res += gcds[m[p][q]-1][m[i/SIZE][i%SIZE]-1]*distances[p*SIZE+q][i];
  }
  return res;
}


/*
    0,0   0,1   0,2   0,3   0,4
    1,0   1,1   1,2   1,3   1,4
    2,0   2,1   2,2   2,3   2,4
    3,0   3,1   3,2   3,3   3,4
    4,0   4,1   4,2   4,3   4,4
*/


