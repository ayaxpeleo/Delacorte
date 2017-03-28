#include <stdio.h>

#define SIZE 6


void print(long unsigned int m[SIZE][SIZE], long unsigned int mark1, long unsigned int mark2);
inline void copy(long unsigned int source[SIZE][SIZE], long unsigned int destiny[SIZE][SIZE]);
inline void swap(int i, int j, int p, int q, long unsigned int m[SIZE][SIZE]);
int distance(int i, int j, int p, int q);
int gcd(int m, int n);
inline long unsigned int eval(long unsigned int m[SIZE][SIZE]);
  long unsigned int d[SIZE][SIZE], small[SIZE][SIZE], big[SIZE][SIZE];

  long unsigned int total, current, distances[SIZE*SIZE][SIZE*SIZE], gcds[SIZE*SIZE][SIZE*SIZE];
  
  long unsigned int n[SIZE][SIZE] = { {12,  9, 5,   8, 10},
                                      {6,   7, 2,  11, 21},
                                      { 3,  1, 1,   1,  4},
                                      {8, 24, 1,   1, 15},
                                      {20, 22, 14, 18, 5}};

int main()
{
  int i, j;

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
  
  while(1)
  {
    printf("current matrix: \n");
    print(n, eval(n), 0);
    printf("swap: ");
    int x1, y1, x2, y2;
    scanf("%d %d %d %d", &x1, &y1, &x2, &y2);
    if(x1==-1 || y1==-1 || x2==-1 || y2==-1)
      break;
    swap(x1, y1, x2, y2, n);
  }
  
  return 0;
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


long unsigned int eval(long unsigned int m[SIZE][SIZE])
{
  int i, j;
  long unsigned int res;
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


void print(long unsigned int m[SIZE][SIZE], long unsigned int mark1, long unsigned int mark2)
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

void copy(long unsigned int source[SIZE][SIZE], long unsigned int destiny[SIZE][SIZE])
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

void swap(int i, int j, int p, int q, long unsigned int m[SIZE][SIZE])
{
  long unsigned int tmp;
  tmp = m[i][j];
  m[i][j]=m[p][q];
  m[p][q]=tmp;
}