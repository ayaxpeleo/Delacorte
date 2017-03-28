#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#define SIZE 7

void find(int m[SIZE][SIZE], long int *record, int index, int limit, long int total);
void print(int m[SIZE][SIZE], long int val);
inline void copy(int source[SIZE][SIZE], int destiny[SIZE][SIZE]);
inline void swap(int i, int j, int p, int q, int m[SIZE][SIZE]);
int distance(int i, int j, int p, int q);
int gcd(int m, int n);
inline long int eval(int m[SIZE][SIZE]);
long int eval_single(int m[SIZE][SIZE], int p, int q);
int d[SIZE][SIZE], small[SIZE][SIZE], big[SIZE][SIZE];
clock_t init, init2;  
clock_t totaltime;
long int total, current;
int distances[SIZE*SIZE][SIZE*SIZE], gcds[SIZE*SIZE][SIZE*SIZE];
void shuffle(int *array, size_t n);
  
  int debug;
  
class Int
{
	short int n;
	unsigned char m1;
	struct
	{
		unsigned int f1	:1;
		unsigned int f2	:1;
		unsigned int f3	:1;
		unsigned int f4	:1;
		unsigned int f5	:1;
		unsigned int f6	:1;
		unsigned int f7	:1;
		unsigned int f8	:1;
	};
};
class M
{
	Int *m[SIZE*SIZE];


};


  
  /* define a shuffle function. e.g. decl_shuffle(double).
 * advantage: compiler is free to optimize the swap operation without
 *            indirection with pointers, which could be much faster.
 * disadvantage: each datatype needs a separate instance of the function.
 *            for a small funciton like this, it's not very big a deal.
 */
#define decl_shuffle(type)				\
void shuffle_##type(type *list, size_t len) {		\
	int j;						\
	type tmp;					\
	while(len) {					\
		j = irand(len);				\
		if (j != len - 1) {			\
			tmp = list[j];			\
			list[j] = list[len - 1];	\
			list[len - 1] = tmp;		\
		}					\
		len--;					\
	}						\
}							\

/* random integer from 0 to n-1 */
int irand(int n)
{
	int r, rand_max = RAND_MAX - (RAND_MAX % n);
	/* reroll until r falls in a range that can be evenly
	 * distributed in n bins.  Unless n is comparable to
	 * to RAND_MAX, it's not *that* important really. */
	while ((r = rand()) >= rand_max);
	return r / (rand_max / n);
}
 
/* declare and define int type shuffle function from macro */
decl_shuffle(int);
  
  
int main()
{

  printf("%d\n", sizeof(M));
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
    
  current = 0;
  int m[SIZE][SIZE];
  long int record;
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
      d[i/SIZE][i%SIZE]+=dist;
    }
  }
  print(d, 0);
  
  //we initialize the matrix
  for( i=0; i<SIZE ;i++)
  {
    for( j=0; j<SIZE; j++)
    {
      m[i][j] = i*SIZE+j+1;
    }
  }

  //eliminate numbers with gcd of 1
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
          printf("#%d: %d = %d\n", u_numbers+1, m[i][j], k+1);
          m[i][j]=k+1;
          u_numbers++;
          break;
        }
      }   
    }
  }
  printf("\n%d redundant numbers were modified...\n", u_numbers);
  
  
  record=eval(m);
  record=eval(m);
  printf("\nAt the start, the record value is %ld.", record);
  init = clock();
  init2 = clock();
  totaltime = clock();
  print(m, record);
  find(m, &record, 0, 1,eval(m));
  print(m, record);
  long int r=record;
  while(1)
  {
    shuffle((int*)m, SIZE*SIZE);
	//print(m, record);
    //shuffle_int((int*)m, SIZE*SIZE);
	//print(m, record);
    find(m, &record, 0, 1,eval(m));
    if(record>r)
    {
      clock_t t = clock()-init2;
      init2 = clock();
      float ti  = ((float)t)/CLOCKS_PER_SEC;
      float total_t  = ((float)(clock()-totaltime))/CLOCKS_PER_SEC;
      total_t = total_t/60;
      r = record;
      print(m, record);
      printf ("It took me %.3f s, total: %.3f m (%d)\n", ti, total_t, SIZE);
    }
  }
  return 0;
}

void find(int m[SIZE][SIZE], long int *record, int index, int limit, long int total)
{
  long int swap1, swap2, newtotal, total_tmp, eval1;
  total_tmp=total;
  int best = 0;
  do{
    int i, j, a, b, cur;
    cur=total_tmp;
    a=0;
    b=0;
    best=0;
    for(i=index;i<(SIZE*SIZE-1);i++)
    {
      eval1 = eval_single(m, i/SIZE, i%SIZE);
      for(j=i+1;j<SIZE*SIZE;j++)
      {
        if(m[i/SIZE][i%SIZE]==m[j/SIZE][j%SIZE])// || fixed[i/SIZE][i%SIZE] || fixed[j/SIZE][j%SIZE])
          continue;
        //eval2 = eval_single(m, j/SIZE, j%SIZE);
        swap1 = eval1 + eval_single(m, j/SIZE, j%SIZE);
        //swap1 = eval_single(m, i/SIZE, i%SIZE) + eval_single(m, j/SIZE, j%SIZE);
        swap(i/SIZE, i%SIZE, j/SIZE, j%SIZE, m);
        swap2 = eval_single(m, i/SIZE, i%SIZE) + eval_single(m, j/SIZE, j%SIZE);
        newtotal = total_tmp - swap1 + swap2;
        if(newtotal>cur)
        {
          a = i;
          b = j;
          cur = newtotal;
          best=1;
        }
        swap(i/SIZE, i%SIZE, j/SIZE, j%SIZE, m);
      }
    }
    if(best)
    {
      swap(a/SIZE, a%SIZE, b/SIZE, b%SIZE, m);
      total_tmp=cur;
      *record=cur;
    }
    //clock_t t = clock()-init;
    //init = clock();
    //float ti  = ((float)t)/CLOCKS_PER_SEC;
    //float total_t  = ((float)(clock()-totaltime))/CLOCKS_PER_SEC;
    //total_t = total_t/60;
    //printf ("It took me %.3f s, total: %.3f m (%d)\n", ti, total_t, SIZE);
    //print(m, *record);
  }while(best!=0);
}

void print(int m[SIZE][SIZE], long int val)
{
  int i, j;
  printf("\n");
  for(i=0;i<SIZE;i++)
  {
    for(j=0;j<SIZE;j++)
    {
      if(j==0) printf("(");
      if(m[i][j]<100) printf(" ");
      if(m[i][j]<10) printf(" ");
      printf("%d", m[i][j]);
      if(j!=SIZE-1) printf(",");
    }
    if(i!=SIZE-1) printf("),");
    else  printf(")");
    printf("\n");
  }
  printf("\t val: %lu\n", val);
}

void copy(int source[SIZE][SIZE], int destiny[SIZE][SIZE])
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

void swap(int i, int j, int p, int q, int m[SIZE][SIZE])
{
  int tmp;
  tmp = m[i][j];
  m[i][j]=m[p][q];
  m[p][q]=tmp;
}

long int eval(int m[SIZE][SIZE])
{
  int i, j;
  long int res;
  //printf("a\tb\tgcd(a, b)\tDa,b\n");
  for(res=0, i=0;i<(SIZE*SIZE);i++)
  {
    for(j=(i+1);j<(SIZE*SIZE);j++)
    {
      res += gcds[m[i/SIZE][i%SIZE]-1][m[j/SIZE][j%SIZE]-1]*distances[i][j];
    }
  }
  return res;
}

long int eval_single(int m[SIZE][SIZE], int p, int q)
{
  int i, j, k;
  long int res;
  j = m[p][q]-1;
  k = p*SIZE+q;
  //printf("a\tb\tgcd(a, b)\tDa,b\n");
  for(res=0, i=0;i<(SIZE*SIZE);i++)
  {
    //res += gcds[m[p][q]-1][m[i/SIZE][i%SIZE]-1]*distances[p*SIZE+q][i];
    res += gcds[j][m[i/SIZE][i%SIZE]-1]*distances[k][i];
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
/* Arrange the N elements of ARRAY in random order.
   Only effective if N is much smaller than RAND_MAX;
   if this may not be the case, use a better random
   number generator. */
void shuffle(int *array, size_t n)
{
    if (n > 1) 
    {
        size_t i;
        for (i = 0; i < n - 1; i++) 
        {
          size_t j = i + rand() / (RAND_MAX / (n - i) + 1);
          int t = array[j];
          array[j] = array[i];
          array[i] = t;
        }
    }
}