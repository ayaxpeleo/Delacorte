#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <iostream>
#include <math.h>
#include <thread>
//#include <random>

#define SIZE 8
#define SMALLER 0
#define BIGGER 1
#define MODE BIGGER
#define CONC 8

void print(int m[SIZE][SIZE], long int val);
void print(float m[SIZE][SIZE], long int val);
inline void copy(int source[SIZE][SIZE], int destiny[SIZE][SIZE]);
inline void copy(float source[SIZE][SIZE], float destiny[SIZE][SIZE]);
int distance(int i, int j, int p, int q);
float fdistance(int i, int j, int p, int q);
int gcd(int m, int n);
int d[SIZE][SIZE];
float df[SIZE][SIZE];
float av[SIZE][SIZE];
float region[SIZE][SIZE];
clock_t init, init2;  
clock_t totaltime;
long int total, current;


int distances[SIZE*SIZE][SIZE*SIZE], gcds[SIZE*SIZE][SIZE*SIZE];
void shuffle(int *array, size_t n);
bool compare(long int i, long int j);

class M;
 
class Int
{
  friend M;
  public:
  int n, n2, pos;
  bool used;
  union
  {
    struct
    {
      unsigned int f1  :1;
      unsigned int f2  :1;
      unsigned int f3  :1;
      unsigned int f4  :1;
      unsigned int f5  :1;
      unsigned int f6  :1;
      unsigned int f7  :1;
      unsigned int f8  :1;
    };
    unsigned char flags;
  };
  void copy(Int *d);
};
void Int::copy(Int *d)
{
  d->n = n;
  d->n2 = n2;
  d->pos = pos;
  d->used = used;
  d->flags = flags;
}
class M
{
  public:
    M()
    {
      for(int i=0;i<SIZE*SIZE;i++)
      {
        m[i] = new Int;
        m[i]->pos = i;
        m[i]->n = i+1;
        m[i]->n2 = i+1;
        m[i]->flags = 0;
        lim = SIZE*SIZE;
        m[i]->used=false;
      }
      record=eval();
    }
    ~M()
    {
      for(int i=0;i<SIZE*SIZE;i++)
        delete m[i];
    }
    long int eval();
    void swap(int n1, int n2);
    void swap_pos(int n1, int n2);
    void reduceBigPrimes();
    void reduceRedundantNumbers();
    void print();
    long int eval_single(int p);
    long int move(int index, int limit, long int total, M *d, long int *best);
    long int find_move(int index, int limit, long int total, int *u, int *v);
    void shuffle();
    void shuffle_m();
    void shuffle2(int k);
    void ban(int n);
    void put(int n, int pos);
    bool isValidSwap(int i, int j);
    long int find(int index, int limit, long int total);
    void copy(M *d);
    void query();
  
    Int *m[SIZE*SIZE];
    int lim;
    long int record;
};

void M::query(void)
{
  int n[SIZE][SIZE];
  int i, j;
  int div;
  std::cout << "Enter divisor: ";
  std::cin >> div;
  if (div<=0) return;
  for(i=0;i<SIZE*SIZE;i++)
  {
    int pos = m[i]->pos;
    n[pos/SIZE][pos%SIZE]=m[i]->n2;
  }
  for(i=0;i<SIZE;i++)
  {
    for(j=0;j<SIZE;j++)
    {
      n[i][j]=gcd(n[i][j], div);
    }
  }
  ::print(n, div);
  return;
}

void M::copy(M *d)
{
  d->lim = lim;
  d->record = record;
  for(int i=0;i<SIZE*SIZE;i++)
  {
    m[i]->copy(d->m[i]);
  }
}

void M::put(int n, int pos)
{
  int i, j;
  //search for the element
  for(i=0;i<lim;i++)
  {
    if(m[i]->n2==n)
    {
      //found!!
      //search for the position
      for(j=0;j<lim;j++)
      {
        if(m[j]->pos==pos)
        {
          //found position, now let's swap them
          swap_pos(i, j);
          return;
        }
      }
    }
  }
}

void M::ban(int n)
{
  //check that we have an item within the limits
  if(n<0 || n>=SIZE*SIZE)
  {
    std::cout << "ERROR when trying to ban #" << n << "\n";
    return;
  }
  //find the number
  int i;
  for(i=0;i<lim;i++)
  {
    if(m[i]->n2 == n)
    {
      //send it to the tail
      swap(i, lim-1);
      //and reduce the amount of items to be checked every time
      lim--;
      return;
    }
  }
}

long int M::eval()
{
  int i, j;
  long int res=0;
  for(i=0;i<(SIZE*SIZE-1);i++)
  {
    for(j=i+1;j<(SIZE*SIZE);j++)
    {
      res += gcds[m[i]->n-1][m[j]->n-1]*distances[m[i]->pos][m[j]->pos];
    }
  }
  return res;
}
  
void M::swap(int n1, int n2)
{
  Int *tmp = m[n2];
  m[n2] = m[n1];
  m[n1] = tmp;
}
void M::swap_pos(int n1, int n2)
{
  int tmp = m[n2]->pos;
  m[n2]->pos = m[n1]->pos;
  m[n1]->pos = tmp;
}

void M::reduceBigPrimes()
{
  //eliminate numbers with gcd of 1
  for(int i=0;i<SIZE*SIZE;i++)
  {
    int flag_ok=0;
    for(int j=0;j<SIZE*SIZE;j++)
    {
      if(m[i]->n == j+1)
        continue;
      if(gcd(m[i]->n, j+1)!=1)
      {
        flag_ok=1;
      }
    }
    if(flag_ok==0)
    {
      std::cout << "Prime detected and marked: " << m[i]->n << "\n";
      m[i]->n=1;
    }
  }
}

void M::reduceRedundantNumbers()
{
  int i, k, l;
  //eliminate unnecesary numbers
  int u_numbers=0;
  for( i=0; i<SIZE*SIZE ;i++)
  {
    int check_number = m[i]->n-1;
    int flag_out=0;
    //if the number has a multiple in the array, we better leave
    
    for(k=1;k<=SIZE*SIZE;k++)
    {
      if(k%m[i]->n==0 && k!=m[i]->n)
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
        if(l==check_number)
          continue;
        if(gcds[check_number][l] != gcds[k][l])
        {
          flag_out=1;
        }
      }  
      if(flag_out==0 && m[i]->n>(k+1))
      {
        std::cout << "#" << u_numbers+1 << ": " << m[i]->n << " = " << k+1 << "\n";
        m[i]->n=k+1;
        u_numbers++;
        break;
      }
    }   
  }
  std::cout << "\n" << u_numbers << " redundant numbers were modified...\n";
}

void M::print()
{
  int h, i;
  std::cout << "\n";
  for(h=0;h<SIZE*SIZE;h++)
  {
    for(i=0;i<SIZE*SIZE;i++)
    {
      if(m[i]->pos==h)
      {
        int number=m[i]->n2;
          
        if((m[i]->pos%SIZE)==0) std::cout << "(";
        if(number<100) std::cout << " ";
        if(number<10) std::cout << " ";
        std::cout << number;
        if((m[i]->pos%SIZE)!=SIZE-1) std::cout << ",";
        if(m[i]->pos%SIZE==SIZE-1)
        {
          if(m[i]->pos/SIZE!=SIZE-1) std::cout << "),";
            else  std::cout << ")";
          std::cout << "\n";
        }
      }
    }
  }
  std::cout << "\t val: " << eval() << "\n" ;
}

long int M::eval_single(int p)
{
  int i, g1, p1;
  long int res;
  g1=m[p]->n-1;
  p1=m[p]->pos;
  for(res=0, i=0;i<SIZE*SIZE;i++)
  {
    res += gcds[g1][m[i]->n-1]*distances[p1][m[i]->pos];
  }
  return res;
}

long int M::move(int index, int limit, long int total, M *d, long int *best)
{
  if(index==limit)
    return total;
  long int res=0;
  long int eval[lim];
  for(int i=0;i<lim;i++) eval[i]=eval_single(i);
  
  for(int i=index; i<(lim-1) ;i++)
  {
    long int eval1 = eval[i];
    for(int j=i+1; j<lim; j++)
    {
      long int swap1, swap2, newtotal;
      if(m[i]->n==m[j]->n)
        continue;
      swap1 = eval1 + eval[j];
      swap_pos(i, j);
      swap2 = eval_single(i) + eval_single(j);
      newtotal = total - swap1 + swap2;
      
      long int tmpr = move(index+1, limit, newtotal, d, best);
      if(compare(tmpr, res))
      {
        res = tmpr;
        //copy(*d);
      }
      if(compare(newtotal, res))
      {
        res = newtotal;
      }
      if(compare(newtotal, *best))
      {
        *best = newtotal;
        copy(d);
        d->print();
      }
      
      swap_pos(i, j);
    }
  }
  return res;
}

long int M::find_move(int index, int limit, long int total, int *u, int *v)
{
  long int res, eval[lim], i;
  res=total;
  for(i=0;i<lim;i++) eval[i]=eval_single(i);
  
  for(i=0; i<(lim-1) ;i++)
  {
    for(int j=i+1; j<lim; j++)
    {
      long int swap1, swap2, newtotal;
      if(m[i]->n==m[j]->n)// || isValidSwap(i, j)==false)
        continue;
      swap1 = eval[i] + eval[j];
      swap_pos(i, j);
      swap2 = eval_single(i) + eval_single(j);
      newtotal = total - swap1 + swap2;
      
      //newr = move(index, limit, newtotal);
      if(compare(newtotal, res))
      {
        res = newtotal;
        *u=i;
        *v=j;
      }
      
      swap_pos(i, j);
    }
  }
  if(compare(res, total))
    return res;
  return 0;
}

long int M::find(int index, int limit, long int total)
{
  int i, j;
  i=0;
  j=0;
  long int ret, v1;
  ret=total;
  while(1)
  {
    v1 = ret;
    ret = find_move(index, limit, ret, &i, &j);
    if(ret==0)
      break;
    swap_pos(i, j);
  }
  return v1;
}

void M::shuffle()
{
  srand ( time(NULL) );
  if (lim > 1) 
  {
      int i;
      for (i = 0; i < lim - 1; i++) 
      {
        int j = i + rand() / (RAND_MAX / (lim - i) + 1);
        swap_pos(i, j);
      }
  }
}

void M::shuffle_m()
{
  srand ( time(NULL) );
  if (lim > 1) 
  {
      int i;
      for (i = 0; i < lim - 1; i++)
      {
        int j = i + rand() / (RAND_MAX / (lim - i) + 1);
        swap(i, j);
      }
  }
}
void M::shuffle2(int k)
{
  //shuffle_m();
  int i, j;
  for (i = lim-1; i > (lim-k) && i>=0; i--) 
  {
    j = rand() % (i+1);
    swap_pos(i, j);
  }
}

bool M::isValidSwap(int i, int j)
{
  float val1, val2;
  val1 = fabs(region[(m[i]->n2-1)/SIZE][(m[i]->n2-1)%SIZE] - df[m[j]->pos/SIZE][m[j]->pos%SIZE])/region[(m[i]->n2-1)/SIZE][(m[i]->n2-1)%SIZE];
  val2 = fabs(region[(m[j]->n2-1)/SIZE][(m[j]->n2-1)%SIZE] - df[m[i]->pos/SIZE][m[i]->pos%SIZE])/region[(m[j]->n2-1)/SIZE][(m[j]->n2-1)%SIZE];
  if(0 && record>=57134)
  {
    std::cout << "val1: " << val1 << "\n";
    std::cout << "val2: " << val2 << "\n";
  }
  return true;
  float l = 0.5;
  
  if(val1>l && val2>l)
    return false;
  else
    return true;
}

int main()
{

  //srand ( time(NULL) );
  srand ( 1792917929 );
  M mat;
  mat.print();
  int i, j;
  /*
  int n[SIZE][SIZE] ={{ 18, 28, 33, 10, 16, 24},
                      { 35,  4, 13,  2, 34, 15},
                      {  9,  5,  1, 19,  7,  6},
                      {  8,  3, 23, 29, 11, 27},
                      { 12, 17, 31, 25, 26, 20},
                      { 30, 32, 22, 21, 14, 36}};
                        
  for(i=0;i<SIZE*SIZE;i++)
  {
    mat.m[i]->n = n[i/SIZE][i%SIZE];
    mat.m[i]->n2 = n[i/SIZE][i%SIZE];
    mat.m[i]->pos = i;
  }
  */
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
  float smallest;
  for(i=0;i<SIZE*SIZE;i++)
  {
    d[i/SIZE][i%SIZE]=0;
  }
  for(i=0;i<SIZE*SIZE;i++)
  {
    for(j=0;j<SIZE*SIZE;j++)
    {
      float dist = fdistance(i/SIZE, i%SIZE, j/SIZE, j%SIZE);
      df[i/SIZE][i%SIZE]+=dist;
    }
  }
  std::cout << "\n  Distances matrix: \n";
  print(df, 0);
  smallest=1000000000;
  for(i=0;i<SIZE*SIZE;i++)
  {
    if(df[i/SIZE][i%SIZE]<smallest)
      smallest = df[i/SIZE][i%SIZE];
  }
  for(i=0;i<SIZE*SIZE;i++)
  {
    df[i/SIZE][i%SIZE] = ((float)df[i/SIZE][i%SIZE])/((float)smallest);
  }
  std::cout << "\n  Normalized distance matrix: \n";
  print(df, 0);
  
  
  float av_tmp[SIZE][SIZE], df_tmp[SIZE][SIZE];
  /*
  for(i=0;i<SIZE*SIZE;i++) df_tmp[i/SIZE][i%SIZE]=0;
  for(i=0;i<SIZE;i++)
  {
    for(j=0;j<SIZE;j++)
    {
      int cntr=0;
      for(int u=-1; u<2 ; u++)
      {
        for(int v=-1; v<2; v++)
        {
          if((i+u)<0 || (i+u)>=SIZE || (j+v)<0 || (j+v)>=SIZE)
            continue;
          df_tmp[i][j]+=df[i+u][j+v];
          cntr++;
        }
      }
      if(cntr!=0)
        df_tmp[i][j] = df_tmp[i][j]/cntr;
      else std::cout << "\nERROR dividing by zero.\n";
    }
  }
  copy(df_tmp, df);
  
  std::cout << "\n  Normalized and smoothed distance matrix: \n";
  print(df, 0);
  */
  for(i=0;i<SIZE*SIZE;i++)
  {
    av[i/SIZE][i%SIZE] = 0;
  }
  for(i=0;i<SIZE*SIZE;i++)
  {
    for(j=0;j<SIZE*SIZE;j++)
    {
      if(i==j) continue;
      av[i/SIZE][i%SIZE] += ((float)gcds[i][j])/((float)(SIZE*SIZE-1));
    }
  }
  std::cout << "\n  Normalized average matrix: \n";
  print(av, 0);
  
  //initialize array;
  for(i=0;i<SIZE;i++)
  {
    for(j=0;j<SIZE;j++)
    {
      region[i][j] = 0;
    }
  }
  
  
  copy(av, av_tmp);
  copy(df, df_tmp);
  for(i=0;i<SIZE*SIZE;i++)
  {
    int big_d, big_av=0;
    //find the biggest value's index in the average
    for(j=0, big_av=0;j<SIZE*SIZE;j++)
    {
      if(av_tmp[j/SIZE][j%SIZE]>av_tmp[big_av/SIZE][big_av%SIZE])
        big_av = j;
    }
    //find the biggest value's index in the distances
    for(j=0, big_d=0;j<SIZE*SIZE;j++)
    {
      if(df_tmp[j/SIZE][j%SIZE]>df_tmp[big_d/SIZE][big_d%SIZE])
        big_d = j;
    }
    //let's arrange the first matrix
    mat.m[i]->n=big_av+1;
    mat.m[i]->n2=big_av+1;
    mat.m[i]->pos=big_d;
    //put the biggest value's index in the place of the biggest distance
    region[big_av/SIZE][big_av%SIZE] = df_tmp[big_d/SIZE][big_d%SIZE];
    //set both the used biggest value and distance to zero, to mark the places as used
    av_tmp[big_av/SIZE][big_av%SIZE]=0;
    df_tmp[big_d/SIZE][big_d%SIZE]=0;
  }
  mat.record = mat.eval();
  std::cout << "\n  Region matrix: \n";
  print(region, 0);
  
  std::cout << "\n  Rearranged matrix: \n";  
  mat.print();
  std::cout << "\n" << mat.eval() << "\n";
  std::cout << "\nAt the start, the record value is " << mat.eval();
  init = clock();
  init2 = clock();
  totaltime = clock();
  mat.print();
  mat.reduceBigPrimes();
  mat.reduceRedundantNumbers();
  
  long int r=mat.eval();
  mat.record=mat.eval();
  int first_flag=1;
  //hola
  M tmp;
  /*
  M tmpa[CONC];
  for(i=0; i<CONC; i++)
  {
    mat.copy(&tmpa[i]);
  }
  */
  //void M::put(int n, int pos)
  //void M::ban(int n)
  mat.put(30, 0);
  mat.ban(30);
  
  mat.put(42, 7);
  mat.ban(42);
  
  mat.put(48, 56);
  mat.ban(48);
  
  mat.put(60, 63);
  mat.ban(60);
  
  
  
  
  mat.copy(&tmp);
  long int record=137591;
  long int cycle0=0;
  int disturbance_size;
  for(disturbance_size=SIZE*SIZE-1-4;disturbance_size>(SIZE*SIZE)/2;disturbance_size-=1)
  {
    printf("Now trying with disturbance size = %d...\n", disturbance_size);
    for(cycle0=0;cycle0<100000;cycle0++)
    {
      tmp.copy(&mat);
      r = mat.find(0, 1, mat.eval());
      int go=0;
      if(compare(r, record))
        go=1;
      //std::cout << "matrix to be analized:'n";
      //mat.print();
      //mat.query();
      //std::cout << "m1." << r << ", " << mat.record << "\n";
      if(go)
      {
        if(first_flag)
        {
          std::cout << "FIRST RESULT:\n";
          first_flag=0;
        }
        clock_t t = clock()-init2;
        init2 = clock();
        float ti  = ((float)t)/CLOCKS_PER_SEC;
        float total_t  = ((float)(clock()-totaltime))/CLOCKS_PER_SEC;
        total_t = total_t/60;
        record=r;
        mat.print();
        printf ("It took me %.3f s, total: %.3f m (%d)\n", ti, total_t, SIZE);
        //std::cout << "mat.record: " << mat.record << "\n";
        //if(total_t>0.5)
        //  return 0;
      }
      //break;
      tmp.shuffle2(disturbance_size);
      //tmp.shuffle();
    }
  }
  return 0;
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

void copy(float source[SIZE][SIZE], float destiny[SIZE][SIZE])
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
float fdistance(int i, int j, int p, int q)
{
  int tmp1, tmp2;
  tmp1 = (p-i)*(p-i);
  tmp2 = (q-j)*(q-j);
  return (float)sqrt(tmp1+tmp2);
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
void print(int m[SIZE][SIZE], long int val)
{
  int i, j;
  std::cout << "\n";
  for(i=0;i<SIZE;i++)
  {
    for(j=0;j<SIZE;j++)
    {
      if(j==0) std::cout << "(";
      if(m[i][j]<100) std::cout << " ";
      if(m[i][j]<10) std::cout << " ";
      std::cout << m[i][j];
      if(j!=SIZE-1) std::cout << ",";
    }
    if(i!=SIZE-1) std::cout << "),";
    else std::cout << ")";
    std::cout << "\n";
  }
   std::cout << "\t val: " << val << "\n";
}
void print(float m[SIZE][SIZE], long int val)
{
  int i, j;
  std::cout << "\n";
  std::cout.precision(2);
  std::cout.setf( std::ios::fixed, std:: ios::floatfield );
  for(i=0;i<SIZE;i++)
  {
    for(j=0;j<SIZE;j++)
    {
      if(j==0) std::cout << "(";
      if(m[i][j]<100) std::cout << " ";
      if(m[i][j]<10) std::cout << " ";
      std::cout << m[i][j];
      if(j!=SIZE-1) std::cout << ",";
    }
    if(i!=SIZE-1) std::cout << "),";
    else std::cout << ")";
    std::cout << "\n";
  }
   std::cout << "\t val: " << val << "\n";
}

bool compare(long int i, long int j)
{
  if(MODE==SMALLER)
  {
    if(i<j)
      return true;
    else
      return false;
  }
  else
  {
    if(i>j)
      return true;
    else
      return false;
  }
}