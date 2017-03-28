#include "stdio.h"
int gcd(int m, int n);

int main()
{
  int n = 100;
  
  int i, j, cntr;
  
  for(i=1;i<=n;i++)
  {
    cntr=1;
    for(j=1;j<=n;j++)
    {
      if(gcd(i, j)!=1 && gcd(i, j)!=j)
        cntr++;
    }
    printf("%d\t%d\n", i, cntr);
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