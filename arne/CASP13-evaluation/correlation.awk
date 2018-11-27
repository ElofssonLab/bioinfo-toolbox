

BEGIN {
  i=0 ;
    max=v2;
    if ( v1 > v2 )
      max=v1;
}

{
  if (NF >= max ) {
    i++ ;
    dat1[i] = $v1 ;
    dat2[i] = $v2 ;
    }}

END {

    count=i+1e-20

    for(i=1;i<=count;i++)
    {
	sum1 += dat1[i] ;
	sum2 += dat2[i] ;
	sumsq1 += dat1[i]*dat1[i] ;
	sumsq2 += dat2[i]*dat2[i] ;
    }
    
    ave1 = sum1 / count
    ave2 = sum2 / count
    mae = 0
    for(i=1;i<=count;i++)
    {
	diffsq1 += (dat1[i]-ave1)*(dat1[i]-ave1)
	diffsq2 += (dat2[i]-ave2)*(dat2[i]-ave2)
	dev += (dat1[i]-ave1)*(dat2[i]-ave2)
	mae +=sqrt((dat1[i]-dat2[i])*(dat1[i]-dat2[i]))
    }
    
    dev1 = (count*sumsq1-sum1*sum1)/(count*(count-1))
    dev2 = (count*sumsq2-sum2*sum2)/(count*(count-1))
    
    dev1 = sqrt(dev1)
    dev2 = sqrt(dev2)

    printf "%8d\t",count
    printf "%8.2f%8.2f\t", ave1, dev1
    printf "%8.2f%8.2f\t", ave2, dev2
    
    corr = sqrt(diffsq1) * sqrt(diffsq2) + 1.e-20
    corr = dev / corr
    printf "%f\t", corr
    printf "%f\n", mae / count
}

