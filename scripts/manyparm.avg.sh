yy=1948
while [ $yy -le 1976 ]
do
  for parm in PRES TMP UFLX VFLX
  do
    for mo in 01 02 03 04 05 06 07 08 09 10 11 12
    do
      fn=grb2d${yy}${mo}      
      wgrib $fn | grep $parm | wgrib -i $fn -nh -append -o ../${parm}.$yy
    done
  done
  yy=`expr $yy + 1`
done
